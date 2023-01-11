//
//  rpfbwt_algorithm.hpp
//
//  Copyright 2022 Marco Oliva. All rights reserved.
//

#ifndef rpfbwt_algorithm_hpp
#define rpfbwt_algorithm_hpp

#include <filesystem>
#include <sys/stat.h>
#include <omp.h>

#include <pfp/pfp.hpp>

#undef max
#include <rle/rle_string.hpp>

namespace rpfbwt
{

template <typename dict_l1_data_type, typename parse_int_type = uint32_t>
class rpfbwt_algo
{
public:
    
    class l2_colex_comp
    {
    private:
        
        const pfpds::dictionary<dict_l1_data_type>& l1_d;
        parse_int_type int_shift;
    
    public:
        
        l2_colex_comp(const pfpds::dictionary<dict_l1_data_type>& l1_d_ref, parse_int_type shift)
        : l1_d(l1_d_ref) , int_shift(shift) {}
        
        bool operator()(parse_int_type l, parse_int_type r)
        {
            assert(l != r);
            if (l < int_shift) { return true; }
            if (r < int_shift) { return false; }
            return l1_d.colex_id[l - int_shift - 1] < l1_d.colex_id[r - int_shift - 1];
        }
    };
    
private:
    
    const std::size_t int_shift = 10;
    std::string l1_prefix;
    std::string out_rle_name;
    
    std::less<dict_l1_data_type> l1_d_comp;
    pfpds::dictionary<dict_l1_data_type> l1_d;
    std::vector<std::size_t> l1_d_lengths;
    std::vector<uint_t> l1_freq; // here occ has the same size as the integers used for gsacak.
    bool l1_cleared = false;
    
    l2_colex_comp l2_comp;
    pfpds::pf_parsing<parse_int_type, l2_colex_comp, pfpds::pfp_wt_sdsl> l2_pfp;
    std::vector<std::vector<std::size_t>> E_arrays;
    std::size_t l1_n = 0; // also used to initialize l1_n
    bool l2_cleared = false;
    
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> l2_pfp_v_table; // (row in which that char appears, number of times per row)

    static constexpr std::size_t default_num_of_chunks = 1;
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> chunks;
    
    rle::RLEString::RLEncoderMerger rle_chunks;
    std::ofstream sa_values_ofstream;
    std::ofstream sa_positions_ofstream;
    
    
    void init_v_table()
    {
        spdlog::info("Creating v table");
        
        // now build V
        for (std::size_t row = 0; row < l2_pfp.M.size(); row++)
        {
            const auto& m =  l2_pfp.M[row];
        
            std::vector<std::pair<parse_int_type, std::size_t>> phrase_counts;
            for (std::size_t r = m.left; r <= m.right; r++)
            {
                auto phrase = l2_pfp.dict.colex_id[r];
                std::size_t phrase_start = l2_pfp.dict.select_b_d(phrase + 1);
                std::size_t phrase_length = l2_pfp.dict.length_of_phrase(phrase + 1);
                parse_int_type c = l2_pfp.dict.d[phrase_start + (phrase_length - m.len - 1)];
            
                if (phrase_counts.empty() or phrase_counts.back().first != c)
                {
                    phrase_counts.push_back(std::make_pair(c, l2_pfp.freq[phrase + 1]));
                }
                else
                {
                    phrase_counts.back().second += l2_pfp.freq[phrase + 1];// 1;
                }
            }
        
            // update
            for (auto const& c_count : phrase_counts)
            {
                std::pair<std::size_t, std::size_t> v_element = std::make_pair(row, c_count.second);
                l2_pfp_v_table[c_count.first].emplace_back(v_element);
            }
        }
    }
    
    std::size_t init_E_arrays()
    {
        spdlog::info("Creating E arrays");
        
        E_arrays.resize(l2_pfp.dict.n_phrases());
        std::size_t end_pos_from_start = 0;
        
        // read P2 and expand each phrase to get the end position in the text order
        std::vector<std::size_t> end_positions;
        for (std::size_t i = 0; i < l2_pfp.pars.p.size() - 1; i++)
        {
            std::size_t l2_pid = l2_pfp.pars.p[i];
            std::size_t l2_phrase_start = l2_pfp.dict.select_b_d(l2_pid);
            if (l2_phrase_start < l2_pfp.w) { l2_phrase_start = l2_pfp.w; }
            std::size_t l2_phrase_end = l2_pfp.dict.select_b_d(l2_pid + 1) - 2;
            
            for (std::size_t j = l2_phrase_start; j <= (l2_phrase_end - l2_pfp.w); j++)
            {
                std::size_t l1_pid = l2_pfp.dict.d[j] - int_shift;
                std::size_t l1_phrase_start = l1_d.select_b_d(l1_pid);
                std::size_t l1_phrase_end = l1_d.select_b_d(l1_pid + 1) - 2;
                std::size_t l1_phrase_length = (l1_phrase_end - l1_phrase_start) + 1;
    
                end_pos_from_start += (l1_phrase_length - l1_d.w);
            }
            end_positions.push_back(end_pos_from_start);
        }
        
        // now reorder according to P2's SA
        for (std::size_t i = 0; i < l2_pfp.pars.saP.size(); i++)
        {
            std::size_t sa_value = l2_pfp.pars.saP[i];
            
            if (sa_value == 0) { continue; } // { E_arrays[l2_pfp.pars.p[l2_pfp.pars.p.size() - 2]].push_back(end_positions.back()); }
            else
            {
                E_arrays[l2_pfp.pars.p[sa_value - 1] - 1].push_back(end_positions[sa_value - 1]);
            }
        }
        
        return end_pos_from_start;
    }
    
    // Computes ranges for parallel computation
    // suffix start, suffix end, this_left, this_row
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>
    compute_chunks(std::size_t num_of_chunks)
    {
        spdlog::info("Computing chunks for parallel execution. Total input size: {}. Requested chunks: {}", l1_n, num_of_chunks);
        
        std::size_t chunk_size = (num_of_chunks > 1) ? (l1_n / (num_of_chunks - 1)) : l1_n + 1;
        
        std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> out;
        
        // Go through the suffixes of D and compute chunks
        std::size_t i = 1; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
        
        std::size_t l_left  = 0;
        std::size_t l_right = 0;
        std::size_t chunk_suffix_start = i;
        std::size_t chunk_start = l_left;
        std::size_t chunk_row_start = 0;
        std::size_t table_row = 0;
        while (i < l1_d.saD.size())
        {
            auto sn = l1_d.saD[i];
            // Check if the suffix has length at least w and is not the complete phrase.
            auto phrase = l1_d.daD[i] + 1;
            assert(phrase > 0 && phrase < l1_freq.size()); // + 1 because daD is 0-based
            uint_t suffix_length = l1_d.select_b_d(l1_d.rank_b_d(sn + 1) + 1) - sn - 1;
            if (l1_d.b_d[sn] || suffix_length < l1_d.w) // skip full phrases or suffixes shorter than w
            {
                ++i; // Skip
            }
            else
            {
                i++;
                l_right += l1_freq[phrase] - 1;
                
                if (i < l1_d.saD.size())
                {
                    auto new_sn = l1_d.saD[i];
                    auto new_phrase = l1_d.daD[i] + 1;
                    assert(new_phrase > 0 && new_phrase < l1_freq.size()); // + 1 because daD is 0-based
                    uint_t new_suffix_length = l1_d.select_b_d(l1_d.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
                    
                    while (i < l1_d.saD.size() && (l1_d.lcpD[i] >= suffix_length) && (suffix_length == new_suffix_length))
                    {
                        ++i;
                        l_right += l1_freq[new_phrase];
                        
                        if (i < l1_d.saD.size())
                        {
                            new_sn = l1_d.saD[i];
                            new_phrase = l1_d.daD[i] + 1;
                            assert(new_phrase > 0 && new_phrase < l1_freq.size()); // + 1 because daD is 0-based
                            new_suffix_length = l1_d.select_b_d(l1_d.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
                        }
                    }
                    
                }
                
                if (l_right - chunk_start > chunk_size)
                {
                    // store SA_d range for this chunk
                    out.emplace_back(chunk_suffix_start, i, chunk_start, chunk_row_start);
                    
                    // prepare for next iteration
                    chunk_suffix_start = i;
                    chunk_start = l_right + 1;
                    chunk_row_start = table_row + 1;
                }
                
                l_left = l_right + 1;
                l_right = l_left;
                table_row++;
            }
        }
        
        // add last_chunk
        if ((out.empty()) or (std::get<1>(out.back()) < i))
        {
            // store SA_d range for last chunk
            out.emplace_back(chunk_suffix_start, i, chunk_start, chunk_row_start);
        }
        
        spdlog::info("BWT ranges divided in {} chunks.", out.size());
        if (out.size() > 500)
        {
            spdlog::error("Each chunk creates a temporary file on disk. To now overwhelm the filesystem increase the chunk size.");
            std::exit(1);
        }
        return out;
    }
    
    std::vector<uint_t> read_l1_freq(std::string& l1_prefix, std::size_t l1_d_phrases)
    {
        spdlog::info("Reading in l1 frequencies");
        
        // get l1 parse size to get appropriate int size for occurrences
    
        std::filesystem::path parse_path(std::string(l1_prefix + ".parse"));
        std::size_t parse_size = std::filesystem::file_size(parse_path) / sizeof(uint32_t);

        spdlog::info("l1 parse size: {} bytes", parse_size * sizeof(uint32_t));
    
        std::size_t occ_bytes;
        if (parse_size < (std::numeric_limits<uint32_t>::max() - 1))
        { spdlog::info("Reading in 32 bits occ file"); occ_bytes = 4; }
        else
        { spdlog::info("Reading in 64 bits occ file"); occ_bytes = 8; }
    
        std::vector<uint_t> out;
        out.push_back(0);
        std::ifstream occ_file(std::string(l1_prefix + ".occ"), std::ios::binary);
        for (std::size_t i = 0; i < l1_d_phrases; i++)
        {
            uint_t occ_value;
            occ_file.read((char*)&occ_value, occ_bytes);
            out.push_back(occ_value);
        }
        
        return out;
    }
    
    std::vector<std::size_t> init_d1_lengths(pfpds::dictionary<dict_l1_data_type> l1_d)
    {
        spdlog::info("Precomputing l1 d lengths");
        
        std::vector<std::size_t> out;
        
        std::size_t curr_length = 0;
        for (std::size_t i = 0; i < l1_d.d.size(); i++)
        {
            curr_length += 1;
            if (l1_d.d[i] == EndOfWord) { out.push_back(curr_length - 1); curr_length = 0; }
        }
        
        return out;
    }
    
    
    bool clear_L1_unused_data_structures()
    {
        spdlog::info("Removing unused L1 data structures.");
        l1_d.inv_colex_id.resize(0);
        l1_d.colex_daD.resize(0);
        l1_d.isaD.resize(0);
        return true;
    }
    
    bool clear_L2_unused_data_structures()
    {
        spdlog::info("Removing unused L2 data structures.");
        l2_pfp.dict.inv_colex_id.resize(0);
        l2_pfp.dict.colex_daD.resize(0);
        l2_pfp.dict.lcpD.resize(0);
        l2_pfp.dict.daD.resize(0);
        l2_pfp.dict.isaD.resize(0);
        l2_pfp.pars.isaP.resize(0);
        l2_pfp.pars.saP.resize(0);
        l2_pfp.pars.lcpP.resize(0);
        return true;
    }
    
public:
    
    rpfbwt_algo(std::vector<dict_l1_data_type>& l1_d_v,
                std::vector<uint_t>& l1_freq_v,
                std::size_t l1_w,
                std::vector<parse_int_type>& l2_d_v,
                std::vector<parse_int_type>& l2_p_v,
                std::vector<uint_t>& l2_freq_v,
                std::size_t l2_w,
                std::size_t bwt_chunks = default_num_of_chunks)
    : l1_d(l1_d_v, l1_w, l1_d_comp, true, false, true, true, false, true, false), l1_freq(l1_freq_v),
      l1_cleared(clear_L1_unused_data_structures()),
      l1_d_lengths(init_d1_lengths(l1_d)),
      l1_prefix("mem"), out_rle_name(l1_prefix + ".rlebwt"),
      l2_comp(l1_d, int_shift),
      l2_pfp(l2_d_v, l2_comp, l2_p_v, l2_freq_v, l2_w, int_shift, false, true),
      l2_pfp_v_table(l2_pfp.dict.alphabet_size),
      l1_n(init_E_arrays()),
      l2_cleared(clear_L2_unused_data_structures()),
      chunks(compute_chunks(bwt_chunks)),
      rle_chunks(out_rle_name, chunks.size()),
      sa_values_ofstream(l1_prefix + ".ssa", std::ios::out | std::ios::binary),
      sa_positions_ofstream(l1_prefix + ".pssa", std::ios::out | std::ios::binary)
    { init_v_table(); }
    
    rpfbwt_algo(std::string& l1_prefix, std::size_t l1_w, std::size_t l2_w,  std::size_t bwt_chunks = default_num_of_chunks)
    : l1_d(l1_prefix, l1_w, l1_d_comp, true, false, true, true, false, true, false),
      l1_cleared(clear_L1_unused_data_structures()),
      l1_d_lengths(init_d1_lengths(l1_d)),
      l1_prefix(l1_prefix), out_rle_name(l1_prefix + ".rlebwt"),
      l1_freq(read_l1_freq(l1_prefix, l1_d.n_phrases())),
      l2_comp(l1_d, int_shift),
      l2_pfp(l1_prefix + ".parse", l2_w, l2_comp, int_shift, false, true),
      l2_pfp_v_table(l2_pfp.dict.alphabet_size),
      l1_n(init_E_arrays()),
      l2_cleared(clear_L2_unused_data_structures()),
      chunks(compute_chunks(bwt_chunks)),
      rle_chunks(out_rle_name, chunks.size()),
      sa_values_ofstream(l1_prefix + ".ssa", std::ios::out | std::ios::binary),
      sa_positions_ofstream(l1_prefix + ".pssa", std::ios::out | std::ios::binary)
    { init_v_table(); }
    
    void
    l1_bwt_chunk(
    const std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>& chunk,
    rle::RLEString::RLEncoder& rle_out,
    std::vector<std::size_t>& runs_heads_positions)
    {
        std::size_t i = std::get<0>(chunk);
        
        dict_l1_data_type prev_char = 0;
    
        std::size_t l_left  = std::get<2>(chunk);
        std::size_t l_right = l_left;
        std::size_t easy_chars = 0;
        std::size_t hard_easy_chars = 0;
        std::size_t hard_hard_chars = 0;
        std::size_t row = std::get<3>(chunk);
        while (i < std::get<1>(chunk))
        {
            auto sn = l1_d.saD[i];
            // Check if the suffix has length at least w and is not the complete phrase.
            auto phrase = l1_d.daD[i] + 1;
            assert(phrase > 0 && phrase < l1_freq.size()); // + 1 because daD is 0-based
            uint_t suffix_length = l1_d.select_b_d(l1_d.rank_b_d(sn + 1) + 1) - sn - 1;
            if (l1_d.b_d[sn] || suffix_length < l1_d.w) // skip full phrases or suffixes shorter than w
            {
                ++i; // Skip
            }
            else
            {
                std::set<dict_l1_data_type> chars;
                chars.insert(l1_d.d[((sn + l1_d.d.size() - 1) % l1_d.d.size())]);
    
                std::set<parse_int_type> pids;
                pids.insert(phrase);
                
                i++;
                l_right += l1_freq[phrase] - 1;
            
                if (i < l1_d.saD.size())
                {
                    auto new_sn = l1_d.saD[i];
                    auto new_phrase = l1_d.daD[i] + 1;
                    assert(new_phrase > 0 && new_phrase < l1_freq.size()); // + 1 because daD is 0-based
                    uint_t new_suffix_length = l1_d.select_b_d(l1_d.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
                
                    while (i < l1_d.saD.size() && (l1_d.lcpD[i] >= suffix_length) && (suffix_length == new_suffix_length))
                    {
                        chars.insert(l1_d.d[((new_sn + l1_d.d.size() - 1) % l1_d.d.size())]);
                        pids.insert(new_phrase);
                        ++i;
                    
                        l_right += l1_freq[new_phrase];
                    
                        if (i < l1_d.saD.size())
                        {
                            new_sn = l1_d.saD[i];
                            new_phrase = l1_d.daD[i] + 1;
                            assert(new_phrase > 0 && new_phrase < l1_freq.size()); // + 1 because daD is 0-based
                            new_suffix_length = l1_d.select_b_d(l1_d.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
                        }
                    }
                
                }
                
                if (chars.size() == 1) // easy-easy suffixes
                {
                    // easy suffixes
                    rle_out(*(chars.begin()), (l_right - l_left) + 1);
                    
                    if (prev_char == 0 or prev_char != *(chars.begin()))
                    {
                        std::size_t run_start_pos = easy_chars + hard_easy_chars + hard_hard_chars;
                        runs_heads_positions.push_back(easy_chars + hard_easy_chars + hard_hard_chars + std::get<2>(chunk));
                        prev_char = *(chars.begin());
                    }
                    
                    easy_chars += (l_right - l_left) + 1;
                }
                else // easy-hard and hard-hard suffixes
                {
                    assert(l_right - l_left > 0);
                    std::size_t hard_chars_before = hard_easy_chars + hard_hard_chars;
                    
                    // shorthands
                    typedef std::reference_wrapper<std::vector<std::pair<std::size_t, std::size_t>>> ve_t; // (row in which that char appears, number of times per row)
                    typedef std::size_t ve_pq_t; // for the priority queue the row is enough, the other ingo can be retrieved from ve_t
                    typedef std::pair<ve_pq_t, std::pair<std::size_t, std::size_t>> pq_t; // .first : row, this element is the .second.second-th element of the .second.first array
                    
                    // go in second level and iterate through list of positions for the meta-phrases we are interested into
                    std::vector<ve_t> v;
                    std::vector<parse_int_type> pids_v;
                    for (const auto& pid : pids)
                    {
                        parse_int_type adj_pid = pid + int_shift;
                        v.push_back(std::ref(l2_pfp_v_table[adj_pid]));
                        pids_v.push_back(adj_pid);
                    }

                    // make a priority queue and add elements to it
                    std::priority_queue<pq_t, std::vector<pq_t>, std::greater<>> pq;
                    for (std::size_t vi = 0; vi < v.size(); vi++) { pq.push({ v[vi].get()[0].first, { vi, 0 } }); }

                    while (not pq.empty())
                    {
                        // get all chars from this row entry at l2
                        std::vector<pq_t> from_same_l2_suffix;
                        auto first_suffix = pq.top();
                        while ((not pq.empty()) and (pq.top().first == first_suffix.first))
                        {
                            auto curr = pq.top(); pq.pop();
                            from_same_l2_suffix.push_back(curr);

                            std::size_t arr_i_c = curr.second.first;  // ith array
                            std::size_t arr_x_c = curr.second.second; // index in i-th array
                            if (arr_x_c + 1 < v[arr_i_c].get().size())
                            {
                                pq.push({ v[arr_i_c].get()[arr_x_c + 1].first, { arr_i_c, arr_x_c + 1 } });
                            }
                        }

                        if (from_same_l2_suffix.empty())
                        {
                            spdlog::error("Something went wrong.");
                        }
                        else if (from_same_l2_suffix.size() == 1)
                        {
                            // hard-easy suffix, we can fill in the character in pid
                            auto pq_entry = from_same_l2_suffix[0];
                            parse_int_type adj_pid = pids_v[pq_entry.second.first];
                            uint_t freq = v[pq_entry.second.first].get()[pq_entry.second.second].second;
                            dict_l1_data_type c = l1_d.d[l1_d.select_b_d(adj_pid - int_shift + 1) - (suffix_length + 2)]; // end of next phrases
                            
                            rle_out(c, freq);
                            
                            if (prev_char == 0 or prev_char != c)
                            {
                                std::size_t run_start_pos = easy_chars + hard_easy_chars + hard_hard_chars;
                                runs_heads_positions.push_back(easy_chars + hard_easy_chars + hard_hard_chars + std::get<2>(chunk));
                                prev_char = c;
                            }
                            
                            hard_easy_chars += freq;
                        }
                        else
                        {
                            // hard-hard suffix, need to look at the grid in L2
                            std::size_t curr_l2_row = from_same_l2_suffix[0].first;
                            auto l2_M_entry = l2_pfp.M[curr_l2_row];
                            
                            // get inverted lists of corresponding phrases
                            std::vector<std::reference_wrapper<std::vector<uint_t>>> ilists;
                            std::vector<dict_l1_data_type> ilist_corresponding_chars;
                            for (std::size_t c_it = l2_M_entry.left; c_it <= l2_M_entry.right; c_it ++)
                            {
                                // check if we need that pid
                                std::size_t l2_pid = l2_pfp.dict.colex_id[c_it] + 1;
                                parse_int_type adj_l1_pid = l2_pfp.dict.d[l2_pfp.dict.select_b_d(l2_pid + 1) - (l2_M_entry.len + 2)];
                                assert(adj_l1_pid >= int_shift);
                                parse_int_type l1_pid = adj_l1_pid - int_shift;
    
                                if (pids.contains(l1_pid))
                                {
                                    ilists.push_back(std::ref(l2_pfp.bwt_p_ilist[l2_pfp.dict.colex_id[c_it] + 1]));
                                    dict_l1_data_type c = l1_d.d[l1_d.select_b_d(l1_pid + 1) - (suffix_length + 2)];
                                    ilist_corresponding_chars.push_back(c);
                                }
                            }
    
                            // make a priority queue from the inverted lists
                            typedef std::pair<uint_t, std::pair<std::size_t, std::size_t>> ilist_pq_t;
                            std::priority_queue<ilist_pq_t, std::vector<ilist_pq_t>, std::greater<>> ilist_pq;
                            for (std::size_t il_i = 0; il_i < ilists.size(); il_i++) { ilist_pq.push({ ilists[il_i].get()[0], { il_i, 0 } }); }
                            
                            // now pop elements from the priority queue and write out the corresponding character
                            while (not ilist_pq.empty())
                            {
                                auto curr = ilist_pq.top(); ilist_pq.pop();
                                
                                // output corresponding char
                                dict_l1_data_type c = ilist_corresponding_chars[curr.second.first];

                                rle_out(c, 1);
    
                                if (prev_char == 0 or prev_char != c)
                                {
                                    std::size_t run_start_pos = easy_chars + hard_easy_chars + hard_hard_chars;
                                    runs_heads_positions.push_back(easy_chars + hard_easy_chars + hard_hard_chars + std::get<2>(chunk));
                                    prev_char = c;
                                }
                                
                                hard_hard_chars += 1;
                                
                                // keep iterating
                                std::size_t arr_i_c_il = curr.second.first;  // ith array
                                std::size_t arr_x_c_il = curr.second.second; // index in i-th array
                                if (arr_x_c_il + 1 < ilists[arr_i_c_il].get().size())
                                {
                                    ilist_pq.push({ ilists[arr_i_c_il].get()[arr_x_c_il + 1], { arr_i_c_il, arr_x_c_il + 1 } });
                                }
                            }
                        }
                    }
                    
                    // check that we covered the range we were working on
                    assert( ((hard_easy_chars + hard_hard_chars) - hard_chars_before) == ((l_right - l_left) + 1) );
                }
                
                l_left = l_right + 1;
                l_right = l_left;
                row++;
            }
        }
    }
    
    //------------------------------------------------------------
    
    void
    l1_sa_values_chunk(
    const std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>& chunk,
    const std::vector<std::size_t>& runs_heads_positions,
    std::vector<std::size_t>& out_sa_values)
    {
        std::size_t i = std::get<0>(chunk);
    
        std::size_t sa_pos_iterator = 0;
        
        std::size_t l_left  = std::get<2>(chunk);
        std::size_t l_right = l_left;
        std::size_t easy_chars = 0;
        std::size_t hard_easy_chars = 0;
        std::size_t hard_hard_chars = 0;
        std::size_t row = std::get<3>(chunk);
        while (i < std::get<1>(chunk))
        {
            auto sn = l1_d.saD[i];
            // Check if the suffix has length at least w and is not the complete phrase.
            auto phrase = l1_d.daD[i] + 1;
            assert(phrase > 0 && phrase < l1_freq.size()); // + 1 because daD is 0-based
            uint_t suffix_length = l1_d.select_b_d(l1_d.rank_b_d(sn + 1) + 1) - sn - 1;
            if (l1_d.b_d[sn] || suffix_length < l1_d.w) // skip full phrases or suffixes shorter than w
            {
                ++i; // Skip
            }
            else
            {
                std::set<dict_l1_data_type> chars;
                chars.insert(l1_d.d[((sn + l1_d.d.size() - 1) % l1_d.d.size())]);
            
                std::set<parse_int_type> pids;
                pids.insert(phrase);
            
                i++;
                l_right += l1_freq[phrase] - 1;
            
                if (i < l1_d.saD.size())
                {
                    auto new_sn = l1_d.saD[i];
                    auto new_phrase = l1_d.daD[i] + 1;
                    assert(new_phrase > 0 && new_phrase < l1_freq.size()); // + 1 because daD is 0-based
                    uint_t new_suffix_length = l1_d.select_b_d(l1_d.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
                
                    while (i < l1_d.saD.size() && (l1_d.lcpD[i] >= suffix_length) && (suffix_length == new_suffix_length))
                    {
                        chars.insert(l1_d.d[((new_sn + l1_d.d.size() - 1) % l1_d.d.size())]);
                        pids.insert(new_phrase);
                        ++i;
                    
                        l_right += l1_freq[new_phrase];
                    
                        if (i < l1_d.saD.size())
                        {
                            new_sn = l1_d.saD[i];
                            new_phrase = l1_d.daD[i] + 1;
                            assert(new_phrase > 0 && new_phrase < l1_freq.size()); // + 1 because daD is 0-based
                            new_suffix_length = l1_d.select_b_d(l1_d.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
                        }
                    }
                
                }
                
                std::size_t sa_values_in_this_range = 0;
                std::size_t sa_pos_f_it = sa_pos_iterator;
                while (sa_pos_f_it < runs_heads_positions.size() and runs_heads_positions[sa_pos_f_it] <= l_right)
                {
                    sa_values_in_this_range += 1;
                    sa_pos_f_it++;
                }
                
                if (sa_values_in_this_range != 0) // do the whole thing
                {
                    std::size_t sa_values_c = 0; // number of sa values computed in this range
                    
                    // shorthands
                    typedef std::reference_wrapper<std::vector<std::pair<std::size_t, std::size_t>>> ve_t; // (row in which that char appears, number of times per row)
                    typedef std::size_t ve_pq_t; // for the priority queue the row is enough, the other ingo can be retrieved from ve_t
                    typedef std::pair<ve_pq_t, std::pair<std::size_t, std::size_t>> pq_t; // .first : row, this element is the .second.second-th element of the .second.first array
    
                    // go in second level and iterate through list of positions for the meta-phrases we are interested into
                    std::vector<ve_t> v;
                    std::vector<parse_int_type> pids_v;
                    for (const auto& pid : pids)
                    {
                        parse_int_type adj_pid = pid + int_shift;
                        v.push_back(std::ref(l2_pfp_v_table[adj_pid]));
                        pids_v.push_back(adj_pid);
                    }
    
                    // make a priority queue and add elements to it
                    std::priority_queue<pq_t, std::vector<pq_t>, std::greater<>> pq;
                    for (std::size_t vi = 0; vi < v.size(); vi++) { pq.push({ v[vi].get()[0].first, { vi, 0 } }); }
    
                    while (not pq.empty())
                    {
                        // get all chars from this row entry at l2
                        std::vector<pq_t> from_same_l2_suffix;
                        auto first_suffix = pq.top();
                        while ((not pq.empty()) and (pq.top().first == first_suffix.first))
                        {
                            auto curr = pq.top(); pq.pop();
                            from_same_l2_suffix.push_back(curr);
            
                            std::size_t arr_i_c = curr.second.first;  // ith array
                            std::size_t arr_x_c = curr.second.second; // index in i-th array
                            if (arr_x_c + 1 < v[arr_i_c].get().size())
                            {
                                pq.push({ v[arr_i_c].get()[arr_x_c + 1].first, { arr_i_c, arr_x_c + 1 } });
                            }
                        }
        
                        // hard-hard suffix, need to look at the grid in L2
                        std::size_t curr_l2_row = from_same_l2_suffix[0].first;
                        auto& l2_M_entry = l2_pfp.M[curr_l2_row];
        
                        // get inverted lists of corresponding phrases
                        std::vector<std::reference_wrapper<std::vector<uint_t>>> ilists;
                        std::vector<std::reference_wrapper<std::vector<std::size_t>>> ilists_e_arrays;
                        std::vector<dict_l1_data_type> ilist_corresponding_chars;
                        std::vector<std::size_t> ilist_corresponding_sa_expanded_values;
                        for (std::size_t c_it = l2_M_entry.left; c_it <= l2_M_entry.right; c_it++)
                        {
                            // check if we need that pid
                            std::size_t l2_pid = l2_pfp.dict.colex_id[c_it] + 1;
                            parse_int_type adj_l1_pid = l2_pfp.dict.d[l2_pfp.dict.select_b_d(l2_pid + 1) - (l2_M_entry.len + 2)];
                            assert(adj_l1_pid >= int_shift);
                            parse_int_type l1_pid = adj_l1_pid - int_shift;
            
                            if (pids.contains(l1_pid))
                            {
                                ilists.push_back(std::ref(l2_pfp.bwt_p_ilist[l2_pid]));
                                ilists_e_arrays.push_back(std::ref(E_arrays[l2_pid - 1]));
                
                                dict_l1_data_type c = l1_d.d[l1_d.select_b_d(l1_pid + 1) - (suffix_length + 2)];
                                ilist_corresponding_chars.push_back(c);
                
                                // get the length of the current l2_suffix by expanding phrases
                                std::size_t l2_suff_start = l2_pfp.dict.select_b_d(l2_pid + 1) - 2 - (l2_M_entry.len - 1);
                                std::size_t l2_suff_end = l2_pfp.dict.select_b_d(l2_pid + 1) - 2 - (l2_pfp.w - 1) - 1;
                                std::size_t sa_r = 0;
                                for (std::size_t p_it_b = l2_suff_start; p_it_b <= l2_suff_end; p_it_b++)
                                {
                                    std::size_t adj_pid = l2_pfp.dict.d[p_it_b] - int_shift;
                                    sa_r += l1_d_lengths[adj_pid - 1] - l1_d.w; // l1_d.select_b_d(adj_pid + 1) - l1_d.select_b_d(adj_pid) - 1 - l1_d.w;
                                }
                                ilist_corresponding_sa_expanded_values.push_back(sa_r); // the sum of the lengsh of each l1 pid in the l2 phrase suffix
                            }
                        }
        
                        // make a priority queue from the inverted lists
                        typedef std::pair<uint_t, std::pair<std::size_t, std::size_t>> ilist_pq_t;
                        std::priority_queue<ilist_pq_t, std::vector<ilist_pq_t>, std::greater<>> ilist_pq;
                        for (std::size_t il_i = 0; il_i < ilists.size(); il_i++) { ilist_pq.push({ ilists[il_i].get()[0], { il_i, 0 } }); }
        
                        // now pop elements from the priority queue and write out the corresponding character
                        while (not ilist_pq.empty())
                        {
                            auto curr = ilist_pq.top(); ilist_pq.pop();
            
                            // compute corresponding sa value
                            std::size_t phrase_end = ilists_e_arrays[curr.second.first].get()[curr.second.second];
                            std::size_t out_sa_value = phrase_end - (ilist_corresponding_sa_expanded_values[curr.second.first]);
            
                            // adjust sa value for circular representation
                            if (out_sa_value > suffix_length) { out_sa_value -= suffix_length; }
                            else { out_sa_value = (l1_n + out_sa_value - suffix_length) % l1_n; }
                            
                            // output if run head
                            if (l_left + sa_values_c == runs_heads_positions[sa_pos_iterator])
                            {
                                out_sa_values.push_back(out_sa_value);
                                sa_pos_iterator += 1;
                            }
            
                            // keep iterating
                            std::size_t arr_i_c_il = curr.second.first;  // ith array
                            std::size_t arr_x_c_il = curr.second.second; // index in i-th array
                            if (arr_x_c_il + 1 < ilists[arr_i_c_il].get().size())
                            {
                                ilist_pq.push({ ilists[arr_i_c_il].get()[arr_x_c_il + 1], { arr_i_c_il, arr_x_c_il + 1 } });
                            }
                            
                            sa_values_c++;
                        }
                    }
                }
                else
                { sa_pos_iterator += sa_values_in_this_range; } // else skip to the next range
            
                l_left = l_right + 1;
                l_right = l_left;
                row++;
            }
        }
    }
    
    //------------------------------------------------------------
    
    void l1_refined_rindex(std::size_t threads = 1)
    {
        // Set threads accordingly to configuration
        omp_set_num_threads(threads);
        
        // Compute run length encoded bwt and run heads positions
        std::vector<std::vector<std::size_t>> run_start_positions(chunks.size());

        #pragma omp parallel for schedule(static) default(none) shared(run_start_positions)
        for (std::size_t i = 0; i < chunks.size(); i++)
        {
            l1_bwt_chunk(chunks[i], rle_chunks.get_encoder(i), run_start_positions[i]);
            spdlog::info("Chunk {}/{} completed", i + 1, chunks.size());
        }
        rle_chunks.close();
    
        // Compute sa values at run heads
        std::vector<std::vector<std::size_t>> run_start_sa_values(chunks.size());
        
        #pragma omp parallel for schedule(static) default(none) shared(run_start_positions, run_start_sa_values)
        for (std::size_t i = 0; i < chunks.size(); i++)
        {
            l1_sa_values_chunk(chunks[i], run_start_positions[i], run_start_sa_values[i]);
            spdlog::info("Chunk {}/{} completed", i + 1, chunks.size());
        }
        
        spdlog::info("Writing out sa values and positions");
        std::size_t start_pos_size = 0;
        for (std::size_t i = 0; i < chunks.size(); i++) { start_pos_size += run_start_positions[i].size(); }
    
        std::size_t sa_values_size = 0;
        for (std::size_t i = 0; i < chunks.size(); i++) { sa_values_size += run_start_sa_values[i].size(); }
        
        assert(start_pos_size == sa_values_size);
        
        // write out first size, then content
        sa_positions_ofstream.write((char*)&start_pos_size, sizeof(start_pos_size));
        for (std::size_t i = 0; i < chunks.size(); i++)
        {
            sa_positions_ofstream.write((char*)&(run_start_positions[i][0]), run_start_positions[i].size() * sizeof(std::size_t));
        }
        sa_positions_ofstream.close();
    
        // write out first size, then content
        sa_values_ofstream.write((char*)&sa_values_size, sizeof(sa_values_size));
        for (std::size_t i = 0; i < chunks.size(); i++)
        {
            sa_values_ofstream.write((char*)&(run_start_sa_values[i][0]), run_start_sa_values[i].size() * sizeof(std::size_t));
        }
        sa_values_ofstream.close();
    }
    
};

}


#endif //rpfbwt_algorithm_hpp
