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
        
        bool operator()(parse_int_type l, parse_int_type r) const
        {
            assert(l != r);
            if (l < int_shift) { return true; }
            if (r < int_shift) { return false; }
            return l1_d.colex_id[l - int_shift - 1] < l1_d.colex_id[r - int_shift - 1];
        }
    };
    
private:
    
    pfpds::long_type int_shift = 10;
    std::string l1_prefix;
    std::string out_rle_name;
    
    std::less<dict_l1_data_type> l1_d_comp;
    const pfpds::dictionary<dict_l1_data_type>& l1_d;
    std::vector<pfpds::long_type> l1_d_lengths;
    std::vector<pfpds::long_type> l1_freq;
    
    l2_colex_comp l2_comp;
    const pfpds::dictionary<parse_int_type, l2_colex_comp>& l2_d;
    const pfpds::parse& l2_p;
    const pfpds::pf_parsing<parse_int_type, l2_colex_comp, pfpds::pfp_wt_sdsl>& l2_pfp;
    std::vector<pfpds::long_type> l2_d_lengths;
    std::vector<std::vector<pfpds::long_type>> E_arrays;
    pfpds::long_type l1_n = 0; // also used to initialize l1_n
    
    std::vector<std::vector<std::pair<pfpds::long_type, pfpds::long_type>>> l2_pfp_v_table; // (row in which that char appears, number of times per row)

    static constexpr pfpds::long_type default_num_of_chunks = 1;
    std::vector<std::tuple<pfpds::long_type, pfpds::long_type, pfpds::long_type, pfpds::long_type>> chunks;
    
    rle::RLEString::RLEncoderMerger rle_chunks;
    
    template <typename dict_char_type>
    void
    init_d_lengths(const std::vector<dict_char_type>& dict_array, std::vector<pfpds::long_type>& out)
    {
        spdlog::info("Precomputing d lengths");
        
        pfpds::long_type curr_length = 0;
        for (pfpds::long_type i = 0; i < dict_array.size(); i++)
        {
            curr_length += 1;
            if (dict_array[i] == EndOfWord) { out.push_back(curr_length - 1); curr_length = 0; }
        }
    }
    
    void
    init_v_table()
    {
        spdlog::info("Creating v table");
        
        // now build V
        for (pfpds::long_type row = 0; row < l2_pfp.M.size(); row++)
        {
            const auto& m =  l2_pfp.M[row];
        
            std::vector<std::pair<parse_int_type, pfpds::long_type>> phrase_counts;
            for (pfpds::long_type r = m.left; r <= m.right; r++)
            {
                auto phrase = l2_pfp.dict.colex_id[r];
                pfpds::long_type phrase_start = l2_pfp.dict.select_b_d(phrase + 1);
                pfpds::long_type phrase_length =  l2_d_lengths[phrase];
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
                std::pair<pfpds::long_type, pfpds::long_type> v_element = std::make_pair(row, c_count.second);
                l2_pfp_v_table[c_count.first].emplace_back(v_element);
            }
        }
    }
    
    pfpds::long_type
    init_E_arrays()
    {
        spdlog::info("Creating E arrays");
        
        E_arrays.resize(l2_pfp.dict.n_phrases());
        pfpds::long_type end_pos_from_start = 0;
        
        // read P2 and expand each phrase to get the end position in the text order
        std::vector<pfpds::long_type> end_positions;
        for (pfpds::long_type i = 0; i < l2_pfp.pars.p.size() - 1; i++)
        {
            pfpds::long_type l2_pid = l2_pfp.pars.p[i];
            pfpds::long_type l2_phrase_start = l2_pfp.dict.select_b_d(l2_pid);
            pfpds::long_type l2_phrase_end = l2_phrase_start + l2_d_lengths[l2_pid - 1] - 1;
            if (l2_phrase_start < l2_pfp.w) { l2_phrase_start = l2_pfp.w; }
            
            for (pfpds::long_type j = l2_phrase_start; j <= (l2_phrase_end - l2_pfp.w); j++)
            {
                pfpds::long_type l1_pid = l2_pfp.dict.d[j] - int_shift;
                pfpds::long_type l1_phrase_length = l1_d_lengths[l1_pid - 1];
    
                end_pos_from_start += (l1_phrase_length - l1_d.w);
            }
            end_positions.push_back(end_pos_from_start);
        }
        
        // now reorder according to P2's SA
        for (pfpds::long_type i = 0; i < l2_pfp.pars.saP.size(); i++)
        {
            pfpds::long_type sa_value = l2_pfp.pars.saP[i];
            
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
    void
    compute_chunks(pfpds::long_type num_of_chunks)
    {
        spdlog::info("Computing chunks for parallel execution. Total input size: {}. Requested chunks: {}", l1_n, num_of_chunks);
        
        pfpds::long_type chunk_size = (num_of_chunks > 1) ? (l1_n / (num_of_chunks)) : l1_n + 1;
        
        // Go through the suffixes of D and compute chunks
        pfpds::long_type i = 1; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
        
        pfpds::long_type l_left  = 0;
        pfpds::long_type l_right = 0;
        pfpds::long_type chunk_suffix_start = i;
        pfpds::long_type chunk_start = l_left;
        pfpds::long_type chunk_row_start = 0;
        pfpds::long_type table_row = 0;
        while (i < l1_d.saD.size())
        {
            auto sn = l1_d.saD[i];
            // Check if the suffix has length at least w and is not the complete phrase.
            auto phrase = l1_d.daD[i] + 1;
            assert(phrase > 0 && phrase < l1_freq.size()); // + 1 because daD is 0-based
            pfpds::long_type suffix_length = l1_d.select_b_d(l1_d.rank_b_d(sn + 1) + 1) - sn - 1;
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
                    pfpds::long_type new_suffix_length = l1_d.select_b_d(l1_d.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
                    
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
                    chunks.emplace_back(chunk_suffix_start, i, chunk_start, chunk_row_start);
                    
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
        if ((chunks.empty()) or (std::get<1>(chunks.back()) < i))
        {
            // store SA_d range for last chunk
            chunks.emplace_back(chunk_suffix_start, i, chunk_start, chunk_row_start);
        }
        
        assert(chunks.size() == num_of_chunks);
        spdlog::info("BWT ranges divided in {} chunks.", chunks.size());
    }
    
public:
    
    rpfbwt_algo(
    std::string prefix,
    const pfpds::dictionary<dict_l1_data_type>& l1_d_r,
    const pfpds::pf_parsing<parse_int_type, l2_colex_comp, pfpds::pfp_wt_sdsl>& l2_pfp_r,
    pfpds::long_type pfp_integer_shift,
    pfpds::long_type bwt_chunks = default_num_of_chunks)
    :
    int_shift(pfp_integer_shift),
    l1_d(l1_d_r),
    l1_freq(l1_d.n_phrases() + 1, 0),
    l1_prefix(prefix),
    out_rle_name(l1_prefix + ".rlebwt"),
    l2_comp(l1_d, int_shift),
    l2_d(l2_pfp_r.dict),
    l2_p(l2_pfp_r.pars),
    l2_pfp(l2_pfp_r),
    l2_pfp_v_table(l2_d.alphabet_size),
    rle_chunks(out_rle_name, bwt_chunks)
    {
        init_d_lengths<dict_l1_data_type>(l1_d.d, l1_d_lengths);
        init_d_lengths<parse_int_type>(l2_d.d, l2_d_lengths);
        
        // compute frequencies
        for (pfpds::long_type i = 0; i < l2_pfp.pars.p.size() - 1; i++)
        {
            pfpds::long_type l2_pid = l2_pfp.pars.p[i];
            pfpds::long_type l2_phrase_start = l2_pfp.dict.select_b_d(l2_pid);
            pfpds::long_type l2_phrase_end = l2_phrase_start + l2_d_lengths[l2_pid - 1] - 1;
            if (l2_phrase_start < l2_pfp.w) { l2_phrase_start = l2_pfp.w; }
        
            for (pfpds::long_type j = l2_phrase_start; j <= (l2_phrase_end - l2_pfp.w); j++)
            {
                pfpds::long_type l1_pid = l2_pfp.dict.d[j] - int_shift;
                l1_freq[l1_pid] += 1;
            }
        }
        
        init_v_table();
        this->l1_n = init_E_arrays();
        compute_chunks(bwt_chunks);
    }
    
    void
    l1_bwt_chunk(
    const std::tuple<pfpds::long_type, pfpds::long_type, pfpds::long_type, pfpds::long_type>& chunk,
    rle::RLEString::RLEncoder& rle_out)
    {
        pfpds::long_type i = std::get<0>(chunk);
        
        dict_l1_data_type prev_char = 0;
    
        pfpds::long_type l_left  = std::get<2>(chunk);
        pfpds::long_type l_right = l_left;
        pfpds::long_type easy_chars = 0;
        pfpds::long_type hard_easy_chars = 0;
        pfpds::long_type hard_hard_chars = 0;
        pfpds::long_type row = std::get<3>(chunk);
        while (i < std::get<1>(chunk))
        {
            pfpds::long_type sn = l1_d.saD[i];
            // Check if the suffix has length at least w and is not the complete phrase.
            pfpds::long_type phrase = l1_d.daD[i] + 1;
            assert(phrase > 0 && phrase < l1_freq.size()); // + 1 because daD is 0-based
            pfpds::long_type suffix_length = l1_d.select_b_d(l1_d.rank_b_d(sn + 1) + 1) - sn - 1;
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
                    pfpds::long_type new_sn = l1_d.saD[i];
                    pfpds::long_type new_phrase = l1_d.daD[i] + 1;
                    assert(new_phrase > 0 && new_phrase < l1_freq.size()); // + 1 because daD is 0-based
                    pfpds::long_type new_suffix_length = l1_d.select_b_d(l1_d.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
                
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
                    easy_chars += (l_right - l_left) + 1;
                }
                else // easy-hard and hard-hard suffixes
                {
                    assert(l_right - l_left > 0);
                    pfpds::long_type hard_chars_before = hard_easy_chars + hard_hard_chars;
                    
                    // shorthands
                    typedef std::reference_wrapper<std::vector<std::pair<pfpds::long_type, pfpds::long_type>>> ve_t; // (row in which that char appears, number of times per row)
                    typedef pfpds::long_type ve_pq_t; // for the priority queue the row is enough, the other ingo can be retrieved from ve_t
                    typedef std::pair<ve_pq_t, std::pair<pfpds::long_type, pfpds::long_type>> pq_t; // .first : row, this element is the .second.second-th element of the .second.first array
                    
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
                    std::priority_queue<pq_t, std::vector<pq_t>, std::greater<pq_t>> pq;
                    for (pfpds::long_type vi = 0; vi < v.size(); vi++) { pq.push({ v[vi].get()[0].first, { vi, 0 } }); }

                    while (not pq.empty())
                    {
                        // get all chars from this row entry at l2
                        std::vector<pq_t> from_same_l2_suffix;
                        auto first_suffix = pq.top();
                        while ((not pq.empty()) and (pq.top().first == first_suffix.first))
                        {
                            auto curr = pq.top(); pq.pop();
                            from_same_l2_suffix.push_back(curr);

                            pfpds::long_type arr_i_c = curr.second.first;  // ith array
                            pfpds::long_type arr_x_c = curr.second.second; // index in i-th array
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
                            pfpds::long_type freq = v[pq_entry.second.first].get()[pq_entry.second.second].second;
                            dict_l1_data_type c = l1_d.d[l1_d.select_b_d(adj_pid - int_shift + 1) - (suffix_length + 2)]; // end of next phrases
                            
                            rle_out(c, freq);
                            hard_easy_chars += freq;
                        }
                        else
                        {
                            // hard-hard suffix, need to look at the grid in L2
                            pfpds::long_type curr_l2_row = from_same_l2_suffix[0].first;
                            auto l2_M_entry = l2_pfp.M[curr_l2_row];
                            
                            // get inverted lists of corresponding phrases
                            std::vector<std::reference_wrapper<const std::vector<pfpds::long_type>>> ilists;
                            std::vector<dict_l1_data_type> ilist_corresponding_chars;
                            for (pfpds::long_type c_it = l2_M_entry.left; c_it <= l2_M_entry.right; c_it ++)
                            {
                                // check if we need that pid
                                pfpds::long_type l2_pid = l2_pfp.dict.colex_id[c_it] + 1;
                                parse_int_type adj_l1_pid = l2_pfp.dict.d[l2_pfp.dict.select_b_d(l2_pid + 1) - (l2_M_entry.len + 2)];
                                assert(adj_l1_pid >= int_shift);
                                parse_int_type l1_pid = adj_l1_pid - int_shift;
    
                                if (pids.find(l1_pid) != pids.end())
                                {
                                    ilists.push_back(std::cref(l2_pfp.bwt_p_ilist[l2_pfp.dict.colex_id[c_it] + 1]));
                                    dict_l1_data_type c = l1_d.d[l1_d.select_b_d(l1_pid + 1) - (suffix_length + 2)];
                                    ilist_corresponding_chars.push_back(c);
                                }
                            }
    
                            // make a priority queue from the inverted lists
                            typedef std::pair<pfpds::long_type, std::pair<pfpds::long_type, pfpds::long_type>> ilist_pq_t;
                            std::priority_queue<ilist_pq_t, std::vector<ilist_pq_t>, std::greater<ilist_pq_t>> ilist_pq;
                            for (pfpds::long_type il_i = 0; il_i < ilists.size(); il_i++) { ilist_pq.push({ ilists[il_i].get()[0], { il_i, 0 } }); }
                            
                            // now pop elements from the priority queue and write out the corresponding character
                            while (not ilist_pq.empty())
                            {
                                auto curr = ilist_pq.top(); ilist_pq.pop();
                                
                                // output corresponding char
                                dict_l1_data_type c = ilist_corresponding_chars[curr.second.first];

                                rle_out(c, 1);
                                hard_hard_chars += 1;
                                
                                // keep iterating
                                pfpds::long_type arr_i_c_il = curr.second.first;  // ith array
                                pfpds::long_type arr_x_c_il = curr.second.second; // index in i-th array
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
    const std::tuple<pfpds::long_type, pfpds::long_type, pfpds::long_type, pfpds::long_type>& chunk,
    const rle::bitvector& run_heads_bitvector,
    const std::string& tmp_file_name)
    {
        pfpds::long_type skipped_portion = 0;
        pfpds::long_type processed_portion = 0;
        
        pfpds::long_type i = std::get<0>(chunk);
        
        std::ofstream out_sa_tmp_fstream(tmp_file_name, std::ios::out | std::ios::binary);
    
        pfpds::long_type sa_pos_iterator = 0;
        
        pfpds::long_type l_left  = std::get<2>(chunk);
        pfpds::long_type l_right = l_left;
        pfpds::long_type easy_chars = 0;
        pfpds::long_type hard_easy_chars = 0;
        pfpds::long_type hard_hard_chars = 0;
        pfpds::long_type row = std::get<3>(chunk);
        while (i < std::get<1>(chunk))
        {
            auto sn = l1_d.saD[i];
            // Check if the suffix has length at least w and is not the complete phrase.
            auto phrase = l1_d.daD[i] + 1;
            assert(phrase > 0 && phrase < l1_freq.size()); // + 1 because daD is 0-based
            pfpds::long_type suffix_length = l1_d.select_b_d(l1_d.rank_b_d(sn + 1) + 1) - sn - 1;
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
                    pfpds::long_type new_suffix_length = l1_d.select_b_d(l1_d.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
                
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
                
                pfpds::long_type run_heads_in_range = 0;
                if (run_heads_bitvector[l_left]) { run_heads_in_range = 1; } // at least
                else { run_heads_in_range = run_heads_bitvector.rank(l_right + 1) - run_heads_bitvector.rank(l_left); } ; // include the extremes
                
                if (run_heads_in_range != 0) // do the whole thing
                {
                    pfpds::long_type sa_values_c = 0; // number of sa values computed in this range
                    
                    // shorthands
                    typedef std::reference_wrapper<std::vector<std::pair<pfpds::long_type, pfpds::long_type>>> ve_t; // (row in which that char appears, number of times per row)
                    typedef pfpds::long_type ve_pq_t; // for the priority queue the row is enough, the other ingo can be retrieved from ve_t
                    typedef std::pair<ve_pq_t, std::pair<pfpds::long_type, pfpds::long_type>> pq_t; // .first : row, this element is the .second.second-th element of the .second.first array
    
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
                    std::priority_queue<pq_t, std::vector<pq_t>, std::greater<pq_t>> pq;
                    for (pfpds::long_type vi = 0; vi < v.size(); vi++) { pq.push({ v[vi].get()[0].first, { vi, 0 } }); }
    
                    while (not pq.empty())
                    {
                        // get all chars from this row entry at l2
                        std::vector<pq_t> from_same_l2_suffix;
                        auto first_suffix = pq.top();
                        while ((not pq.empty()) and (pq.top().first == first_suffix.first))
                        {
                            auto curr = pq.top(); pq.pop();
                            from_same_l2_suffix.push_back(curr);
            
                            pfpds::long_type arr_i_c = curr.second.first;  // ith array
                            pfpds::long_type arr_x_c = curr.second.second; // index in i-th array
                            if (arr_x_c + 1 < v[arr_i_c].get().size())
                            {
                                pq.push({ v[arr_i_c].get()[arr_x_c + 1].first, { arr_i_c, arr_x_c + 1 } });
                            }
                        }
        
                        // hard-hard suffix, need to look at the grid in L2
                        pfpds::long_type curr_l2_row = from_same_l2_suffix[0].first;
                        auto& l2_M_entry = l2_pfp.M[curr_l2_row];
        
                        // get inverted lists of corresponding phrases
                        std::vector<std::reference_wrapper<const std::vector<pfpds::long_type>>> ilists;
                        std::vector<std::reference_wrapper<const std::vector<pfpds::long_type>>> ilists_e_arrays;
                        std::vector<pfpds::long_type> ilist_corresponding_sa_expanded_values;
                        for (pfpds::long_type c_it = l2_M_entry.left; c_it <= l2_M_entry.right; c_it++)
                        {
                            // check if we need that pid
                            pfpds::long_type l2_pid = l2_pfp.dict.colex_id[c_it] + 1;
                            parse_int_type adj_l1_pid = l2_pfp.dict.d[l2_pfp.dict.select_b_d(l2_pid + 1) - (l2_M_entry.len + 2)];
                            assert(adj_l1_pid >= int_shift);
                            parse_int_type l1_pid = adj_l1_pid - int_shift;
    
                            if (pids.find(l1_pid) != pids.end())
                            {
                                ilists.push_back(std::cref(l2_pfp.bwt_p_ilist[l2_pid]));
                                ilists_e_arrays.push_back(std::ref(E_arrays[l2_pid - 1]));
                                
                                // get the length of the current l2_suffix by expanding phrases
                                pfpds::long_type l2_suff_start = l2_pfp.dict.select_b_d(l2_pid + 1) - 2 - (l2_M_entry.len - 1);
                                pfpds::long_type l2_suff_end = l2_pfp.dict.select_b_d(l2_pid + 1) - 2 - (l2_pfp.w - 1) - 1;
                                pfpds::long_type sa_r = 0;
                                for (pfpds::long_type p_it_b = l2_suff_start; p_it_b <= l2_suff_end; p_it_b++)
                                {
                                    pfpds::long_type adj_pid = l2_pfp.dict.d[p_it_b] - int_shift;
                                    sa_r += l1_d_lengths[adj_pid - 1] - l1_d.w; // l1_d.select_b_d(adj_pid + 1) - l1_d.select_b_d(adj_pid) - 1 - l1_d.w;
                                }
                                ilist_corresponding_sa_expanded_values.push_back(sa_r); // the sum of the lengsh of each l1 pid in the l2 phrase suffix
                            }
                        }
        
                        // make a priority queue from the inverted lists
                        typedef std::pair<pfpds::long_type, std::pair<pfpds::long_type, pfpds::long_type>> ilist_pq_t;
                        std::priority_queue<ilist_pq_t, std::vector<ilist_pq_t>, std::greater<ilist_pq_t>> ilist_pq;
                        for (pfpds::long_type il_i = 0; il_i < ilists.size(); il_i++) { ilist_pq.push({ ilists[il_i].get()[0], { il_i, 0 } }); }
        
                        // now pop elements from the priority queue and write out the corresponding character
                        while (not ilist_pq.empty())
                        {
                            auto curr = ilist_pq.top(); ilist_pq.pop();
            
                            // compute corresponding sa value
                            pfpds::long_type phrase_end = ilists_e_arrays[curr.second.first].get()[curr.second.second];
                            pfpds::long_type out_sa_value = phrase_end - (ilist_corresponding_sa_expanded_values[curr.second.first]);
            
                            // adjust sa value for circular representation
                            if (out_sa_value > suffix_length) { out_sa_value -= suffix_length; }
                            else { out_sa_value = (l1_n + out_sa_value - suffix_length) % l1_n; }
                            
                            // output if run head
                            if (run_heads_bitvector[l_left + sa_values_c])
                            {
                                out_sa_tmp_fstream.write((char*)&out_sa_value, sizeof(pfpds::long_type));
                                sa_pos_iterator += 1;
                            }
            
                            // keep iterating
                            pfpds::long_type arr_i_c_il = curr.second.first;  // ith array
                            pfpds::long_type arr_x_c_il = curr.second.second; // index in i-th array
                            if (arr_x_c_il + 1 < ilists[arr_i_c_il].get().size())
                            {
                                ilist_pq.push({ ilists[arr_i_c_il].get()[arr_x_c_il + 1], { arr_i_c_il, arr_x_c_il + 1 } });
                            }
                            
                            sa_values_c++;
                        }
                    }
                    processed_portion += l_right - l_left;
                }
                else // we skipped this portions
                {
                    skipped_portion += l_right - l_left;
                }
            
                l_left = l_right + 1;
                l_right = l_left;
                row++;
            }
        }
        spdlog::info("Processed: {} Skipped: {}", processed_portion, skipped_portion);
    
        // close tmp file stream
        out_sa_tmp_fstream.close();
    }
    
    //------------------------------------------------------------
    
    void l1_rlebwt(pfpds::long_type threads = 1)
    {
        // Set threads accordingly to configuration
        omp_set_num_threads(threads);
    
        // ----------
        // Compute run length encoded bwt and run heads positions
        #pragma omp parallel for schedule(dynamic) default(none)
        for (pfpds::long_type i = 0; i < chunks.size(); i++)
        {
            l1_bwt_chunk(chunks[i], rle_chunks.get_encoder(i));
            spdlog::info("Chunk {}/{} completed", i + 1, chunks.size());
        }
        rle_chunks.close();
    }
    
    
    //------------------------------------------------------------
    
    void l1_refined_rindex(pfpds::long_type threads = 1)
    {
        // Set threads accordingly to configuration
        omp_set_num_threads(threads);
        
        // ----------
        // Compute run length encoded bwt and run heads positions
        #pragma omp parallel for schedule(dynamic) default(none)
        for (pfpds::long_type i = 0; i < chunks.size(); i++)
        {
            l1_bwt_chunk(chunks[i], rle_chunks.get_encoder(i));
            spdlog::info("Chunk {}/{} completed", i + 1, chunks.size());
        }
        rle_chunks.close();
    
        // ----------
        // Get bitvector marking run heads from the rle bwt
        rle::RLEString::RLEDecoder rle_decoder(out_rle_name);
        sdsl::sd_vector_builder run_heads_bv_builder(rle_decoder.metadata.size, rle_decoder.metadata.runs);
        pfpds::long_type runs_bv_it = 0;
    
        while (not rle_decoder.end())
        {
            rle::RunType run = rle_decoder.next();
            dict_l1_data_type c = rle::RunTraits::charachter(run);
            pfpds::long_type len = rle::RunTraits::length(run);
        
            if (len != 0)
            { run_heads_bv_builder.set(runs_bv_it); runs_bv_it += len; }
        }
    
        rle::bitvector run_heads_bv = rle::bitvector(run_heads_bv_builder);
    
        // ----------
        // Compute sa values at run heads
        std::vector<std::string> sa_chunks_tmp_files;
        for (pfpds::long_type i = 0; i < chunks.size(); i++) { sa_chunks_tmp_files.push_back(rle::TempFile::getName("sa_chunk")); }
        
        #pragma omp parallel for schedule(dynamic) default(none) shared(run_heads_bv, sa_chunks_tmp_files)
        for (pfpds::long_type i = 0; i < chunks.size(); i++)
        {
            l1_sa_values_chunk(chunks[i], run_heads_bv, sa_chunks_tmp_files[i]);
            spdlog::info("Chunk {}/{} completed", i + 1, chunks.size());
        }
    
        // ----------
        // Merge sa values tmp files
        spdlog::info("Writing out sa values at run heads");
        std::ofstream out_sa_values(l1_prefix + ".ssa", std::ios::out | std::ios::binary);
        
        // write number of values first
        pfpds::long_type tot_sa_values = run_heads_bv.number_of_1();
        out_sa_values.write((char*) &tot_sa_values, sizeof(tot_sa_values));
        
        for (pfpds::long_type i = 0; i < chunks.size(); i++)
        {
            std::ifstream if_sa_chunk(sa_chunks_tmp_files[i], std::ios_base::binary);
            out_sa_values << if_sa_chunk.rdbuf();
        }
        out_sa_values.close();
    
        spdlog::info("Done");
    }
    
};

}


#endif //rpfbwt_algorithm_hpp
