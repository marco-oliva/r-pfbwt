//
//  rpfbwt_algorithm.hpp
//
//  Copyright 2022 Marco Oliva. All rights reserved.
//

#ifndef rpfbwt_algorithm_hpp
#define rpfbwt_algorithm_hpp

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
    std::vector<uint_t> l1_freq; // here occ has the same size as the integers used for gsacak.
    bool l1_cleared = false;
    
    l2_colex_comp l2_comp;
    pfpds::pf_parsing<parse_int_type, l2_colex_comp, pfpds::pfp_wt_sdsl> l2_pfp;
    bool l2_cleared = false;
    
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> l2_pfp_v_table; // (row in which that char appears, number of times per row)
    
    static const std::size_t chunk_size_default = 50;
    std::size_t chunk_size;
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> chunks;
    
    rle::RLEString::RLEncoderMerger rle_chunks;
    
    void init_v_table()
    {
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
                if (c > int_shift) { c -= int_shift; }
            
                if (phrase_counts.empty() or phrase_counts.back().first != c)
                {
                    phrase_counts.push_back(std::make_pair(c, 1));
                }
                else { phrase_counts.back().second += 1; }
            }
        
            // update
            for (auto const& c_count : phrase_counts)
            {
                std::pair<std::size_t, std::size_t> v_element = std::make_pair(row, c_count.second);
                l2_pfp_v_table[c_count.first].emplace_back(v_element);
            }
        }
    }
    
    // Computes ranges for parallel computation
    // suffix start, suffix end, this_left, this_row
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>
    compute_chunks(std::size_t chunk_size)
    {
        spdlog::info("Computing chunks for parallel execution.");
        
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
        
        
        if (out.empty()) { out.emplace_back(1, l1_d.saD.size(), 0, 0); }
        spdlog::info("BWT ranges divided in {} chunks.", out.size());
        if (out.size() > 500)
        {
            spdlog::error("Each chunk creates a temporary file on disk. To now overwhelm the filesystem increase the chunk size.");
            std::exit(1);
        }
        return out;
    }
    
    std::vector<uint_t> read_l1_freq(std::string& l1_prefix)
    {
        std::vector<uint_t> out;
        std::size_t d1_words; uint_t * occ;
        pfpds::read_file<uint_t> (std::string(l1_prefix + ".occ").c_str(), occ, d1_words);
        out.push_back(0);
        out.insert(out.end(),occ, occ + d1_words);
        delete[] occ;
        
        return out;
    }
    
    
    bool clear_L1_unused_data_structures()
    {
        spdlog::info("Removing unused L1 data structures.");
        l1_d.inv_colex_id.clear();
        l1_d.colex_daD.clear();
        l1_d.isaD.clear();
        return true;
    }
    
    bool clear_L2_unused_data_structures()
    {
        spdlog::info("Removing unused L2 data structures.");
        l2_pfp.dict.inv_colex_id.clear();
        l2_pfp.dict.colex_daD.clear();
        l2_pfp.dict.lcpD.clear();
        l2_pfp.dict.daD.clear();
        l2_pfp.dict.isaD.clear();
        l2_pfp.pars.isaP.clear();
        l2_pfp.pars.saP.clear();
        l2_pfp.pars.lcpP.clear();
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
                std::size_t bwt_chunk_size = chunk_size_default)
    : l1_d(l1_d_v, l1_w, l1_d_comp, true, false, true, true, false, true), l1_freq(l1_freq_v),
      l1_cleared(clear_L1_unused_data_structures()),
      l1_prefix("mem"), out_rle_name(l1_prefix + ".rlebwt"),
      l2_comp(l1_d, int_shift), l2_pfp(l2_d_v, l2_comp, l2_p_v, l2_freq_v, l2_w, int_shift),
      l2_cleared(clear_L2_unused_data_structures()),
      l2_pfp_v_table(l2_pfp.dict.alphabet_size),
      chunk_size(bwt_chunk_size), chunks(compute_chunks(chunk_size)), rle_chunks(out_rle_name, chunks.size())
    { init_v_table(); }
    
    rpfbwt_algo(std::string& l1_prefix, std::size_t l1_w, std::size_t l2_w,  std::size_t bwt_chunk_size = chunk_size_default)
    : l1_d(l1_prefix, l1_w, l1_d_comp, true, false, true, true, false, true),
      l1_prefix(l1_prefix), out_rle_name(l1_prefix + ".rlebwt"),
      l1_freq(read_l1_freq(l1_prefix)),
      l1_cleared(clear_L1_unused_data_structures()),
      l2_comp(l1_d, int_shift),
      l2_pfp(l1_prefix + ".parse", l2_w, l2_comp, int_shift), l2_pfp_v_table(l2_pfp.dict.alphabet_size),
      l2_cleared(clear_L2_unused_data_structures()),
      chunk_size(bwt_chunk_size), chunks(compute_chunks(chunk_size)), rle_chunks(out_rle_name, chunks.size())
    {
        init_v_table();
    }
    
    std::vector<dict_l1_data_type> l1_bwt_chunk(
        std::tuple<std::size_t, std::size_t, std::size_t, std::size_t> chunk,
        rle::RLEString::RLEncoder& rle_out,
        bool out_vector = false)
    {
        std::vector<dict_l1_data_type> out;
        
        std::size_t i = std::get<0>(chunk);
    
        std::size_t l_left  = std::get<2>(chunk);
        std::size_t l_right = l_left;
        std::size_t easy_chars = 0;
        std::size_t hard_easy_chars = 0;
        std::size_t hard_hard_chars = 0;
        std::size_t row = std::get<3>(chunk);
        while (i < std::get<1>(chunk))
        {
            std::size_t left = i;
        
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
                    if (out_vector) { out.insert(out.end(), (l_right - l_left) + 1, *(chars.begin())); }
                    rle_out(*(chars.begin()), (l_right - l_left) + 1);
                    easy_chars += (l_right - l_left) + 1;
                }
                else // easy-hard and hard-hard suffixes
                {
                    // shorthands
                    typedef std::reference_wrapper<std::vector<std::pair<std::size_t, std::size_t>>> ve_t; // (row in which that char appears, number of times per row)
                    typedef std::size_t ve_pq_t; // for the priprity queue the row is enough, the other ingo can be retrieved from ve_t
                    typedef std::pair<ve_pq_t, std::pair<std::size_t, std::size_t>> pq_t; // .first : row, this element is the .second.second-th element of the .second.first array
                    
                    // go in second level and iterate trough list of positions for the meta-phrases we are interested into
                    std::vector<ve_t> v;
                    std::vector<parse_int_type> pids_v;
                    for (const auto& pid : pids)
                    {
                       v.push_back(std::ref(l2_pfp_v_table[pid]));
                       pids_v.push_back(pid);
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
                            parse_int_type pid = pids_v[pq_entry.second.first];
                            uint_t freq = v[pq_entry.second.first].get()[pq_entry.second.second].second;
                            dict_l1_data_type c = l1_d.d[l1_d.select_b_d(pid + 1) - (suffix_length + 2)]; // end of next phrases
                            if (out_vector) { out.insert(out.end(), freq, c); }
                            rle_out(c, freq);
                            hard_easy_chars += freq;
                        }
                        else
                        {
                            // hard-hard suffix, need to look at the grid in L2
                            std::size_t curr_l2_row = from_same_l2_suffix[0].first;
                            auto l2_M_entry = l2_pfp.M[curr_l2_row];
                            
                            auto points = l2_pfp.w_wt.range_search_2d(l2_M_entry.left, l2_M_entry.right, l2_pfp.w_wt.size());
                            std::vector<dict_l1_data_type> hard_hard_chars_v(points.size(), 0);
                            
                            for (auto& point : points)
                            {
                                std::size_t colex_id = point.second;
                                std::size_t l2_pid = l2_pfp.dict.colex_id[colex_id] + 1;

                                parse_int_type l1_pid = l2_pfp.dict.d[l2_pfp.dict.select_b_d(l2_pid + 1) - (l2_M_entry.len + 2)];
                                
                                // check if l1_pid is among the ones we are looking for
                                if (l1_pid >= l2_pfp.shift) { l1_pid -= l2_pfp.shift; }
                                if (pids.contains(l1_pid))
                                {
                                    dict_l1_data_type c = l1_d.d[l1_d.select_b_d(l1_pid + 1) - (suffix_length + 2)];
                                    if (out_vector) { hard_hard_chars_v[point.first - l2_M_entry.left] = c; }
                                    rle_out(c, 1);
                                    hard_hard_chars += 1;
                                }
                            }
                            
                            // output chars in order
                            if (out_vector) {  for (auto& hc : hard_hard_chars_v) { if (hc != 0) { out.push_back(hc); } } }
                        }
                    }

                }
                
                l_left = l_right + 1;
                l_right = l_left;
                row++;
            }
        }

        return out;
    }
    
    std::vector<dict_l1_data_type> l1_bwt(bool out_vector = false)
    {
        std::vector<dict_l1_data_type> out;
        
        for (std::size_t i = 0; i < chunks.size(); i++)
        {
            std::vector<dict_l1_data_type> bwt_chunk = l1_bwt_chunk(chunks[i], rle_chunks.get_encoder(i), out_vector);
            if (not bwt_chunk.empty()) { out.insert(out.end(), bwt_chunk.begin(), bwt_chunk.end()); }
            
        }
        
        return out;
    }
    
    
    void l1_bwt_parallel()
    {
        #pragma omp parallel for schedule(static) default(none)
        for (std::size_t i = 0; i < chunks.size(); i++)
        {
            l1_bwt_chunk(chunks[i], rle_chunks.get_encoder(i));
            spdlog::info("Chunk {}/{} completed", i, chunks.size());
        }
    }
};

}


#endif //rpfbwt_algorithm_hpp
