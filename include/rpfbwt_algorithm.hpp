//
//  rpfbwt_algorithm.hpp
//
//  Copyright 2022 Marco Oliva. All rights reserved.
//

#ifndef rpfbwt_algorithm_hpp
#define rpfbwt_algorithm_hpp

#include <pfp/pfp.hpp>

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
            if (l >= int_shift) { l -= int_shift; }
            if (r >= int_shift) { r -= int_shift; }
            return l1_d.colex_id[l - 1] < l1_d.colex_id[r - 1];
        }
    };
    
public: // TODO: back to private
    
    const std::size_t int_shift = 10;
    std::string l1_prefix;
    
    std::less<dict_l1_data_type> l1_d_comp;
    pfpds::dictionary<dict_l1_data_type> l1_d;
    std::vector<uint_t> l1_freq; // here occ has the same size as the integers used for gsacak.
    
    l2_colex_comp l2_comp;
    pfpds::pf_parsing<parse_int_type, l2_colex_comp, pfpds::pfp_wt_sdsl> l2_pfp;
    
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> l2_pfp_v_table; // (row in which that char appears, number of times per row)
    
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

public:
    
    rpfbwt_algo(std::vector<dict_l1_data_type>& l1_d_v,
                std::vector<uint_t>& l1_freq_v,
                std::size_t l1_w,
                std::vector<parse_int_type>& l2_d_v,
                std::vector<parse_int_type>& l2_p_v,
                std::vector<uint_t>& l2_freq_v,
                std::size_t l2_w)
    : l1_d(l1_d_v, l1_w, l1_d_comp, true, false, true, true, false, true), l1_freq(l1_freq_v), l1_prefix("mem"),
      l2_comp(l1_d, int_shift), l2_pfp(l2_d_v, l2_comp, l2_p_v, l2_freq_v, l2_w, int_shift),
      l2_pfp_v_table(l2_pfp.dict.alphabet_size)
    { init_v_table(); }
    
    rpfbwt_algo(const std::string& l1_prefix, std::size_t l1_w, std::size_t l2_w)
    : l1_d(l1_prefix, l1_w, l1_d_comp, true, true, true, true, true, true), l1_prefix(l1_prefix),
      l2_comp(l1_d, int_shift),
      l2_pfp(l1_prefix + ".parse", l2_w, l2_comp, int_shift), l2_pfp_v_table(l2_pfp.dict.alphabet_size)
    {
        size_t d1_words; uint_t * occ;
        pfpds::read_file<uint_t> (std::string(l1_prefix + ".occ").c_str(), occ, d1_words);
        l1_freq.push_back(0);
        l1_freq.insert(l1_freq.end(),occ, occ + d1_words);
        delete[] occ;
        
        init_v_table();
    }
    
    std::vector<dict_l1_data_type> l1_bwt(bool out_vector = false)
    {
        std::vector<dict_l1_data_type> out;
        
        size_t i = 1; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
        size_t j = 0;
    
        size_t l_left  = 0;
        size_t l_right = 0;
        std::size_t easy_chars = 0;
        std::size_t hard_easy_chars = 0;
        std::size_t hard_hard_chars = 0;
        std::size_t row = 0;
        while (i < l1_d.saD.size())
        {
            size_t left = i;
        
            auto sn = l1_d.saD[i];
            // Check if the suffix has length at least w and is not the complete phrase.
            auto phrase = l1_d.daD[i] + 1;
            assert(phrase > 0 && phrase < l1_freq.size()); // + 1 because daD is 0-based
            size_t suffix_length = l1_d.select_b_d(l1_d.rank_b_d(sn + 1) + 1) - sn - 1;
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
                
                j += l1_freq[phrase] - 1; // the next bits are 0s
                i++;
            
                l_right += l1_freq[phrase] - 1;
            
                if (i < l1_d.saD.size())
                {
                    auto new_sn = l1_d.saD[i];
                    auto new_phrase = l1_d.daD[i] + 1;
                    assert(new_phrase > 0 && new_phrase < l1_freq.size()); // + 1 because daD is 0-based
                    size_t new_suffix_length = l1_d.select_b_d(l1_d.rank_b_d(new_sn + 1) + 1) - new_sn - 1;
                
                    while (i < l1_d.saD.size() && (l1_d.lcpD[i] >= suffix_length) && (suffix_length == new_suffix_length))
                    {
                        chars.insert(l1_d.d[((new_sn + l1_d.d.size() - 1) % l1_d.d.size())]);
                        pids.insert(new_phrase);
                        j += l1_freq[new_phrase];
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
                            hard_easy_chars += freq;
                        }
                        else
                        {
                            std::size_t curr_l2_row = from_same_l2_suffix[0].first;
                            auto l2_M_entry = l2_pfp.M[curr_l2_row];
                            
                            auto points = l2_pfp.w_wt.range_search_2d(l2_M_entry.left, l2_M_entry.right - 1, l2_pfp.w_wt.size());
                            std::sort(points.begin(), points.end());
                            
                            std::vector<parse_int_type> test_pids;
                            std::vector<uint8_t> test_cs;
                            for (auto& point : points)
                            {
                                std::size_t colex_id = point.second - 1;
                                std::size_t l2_pid = l2_pfp.dict.inv_colex_id[colex_id];
                                parse_int_type l1_pid = l2_pfp.dict.d[l2_pfp.dict.select_b_d(l2_pid + 1) - (l2_M_entry.len + 2)];
                                
                                // check if l1_pid is among the ones we are looking for
                                if (l1_pid >= l2_pfp.shift) { l1_pid -= l2_pfp.shift; }
                                if (pids.contains(l1_pid))
                                {
                                    dict_l1_data_type c = l1_d.d[l1_d.select_b_d(l1_pid + 1) - (suffix_length + 2)];
                                    test_pids.push_back(l1_pid);
                                    test_cs.push_back(c);
                                }
                            }
                            
                            // hard-hard suffix, we hit a row with more than one entry
                            for (auto& pq_entry : from_same_l2_suffix)
                            {
                                uint_t freq = v[pq_entry.second.first].get()[pq_entry.second.second].second;
                                if (out_vector) { out.insert(out.end(), freq, 'H'); }
                                hard_hard_chars += freq;
                            }

                        }
                    }

                }
                
                l_left = l_right + 1;
                l_right = l_left;
                row++;
            }
        }
    
        spdlog::info("Easy: {} Hard-Easy: {} Hard-Hard: {}", easy_chars, hard_easy_chars, hard_hard_chars);
        std::ofstream out_stats(this->l1_prefix + ".stats.csv");
        out_stats << "easy,hard-easy,hard-hard" << std::endl;
        out_stats << easy_chars << "," << hard_easy_chars << ","<< hard_hard_chars << std::endl;
        return out;
    }
    
};

}


#endif //rpfbwt_algorithm_hpp
