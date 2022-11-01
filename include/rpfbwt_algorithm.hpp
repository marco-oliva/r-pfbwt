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

public: // TODO: back to private

    pfpds::dictionary<dict_l1_data_type> l1_d;
    std::vector<uint_t> l1_freq; // here occ has the same size as the integers used for gsacak.
    
    pfpds::pf_parsing<parse_int_type> l2_pfp;
    
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
    
    const std::size_t int_shift = 10;
    
    rpfbwt_algo(std::vector<dict_l1_data_type>& l1_d_v,
                std::vector<uint_t>& l1_freq_v,
                std::size_t l1_w,
                std::vector<parse_int_type>& l2_d_v,
                std::vector<parse_int_type>& l2_p_v,
                std::vector<uint_t>& l2_freq_v,
                std::size_t l2_w)
    : l1_d(l1_d_v, l1_w), l1_freq(l1_freq_v), l2_pfp(l2_d_v, l2_p_v, l2_freq_v, l2_w), l2_pfp_v_table(l2_pfp.dict.alphabet_size)
    { init_v_table(); }
    
    rpfbwt_algo(const std::string& l1_prefix, std::size_t l1_w, std::size_t l2_w)
    : l1_d(l1_prefix, l1_w), l2_pfp(l1_prefix + ".parse", l2_w), l2_pfp_v_table(l2_pfp.dict.alphabet_size)
    {
        size_t d1_words; uint_t * occ;
        pfpds::read_file<uint_t> (std::string(l1_prefix + ".occ").c_str(), occ, d1_words);
        l1_freq.push_back(0);
        l1_freq.insert(l1_freq.end(),occ, occ + d1_words);
        delete[] occ;
        
        init_v_table();
    }
    
    std::vector<dict_l1_data_type> l1_bwt()
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
                
                // use the RMQ data structure to find how many of the following suffixes are the same except for the terminator (so they're the same suffix but in different phrases)
                // use the document array and the table of phrase frequencies to find the phrases frequencies and sum them up
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
                    out.insert(out.end(), (l_right - l_left) + 1, *(chars.begin()));
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
                    
                    // get all chars from this row entry
                    std::vector<pq_t> from_same_l2_suffix;
                    auto first_suffix = pq.top();
                    while ((not pq.empty()) and (pq.top().first == first_suffix.first))
                    {
                        auto curr = pq.top(); pq.pop();
                        from_same_l2_suffix.push_back(curr);
                        
                        std::size_t arr_i_c = curr.second.first;  // ith array
                        std::size_t arr_x_c = curr.second.second; // index in i-th array
                        if (arr_i_c + 1 < v[arr_i_c].get().size())
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
                        out.insert(out.end(), freq, c);
                        hard_easy_chars += freq;
                    }
                    else
                    {
                        // hard-hard suffix
                        for (auto& pq_entry : from_same_l2_suffix)
                        {
                            uint_t freq = v[pq_entry.second.first].get()[pq_entry.second.second].second;
                            out.insert(out.end(), freq, 'H');
                            hard_hard_chars += freq;
                        }
                        
                    }
                }
                
                l_left = l_right + 1;
                l_right = l_left;
                row++;
            }
        }
    
        spdlog::info("Easy: {} Hard {}", easy_chars, hard_hard_chars);
        return out;
    }
    
};

//
//// A pair of pairs, first element is going to
//// store value, second element index of array
//// and third element index in the array.
//    typedef pair<int, pair<int, int> > ppi;
//
//// This function takes an array of arrays as an
//// argument and all arrays are assumed to be
//// sorted. It merges them together and prints
//// the final sorted output.
//    vector<int> mergeKArrays(vector<vector<int> > arr)
//    {
//        vector<int> output;
//
//        // Create a min heap with k heap nodes. Every
//        // heap node has first element of an array
//        priority_queue<ppi, vector<ppi>, greater<ppi> > pq;
//
//        for (int i = 0; i < arr.size(); i++)
//            pq.push({ arr[i][0], { i, 0 } });
//
//        // Now one by one get the minimum element
//        // from min heap and replace it with next
//        // element of its array
//        while (pq.empty() == false) {
//            ppi curr = pq.top();
//            pq.pop();
//
//            // i ==> Array Number
//            // j ==> Index in the array number
//            int i = curr.second.first;
//            int j = curr.second.second;
//
//            output.push_back(curr.first);
//
//            // The next element belongs to same array as
//            // current.
//            if (j + 1 < arr[i].size())
//                pq.push({ arr[i][j + 1], { i, j + 1 } });
//        }
//
//        return output;
//    }
//
//// Driver program to test above functions
//    int main()
//    {
//        // Change n at the top to change number
//        // of elements in an array
//        vector<vector<int> > arr{ { 2, 6, 12 },
//                                  { 1, 9 },
//                                  { 23, 34, 90, 2000 } };
//
//        vector<int> output = mergeKArrays(arr);
//
//        cout << "Merged array is " << endl;
//        for (auto x : output)
//            cout << x << " ";
//
//        return 0;
//    }

}


#endif //rpfbwt_algorithm_hpp
