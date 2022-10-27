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

private:

    pfpds::dictionary<dict_l1_data_type> l1_d;
    std::vector<uint_t> l1_freq; // here occ has the same size as the integers used for gsacak.
    
    pfpds::pf_parsing<parse_int_type> l2_pfp;

public:
    
    const std::size_t int_shift = 10;
    
    rpfbwt_algo(std::vector<dict_l1_data_type>& l1_d_v,
                std::vector<uint_t>& l1_freq_v,
                std::size_t l1_w,
                std::vector<parse_int_type>& l2_d_v,
                std::vector<parse_int_type>& l2_p_v,
                std::vector<uint_t>& l2_freq_v,
                std::size_t l2_w)
    : l1_d(l1_d_v, l1_w), l1_freq(l1_freq_v), l2_pfp(l2_d_v, l2_p_v, l2_freq_v, l2_w)
    {
    
    }
    
    rpfbwt_algo(const std::string& l1_prefix, std::size_t l1_w, std::size_t l2_w)
    : l1_d(l1_prefix, l1_w), l2_pfp(l1_prefix + ".parse", l2_w)
    {
        size_t d1_words; uint_t * occ;
        pfpds::read_file<uint_t> (std::string(l1_prefix + ".occ").c_str(), occ, d1_words);
        l1_freq.insert(l1_freq.end(),occ, occ + d1_words);
    }

    std::vector<dict_l1_data_type> get_easy_and_semi_hard_chars()
    {
        std::vector<dict_l1_data_type> out;
        
        size_t i = 1; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
        size_t j = 0;
    
        size_t l_left  = 0;
        size_t l_right = 0;
        std::size_t easy_chars = 0;
        std::size_t semi_hard_chars = 0;
        std::size_t hard_chars = 0;
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
                std::set<uint8_t> chars;
                chars.insert(l1_d.d[((sn + l1_d.d.size() - 1) % l1_d.d.size())]);
            
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
            
                if (chars.size() == 1)
                {
                    // easy syffixes
                    out.insert(out.end(), (l_right - l_left) + 1, *(chars.begin()));
                    easy_chars += (l_right - l_left) + 1;
                }
                else
                {
                    // check for semi-hard.
                    hard_chars += (l_right - l_left) + 1;
                }
            
                l_left = l_right + 1;
                l_right = l_left;
            }
        }
    
        spdlog::info("Easy: {} Hard {}", easy_chars, hard_chars);
        return out;
    }

    void get_hard_characters()
    {
        // Iterate over the BWT of P1
        for(uint_t i = 0; i < l2_pfp.n; i++)
        {

        }
    }
    
};

}


#endif //rpfbwt_algorithm_hpp
