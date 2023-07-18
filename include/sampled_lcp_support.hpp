//
// sample_lcp_support.hpp
//

#ifndef sample_lcp_support_hpp
#define sample_lcp_support_hpp

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

namespace rpfbwt
{

template <typename dict_l1_data_type, typename parse_int_type = uint32_t>
class recursive_slcp_support
{
private:

    // level 1 dictionary
    const pfpds::dictionary<dict_l1_data_type>& l1_d;

    // level 2 dictionary and parse
    const pfpds::dictionary<parse_int_type>& l2_d;
    const pfpds::parse& l2_p;

    // data structures
    sdsl::int_vector<0> e_lcp_l2_D; // l2_d LCP expanded, i.e., number of matching meta-characters + matching characters
                                    // in first non-matching phrase
    sdsl::rmq_succinct_sct<> rmq_e_lcp_l2_D;

    sdsl::int_vector<0> e_lcp_l2_P; // l2_P LCP expanded, i.e., size of expanded matching meta-characters + matching characters
                                    // in first non-matching phrase
    sdsl::rmq_succinct_sct<> rmq_e_lcp_l2_P;

    template <class dict_type, class rmq_type>
    inline std::size_t longest_common_phrase_prefix(dict_type& d, rmq_type& rmq, std::size_t a, std::size_t b)
    {
        if(a == 0 || b == 0) { return 0; }

        // Compute the lcp between phrases a and b
        auto a_in_sa = d.isaD[d.select_b_d(a)]; // position of the phrase a in saD
        auto b_in_sa = d.isaD[d.select_b_d(b)]; // position of the phrase b in saD

        auto lcp_left = std::min(a_in_sa, b_in_sa) + 1;
        auto lcp_right = std::max(a_in_sa, b_in_sa);

        std::size_t lcp_a_b_i = rmq_lcp_D(lcp_left, lcp_right);
        return d.lcpD[lcp_a_b_i];
    }

public:

    recursive_slcp_support() = delete;

    recursive_slcp_support(pfpds::dictionary<dict_l1_data_type>& l1_d, pfpds::dictionary<parse_int_type>& l2_d, pfpds::parse& l2_p)
    : l1_d(l1_d), l2_d(l2_d), l2_p(l2_p)
    {
        // check for needed data structures
        assert(l1_d.lcpD_flag);
        assert(l1_d.rmq_lcp_D_flag);
        assert(l1_d.saD_flag);
        assert(l1_d.isaD_flag);
        assert(l2_d.lcpD_flag);
        assert(l2_p.saP_flag);
        assert(l2_p.isaP_flag);

        build();
    }

    std::size_t extended_lcp_l2_D(std::size_t i) const { return e_lcp_l2_D[i]; }
    std::size_t extended_lcp_l2_P(std::size_t i) const { return e_lcp_l2_P[i]; }

    std::size_t rmq_eLCP_l2_D(std::size_t a, std::size_t b)
    {
        std::size_t l = std::min(a , b);
        std::size_t r = std::max(a , b);
        return rmq_e_lcp_l2_D(l , r);
    }

    std::size_t rmq_eLCP_l2_P(std::size_t a, std::size_t b)
    {
        std::size_t l = std::min(a , b);
        std::size_t r = std::max(a , b);
        return rmq_e_lcp_l2_P(l , r);
    }

private:

    void build()
    {
        // compute the expanded eLCP for the dictionary of level 2
        spdlog::info("Using {} bits for expanded LCP of the L2 dictionary", l1_d.lcpD.width());
        e_lcp_l2_D = sdsl::int_vector<>(l2_d.lcpD.size(), 0ULL, l1_d.lcpD.width());

        for (std::size_t i = 1; i < l2_d.lcpD.size(); i++)
        {
            e_lcp_l2_D[i] = 0;
            for (std::size_t l = 0; l < l2_d.lcpD[i]; l++)
            {
                e_lcp_l2_D[i] += l1_d.length_of_phrase(l1_d.d[l1_d.saD[i]] + l) - l1_d.w;
            }

            // get the first non-matching meta-character for the two suffixes
            parse_int_type p_id_1 = l2_d.d[l2_d.saD[i - 1] + 1];
            parse_int_type p_id_2 = l2_d.d[l2_d.saD[i] + 1];

            // if for one of the two suffixes the match is as long as the phrase, then the value of the eLCP is the value of the LCP
            if (p_id_1 == EndOfWord or p_id_2 == EndOfWord) { continue; }

            // otherwise, we need to know how long is the match between p_id_1 and p_id_2
            e_lcp_l2_D[i] += longest_common_phrase_prefix(l2_d, l2_d.rmq_lcp_D, p_id_1, p_id_2);
        }

        // associated rmq structure
        rmq_e_lcp_l2_D = sdsl::rmq_succinct_sct<>(e_lcp_l2_D);

        // now compute the sampled SLCP at the phrase boundaries of level 2
        // the computation is the same as it would be for level 1, I just need to use the rmq_e_lcp_l2_D
        spdlog::info("Using {} bits for expanded LCP of the L2 parse", l2_p.saP.width());
        e_lcp_l2_P = sdsl::int_vector<>(l2_p.saP.size(), 0ULL, l2_p.saP.width());

        std::size_t l = 0;
        std::size_t lt = 0;
        for (std::size_t i = 0; i < l2_p.saP.size(); i++)
        {
            // if i is the last character LCP is not defined
            std::size_t k = l2_p.isaP[i];
            if (k > 0)
            {
                std::size_t j = l2_p.saP[k - 1];
                // find the longest common prefix of the i-th suffix and the j-th suffix.
                while (l2_p.p[i + l] == l2_p.p[j + l])
                {
                    std::size_t l2_p_s = l2_d.select_b_d(l2_p.p[i + l]);
                    for (std::size_t l2_p_it = 0; l2_p_it < l2_d.length_of_phrase(l2_p.p[i + l]) - l2_d.w; l2_p_it++)
                    {
                        lt += l1_d.length_of_phrase(l2_d.d[l2_p_s + l2_p_it]) - l1_d.w;
                        l++;
                    }
                }

                std::size_t lcp_p = longest_common_phrase_prefix(l2_d, rmq_e_lcp_l2_D, l2_p.p[i + l], l2_p.p[j + l]);

                // l stores the length of the longest common prefix between the i-th suffix and the j-th suffix
                e_lcp_l2_P[k] = lt + lcp_p;
                if (l > 0)
                {
                    l--;
                    lt -= l2_d.length_of_phrase(l2_p.p[i]) - l2_d.w; // I have to remove the length of the first matching phrase
                }
            }
        }
        rmq_e_lcp_l2_P = sdsl::rmq_succinct_sct<>(&e_lcp_l2_P);
    }

};

}

#endif //sample_lcp_support_hpp
