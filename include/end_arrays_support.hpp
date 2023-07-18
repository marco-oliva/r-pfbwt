//
// end_arrays_support.hpp
//

#ifndef end_arrays_support_hpp
#define end_arrays_support_hpp

namespace rpfbwt
{

template <typename dict_l1_data_type, typename l2_colex_comp, typename parse_int_type = uint32_t>
class end_arrays_support
{
private:

    // level 1 dictionary
    const pfpds::dictionary<dict_l1_data_type>& l1_d;

    // level 2 dictionary and parse
    const pfpds::dictionary<parse_int_type, l2_colex_comp>& l2_d;
    const pfpds::parse& l2_p;

    // E arrays
    std::vector<std::vector<pfpds::long_type>> E_arrays;

    // level 1 frequencies
    std::vector<pfpds::long_type> l1_frequencies;

    pfpds::long_type int_shift;
    pfpds::long_type l1_lenght;

public:

    end_arrays_support() = delete;

    end_arrays_support(const pfpds::dictionary<dict_l1_data_type>& l1_d, const pfpds::dictionary<parse_int_type, l2_colex_comp>& l2_d, const pfpds::parse& l2_p, pfpds::long_type is = 10)
    : l1_d(l1_d), l2_d(l2_d), l2_p(l2_p), int_shift(is)
    {
        // check for needed data structures
        assert(l2_p.saP_flag);

        build();
    }

    const std::vector<pfpds::long_type>& operator[](std::size_t i) const { return  E_arrays[i]; }
    pfpds::long_type l1_length() const { return this->l1_lenght; }
    pfpds::long_type l1_freq(std::size_t i) const { return this->l1_frequencies[i]; }

private:

    void build()
    {
        spdlog::info("Creating E arrays");

        l1_frequencies.resize(l1_d.n_phrases() + 1);

        E_arrays.resize(l2_d.n_phrases());
        pfpds::long_type end_pos_from_start = 0;

        // read P2 and expand each phrase to get the end position in the text order
        std::vector<pfpds::long_type> end_positions;
        for (pfpds::long_type i = 0; i < l2_p.p.size() - 1; i++)
        {
            pfpds::long_type l2_pid = l2_p.p[i];
            pfpds::long_type l2_phrase_start = l2_d.select_b_d(l2_pid);
            pfpds::long_type l2_phrase_end = l2_phrase_start + l2_d.lengths[l2_pid - 1] - 1;
            if (l2_phrase_start < l2_d.w) { l2_phrase_start = l2_d.w; }

            for (pfpds::long_type j = l2_phrase_start; j <= (l2_phrase_end - l2_d.w); j++)
            {
                pfpds::long_type l1_pid = l2_d.d[j] - int_shift;
                pfpds::long_type l1_phrase_length = l1_d.lengths[l1_pid - 1];

                end_pos_from_start += (l1_phrase_length - l1_d.w);

                l1_frequencies[l1_pid] += 1; // compute l1 frequencies
            }
            end_positions.push_back(end_pos_from_start);
        }

        // now reorder according to P2's SA
        for (pfpds::long_type i = 0; i < l2_p.saP.size(); i++)
        {
            pfpds::long_type sa_value = l2_p.saP[i];

            if (sa_value == 0) { continue; } // { E_arrays[l2_p.p[l2_p.p.size() - 2]].push_back(end_positions.back()); }
            else
            {
                E_arrays[l2_p.p[sa_value - 1] - 1].push_back(end_positions[sa_value - 1]);
            }
        }

        this->l1_lenght = end_pos_from_start;
    }
};

}

#endif // end_arrays_support_hpp
