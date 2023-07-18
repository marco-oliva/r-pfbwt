//
// l2_v_table_support.hpp
//

#ifndef l2_v_table_support_hpp
#define l2_v_table_support_hpp


namespace rpfbwt
{

template<typename dict_data_type, typename colex_comparator_type = std::less<dict_data_type>, class wt_t = pfpds::pfp_wt_custom>
class v_table_support
{
private:

    typedef uint32_t parse_int_type;

    // PFP
    const pfpds::pf_parsing<dict_data_type, colex_comparator_type, wt_t>& pfp;

    // V table
    typedef typename std::vector<std::pair<pfpds::long_type, pfpds::long_type>> table_entry_type;

    std::vector<table_entry_type> pfp_v_table; // (row in which that char appears, number of times per row)

public:

    v_table_support() = delete;

    v_table_support(const pfpds::pf_parsing<dict_data_type, colex_comparator_type, wt_t>& pfp_)
    : pfp(pfp_), pfp_v_table(pfp.dict.alphabet_size)
    {
        // check for needed data structures
        assert(pfp.dict.colex_id_flag);

        build();
    }

    const table_entry_type& operator[](std::size_t i) const { return  pfp_v_table[i]; }

private:

    void build()
    {
        spdlog::info("Building V table");

        for (pfpds::long_type row = 0; row < pfp.M.size(); row++)
        {
            const auto& m = pfp.M[row];

            std::vector<std::pair<uint32_t, pfpds::long_type>> phrase_counts;
            for (pfpds::long_type r = m.left; r <= m.right; r++)
            {
                auto phrase = pfp.dict.colex_id[r];
                pfpds::long_type phrase_start = pfp.dict.select_b_d(phrase + 1);
                pfpds::long_type phrase_length =  pfp.dict.lengths[phrase];
                parse_int_type c = pfp.dict.d[phrase_start + (phrase_length - m.len - 1)];

                if (phrase_counts.empty() or phrase_counts.back().first != c)
                {
                    phrase_counts.push_back(std::make_pair(c, pfp.freq[phrase + 1]));
                }
                else
                {
                    phrase_counts.back().second += pfp.freq[phrase + 1];// 1;
                }
            }

            // update
            for (auto const& c_count : phrase_counts)
            {
                std::pair<pfpds::long_type, pfpds::long_type> v_element = std::make_pair(row, c_count.second);
                pfp_v_table[c_count.first].emplace_back(v_element);
            }
        }
    }
};

}

#endif // l2_v_table_support_hpp