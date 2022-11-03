//
//  unit_tests.cpp
//
//  Copyright 2022 Marco Oliva. All rights reserved.
//

#include <iostream>

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>


//------------------------------------------------------------------------------

struct listener : Catch::EventListenerBase
{
  using EventListenerBase::EventListenerBase;
  
  virtual void testCaseStarting(Catch::TestCaseInfo const& testInfo) override
  {
      std::cout << testInfo.tagsAsString() << " " << testInfo.name << std::endl;
  }
};
CATCH_REGISTER_LISTENER(listener)

//------------------------------------------------------------------------------

std::string testfiles_dir = "../tests/files";

//------------------------------------------------------------------------------

#include <pfp/pfp.hpp>
#include <rpfbwt_algorithm.hpp>

TEST_CASE( "pfp<uint8_t> from example", "PFP on example" )
{
    size_t w_l1 = 2;
    
    std::vector<char> text = {'A','C','G','T','T','C','G','C','A','A','C','T','A','G','T','C','C','G','G','G','A','G','T','T','A','C','#',
                              'A','C','G','T','T','C','G','G','A','A','T','A','G','T','T','C','C','G','G','G','A','G','G','T','T','A','A','C','#',
                              'A','A','C','T','T','C','G','T','C','A','C','C','T','A','G','T','C','C','G','G','G','A','G','G','T','A','A','C','#',
                              'A','A','C','T','T','C','G','C','A','C','C','T','A','G','T','C','C','G','G','G','A','C','G','T','T','T','A','C','#',
                              'A','C','G','T','T','C','G','T','C','A','C','T','A','G','T','T','C','C','G','G','G','A','G','T','T','A','C','#',
                              'A','A','C','T','T','C','G','G','A','A','T','A','G','T','C','C','G','G','G','A','C','T','A','A','C','#',
                              'A','C','G','T','T','C','G','T','G','A','C','T','A','G','T','T','C','C','G','G','G','A','C','T','A','A','C','#',
                              'A','C','C','T','T','C','G','C','A','C','C','T','A','G','T','C','C','G','G','G','A','G','T','T','A','C','#','#'
    };
    
    std::vector<std::string> dict_l1_prep
    {
    "##AC", // 0
    "AC##", // 1
    "AC#AAC", // 2
    "AC#AC", // 3
    "ACCTAGT", // 4
    "ACCTTCG", // 5
    "ACG", // 6
    "ACTAAC", // 7
    "ACTAGT", // 8
    "ACTTCG", // 9
    "CGCAAC", // 10
    "CGCAC", // 11
    "CGGAATAGT", // 12
    "CGGGAC", // 13
    "CGGGAGGT", // 14
    "CGGGAGT", // 15
    "CGT", // 16
    "GTAAC", // 17
    "GTCAC", // 18
    "GTCCG", // 19
    "GTGAC", // 20
    "GTTAAC", // 21
    "GTTAC", // 22
    "GTTCCG", // 23
    "GTTCG", // 24
    "GTTTAC", // 25
    };
    
    std::vector<uint8_t> dict_l1;
    for (auto& phrase : dict_l1_prep)
    {
        for (auto& c : phrase)
        { dict_l1.push_back(c); }
        dict_l1.push_back(EndOfWord);
    }
    dict_l1.push_back(EndOfDict);
    
    std::vector<uint32_t> parse_l1
    {
        0, 6, 16, 24, 10, 8, 19, 15, 22, 3, 6, 16, 24, 12, 23, 14, 21, 2, 9, 16, 18, 4, 19,
        14, 17, 2, 9, 11, 4, 19, 13, 6, 16, 25, 3, 6, 16, 24, 16, 18, 8, 23, 15, 22, 2, 9,
        12, 19, 13, 7, 3, 6, 16, 24, 16, 20, 8, 23, 13, 7, 3, 5, 11, 4, 19, 15, 22, 1
    };
    for (auto& p_id : parse_l1) { p_id += 1; }
    std::vector<uint_t> freq_l1(dict_l1_prep.size() + 1, 0);
    for (auto& p_id : parse_l1) { freq_l1[p_id] += 1; }
    
    
    // L2, TS: 6-16, 2,9
    std::size_t w_l2 = 2;
    uint32_t int_shift = 10;
    std::vector<std::vector<uint32_t>> dict_l2_prep =
    {
    { 0, 6, 16 },
    { 6, 16, 24, 10, 8, 19, 15, 22 , 3, 6, 16 },
    { 6, 16, 24, 12, 23, 14, 21, 2, 9 },
    { 2, 9, 16, 18, 4, 19, 14, 17, 2, 9 },
    { 2, 9, 11, 4, 19, 13,  6, 16 },
    { 6, 16, 25, 3, 6, 16 },
    { 6, 16, 24, 16, 18, 8, 23, 15, 22, 2, 9 },
    { 2, 9, 12, 19, 13, 7, 3, 6, 16 },
    { 6, 16, 24, 16, 20, 8, 23, 13, 7, 3, 5, 11, 4, 19, 15, 22, 1 },
    };
    for (auto& v : dict_l2_prep) { for (auto& e : v) { e += int_shift + 1; } }
    dict_l2_prep[0].insert(dict_l2_prep[0].begin(), 1, Dollar);
    dict_l2_prep[8].insert(dict_l2_prep[8].end(), w_l2, Dollar);
    std::sort(dict_l2_prep.begin(), dict_l2_prep.end());
    
    std::vector<uint32_t> dict_l2;
    for (auto& phrase : dict_l2_prep)
    {
        for (auto& c : phrase)
        { dict_l2.push_back(c); }
        dict_l2.push_back(EndOfWord);
    }
    dict_l2.push_back(EndOfDict);
    
    std::vector<uint32_t> parse_l2 = { 1, 5, 6, 4, 2, 9, 7, 3, 8, 0};
    
    std::vector<uint_t> freq_l2(dict_l2_prep.size() + 1, 0);
    for (auto& p_id : parse_l2) { freq_l2[p_id] += 1; }
    freq_l2[0] = 0;
    
    std::string bwt = "CCCCC#CCCTTTTC###GGTATAATTACCC#####GGGAGCAAAGGTTTTTGGGTTTAAAAAAAAAAGGTGTTTTTTTTTAAA"
                      "ATTTTTCCCCCCCCTTTAAAAAAACCCAAAACAAGGGGGTGGGGGCCCCCGGGGGGGGCCCCCCCCAAGCCAAAAACGAAAAA"
                      "ACCCCCGTCCTTTTCACCCCACGGGTGGTGTGTTTTTTTTGGGGGTGGGGCCCGCGGG";
    
    rpfbwt::rpfbwt_algo<uint8_t> rpfbwt_algo(dict_l1, freq_l1, w_l1, dict_l2, parse_l2, freq_l2, w_l2);
    std::vector<uint8_t> easy_easy_and_hard_easy_chars = rpfbwt_algo.l1_bwt(true);

    REQUIRE(not easy_easy_and_hard_easy_chars.empty());
    for (std::size_t bwt_i = 0; bwt_i < bwt.size(); bwt_i++)
    {
        bool all_good = (easy_easy_and_hard_easy_chars[bwt_i] == 'H') or (easy_easy_and_hard_easy_chars[bwt_i] == bwt[bwt_i]);
        REQUIRE(all_good);
    }

    size_t i = 1; // This should be safe since the first entry of sa is always the dollarsign used to compute the sa
    size_t j = 0;
    while (i < rpfbwt_algo.l2_pfp.dict.saD.size())
    {
        auto sn = rpfbwt_algo.l2_pfp.dict.saD[i];
        // Check if the suffix has length at least w and is not the complete phrase.
        auto phrase = rpfbwt_algo.l2_pfp.dict.daD[i] + 1;
        assert(phrase > 0 && phrase < rpfbwt_algo.l2_pfp.freq.size()); // + 1 because daD is 0-based
        size_t suffix_length = rpfbwt_algo.l2_pfp.dict.select_b_d(rpfbwt_algo.l2_pfp.dict.rank_b_d(sn + 1) + 1) - sn - 1;
        if (rpfbwt_algo.l2_pfp.dict.b_d[sn] || suffix_length < rpfbwt_algo.l2_pfp.dict.w) // skip full phrases or suffixes shorter than w
        {
            ++i; // Skip
        }
        else
        {
            uint32_t suff = rpfbwt_algo.l2_pfp.dict.d[((sn + rpfbwt_algo.l2_pfp.dict.d.size() - 1) % rpfbwt_algo.l2_pfp.dict.d.size())];
            if (suff > int_shift) { suff -= int_shift; }
            std::cout <<  suff << "\t-\t";
            std::size_t sit = 0;
            while (sit < suffix_length)
            {
                uint32_t c = rpfbwt_algo.l2_pfp.dict.d[sn + sit];
                if (c > int_shift) { c -= int_shift; }
                std::cout << c << " ";
                sit++;
            }
            std::cout << " - " << rpfbwt_algo.l2_pfp.freq[phrase];
            std::cout << std::endl;
    
            ++i;
        }
    }
}

//------------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    Catch::Session session;
    
    using namespace Catch::Clara;
    
    auto cli = session.cli() |
    Opt( testfiles_dir, "dir" ) ["--test-dir"] ("specify the directory containing the test dna sequences files");
    
    session.cli(cli);
    
    int returnCode = session.applyCommandLine(argc, argv);
    
    if( returnCode != 0 ) return returnCode;
    
    session.run();
}