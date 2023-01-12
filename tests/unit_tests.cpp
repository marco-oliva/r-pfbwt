//
//  unit_tests.cpp
//
//  Copyright 2022 Marco Oliva. All rights reserved.
//

#include <iostream>

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>

#include <spdlog/spdlog.h>


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

#include <kseq.h>
#include <stdio.h>
KSEQ_INIT(int, read)

template <class char_type>
void read_fasta_file(const char *filename, std::vector<char_type>& v, std::size_t w = 10){
    FILE* fd;
    
    if ((fd = fopen(filename, "r")) == nullptr)
        spdlog::error("open() file " + std::string(filename) + " failed" );
    
    v.clear();
    v.push_back('\2'); // Dollar
    
    kseq_t *record;
    record = kseq_init(fileno(fd));
    while(kseq_read(record) >= 0)
    {
        // Add sequence
        for (std::size_t i = 0; i < record->seq.l; i++) { v.push_back(record->seq.s[i]); }
        v.insert(v.end(), w - 1, '\5'); // Dollar Prime
        v.insert(v.end(), 1, '\4'); // Dollar Sequence
    }
    
    kseq_destroy(record);
    fclose(fd);
}

//------------------------------------------------------------------------------

#include <pfp/pfp.hpp>
#include <pfp/sa_support.hpp>

#include <rpfbwt_algorithm.hpp>

TEST_CASE( "pfp<uint8_t> SA from example, no chunks", "PFP on example" )
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
    { 0, 6, 16}, // -
    { 6, 16, 24, 10, 8, 19, 15, 22 , 3, 6, 16 }, // -
    { 6, 16, 24, 12, 23, 14, 21, 2, 9 }, //
    { 2, 9, 16, 18, 4, 19, 14, 17, 2, 9 }, //
    { 2, 9, 11, 4, 19, 13,  6, 16 }, // -
    { 6, 16, 25, 3, 6, 16 }, //
    { 6, 16, 24, 16, 18, 8, 23, 15, 22, 2, 9 }, //
    { 2, 9, 12, 19, 13, 7, 3, 6, 16 }, // -
    { 6, 16, 24, 16, 20, 8, 23, 13, 7, 3, 5, 11, 4, 19, 15, 22, 1 } //
    };
    for (auto& v : dict_l2_prep) { for (auto& e : v) { e += int_shift + 1; } }
    dict_l2_prep[0].insert(dict_l2_prep[0].begin(), w_l2, Dollar);
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

    std::size_t chunks = 5;
    std::size_t threads = 12;
    rpfbwt::rpfbwt_algo<uint8_t> rpfbwt_algo(dict_l1, freq_l1, w_l1, dict_l2, parse_l2, freq_l2, w_l2, chunks);
    rpfbwt_algo.l1_refined_rindex(threads);

    bool all_good = true;
    std::size_t mismatch = 0;

    // Check BWT
    rle::RLEString rle_bwt;
    rle_bwt.load("mem.rlebwt");

    all_good = true; mismatch = 0;
    for (std::size_t i = 0; i < bwt.size(); i++)
    {
        all_good = all_good and (bwt[i] == rle_bwt[i]);
        if (not all_good) { mismatch = i; break; }
    }
    REQUIRE(all_good);

    // Check SA
    std::less<uint8_t> char_comp;
    std::vector<uint32_t> parse_l1_copy(parse_l1);
    parse_l1_copy.push_back(0);
    pfpds::pf_parsing<uint8_t> pfp(dict_l1, char_comp, parse_l1_copy, freq_l1, 2);
    pfpds::pfp_sa_support<uint8_t> sa_support(pfp);

    std::vector<std::size_t> sa_values;
    std::ifstream sa_values_is("mem.ssa");
    std::size_t sa_values_size = 0;
    sa_values_is.read((char*)&sa_values_size, sizeof(std::size_t));
    sa_values.resize(sa_values_size);
    sa_values_is.read((char*)&sa_values[0], sizeof(std::size_t) * sa_values.size());
    
    std::size_t sa_i = 0;
    rle::RLEString::RunIterator rle_iterator(rle_bwt);
    while (not rle_iterator.end())
    {
        rle::RunType run = *rle_iterator;
        std::size_t offset = rle::RunTraits::start(run);
        
        std::size_t from_sa_support = sa_support(offset);
        all_good = all_good and (sa_values[sa_i] == from_sa_support);
        if (not all_good) { mismatch = sa_i; break; }
    
        rle_iterator.operator++(); sa_i++;
    }
    REQUIRE(all_good);
}

TEST_CASE( "Compare r-index with PFP-DS", "PFP on yeast.fasta" )
{
    bool all_good = true;
    std::size_t mismatch = 0;

    std::string yeast_pfp_path = testfiles_dir + "/yeast.fasta";
    std::size_t w_l1 = 10, w_l2 = 5;

    // Build PFP-DS, sa_support and bwt
    std::less<char> char_comp;
    pfpds::pf_parsing<char> pfp(yeast_pfp_path, w_l1, char_comp);
    pfpds::pfp_sa_support<char> sa_support(pfp);

    std::vector<char> pfpds_bwt(pfp.n, '\0');
    for (std::size_t i = 0; i < pfp.n; i++)
    {
        auto sn = (sa_support(i) + pfp.w - 1) % pfp.n;  // suffix number
        auto p_i = pfp.rank_b_p(sn + 1);                // phrase number
        auto id_p_i = pfp.pars.p[p_i - 1];              // phrase_id of the phrase that i belongs to.
        size_t occ_in_p_i_in_D = pfp.dict.select_b_d(id_p_i) + (sn - pfp.select_b_p(p_i));
        auto c = pfp.dict.d[occ_in_p_i_in_D];
        pfpds_bwt[i] = c;
    }
    spdlog::info("Done building test bwt");

    // Build the ri
    std::size_t chunks = 5;
    std::size_t threads = 12;
    rpfbwt::rpfbwt_algo<char> rpfbwt_algo(yeast_pfp_path, w_l1, w_l2, chunks);
    rpfbwt_algo.l1_refined_rindex(threads);

    // Read in the rle bwt
    rle::RLEString rle_bwt;
    rle_bwt.load(yeast_pfp_path + ".rlebwt");

    all_good = true; mismatch = 0;
    for (std::size_t i = 0; i < pfpds_bwt.size(); i++)
    {
        all_good = all_good and (pfpds_bwt[i] == rle_bwt[i]);
        if (not all_good) { mismatch = i; break; }
    }
    REQUIRE(all_good);

    // Check SA
    std::vector<std::size_t> sa_values;
    std::ifstream sa_values_is(yeast_pfp_path + ".ssa");
    std::size_t sa_values_size = 0;
    sa_values_is.read((char*)&sa_values_size, sizeof(std::size_t));
    sa_values.resize(sa_values_size);
    sa_values_is.read((char*)&sa_values[0], sizeof(std::size_t) * sa_values.size());

    std::size_t sa_i = 0;
    rle::RLEString::RunIterator rle_iterator(rle_bwt);
    while (not rle_iterator.end())
    {
        rle::RunType run = *rle_iterator;
        std::size_t offset = rle::RunTraits::start(run);

        std::size_t from_sa_support = sa_support(offset);
        all_good = all_good and (sa_values[sa_i] == from_sa_support);
        if (not all_good) { mismatch = sa_i; break; }

        rle_iterator.operator++(); sa_i++;
    }
    REQUIRE(all_good);
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