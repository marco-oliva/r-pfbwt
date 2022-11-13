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
    
    kseq_t *record;
    record = kseq_init(fileno(fd));
    while(kseq_read(record) >= 0)
    {
        // Add sequence
        for (std::size_t i = 0; i < record->seq.l; i++) { v.push_back(record->seq.s[i]); }
        v.insert(v.end(), w - 1, '\5'); // Dollar Prime
        v.insert(v.end(), 1, '\4'); // Dollar Sequence
    }
    v.insert(v.end(), w - 1, '\2'); // Dollar
    
    kseq_destroy(record);
    fclose(fd);
}

//------------------------------------------------------------------------------

#include <pfp/pfp.hpp>
#include <rpfbwt_algorithm.hpp>

//TEST_CASE( "pfp<uint8_t> from example, no chunks", "PFP on example" )
//{
//    size_t w_l1 = 2;
//
//    std::vector<char> text = {'A','C','G','T','T','C','G','C','A','A','C','T','A','G','T','C','C','G','G','G','A','G','T','T','A','C','#',
//                              'A','C','G','T','T','C','G','G','A','A','T','A','G','T','T','C','C','G','G','G','A','G','G','T','T','A','A','C','#',
//                              'A','A','C','T','T','C','G','T','C','A','C','C','T','A','G','T','C','C','G','G','G','A','G','G','T','A','A','C','#',
//                              'A','A','C','T','T','C','G','C','A','C','C','T','A','G','T','C','C','G','G','G','A','C','G','T','T','T','A','C','#',
//                              'A','C','G','T','T','C','G','T','C','A','C','T','A','G','T','T','C','C','G','G','G','A','G','T','T','A','C','#',
//                              'A','A','C','T','T','C','G','G','A','A','T','A','G','T','C','C','G','G','G','A','C','T','A','A','C','#',
//                              'A','C','G','T','T','C','G','T','G','A','C','T','A','G','T','T','C','C','G','G','G','A','C','T','A','A','C','#',
//                              'A','C','C','T','T','C','G','C','A','C','C','T','A','G','T','C','C','G','G','G','A','G','T','T','A','C','#','#'
//    };
//
//    std::vector<std::string> dict_l1_prep
//    {
//    "##AC", // 0
//    "AC##", // 1
//    "AC#AAC", // 2
//    "AC#AC", // 3
//    "ACCTAGT", // 4
//    "ACCTTCG", // 5
//    "ACG", // 6
//    "ACTAAC", // 7
//    "ACTAGT", // 8
//    "ACTTCG", // 9
//    "CGCAAC", // 10
//    "CGCAC", // 11
//    "CGGAATAGT", // 12
//    "CGGGAC", // 13
//    "CGGGAGGT", // 14
//    "CGGGAGT", // 15
//    "CGT", // 16
//    "GTAAC", // 17
//    "GTCAC", // 18
//    "GTCCG", // 19
//    "GTGAC", // 20
//    "GTTAAC", // 21
//    "GTTAC", // 22
//    "GTTCCG", // 23
//    "GTTCG", // 24
//    "GTTTAC", // 25
//    };
//
//    std::vector<uint8_t> dict_l1;
//    for (auto& phrase : dict_l1_prep)
//    {
//        for (auto& c : phrase)
//        { dict_l1.push_back(c); }
//        dict_l1.push_back(EndOfWord);
//    }
//    dict_l1.push_back(EndOfDict);
//
//    std::vector<uint32_t> parse_l1
//    {
//        0, 6, 16, 24, 10, 8, 19, 15, 22, 3, 6, 16, 24, 12, 23, 14, 21, 2, 9, 16, 18, 4, 19,
//        14, 17, 2, 9, 11, 4, 19, 13, 6, 16, 25, 3, 6, 16, 24, 16, 18, 8, 23, 15, 22, 2, 9,
//        12, 19, 13, 7, 3, 6, 16, 24, 16, 20, 8, 23, 13, 7, 3, 5, 11, 4, 19, 15, 22, 1
//    };
//    for (auto& p_id : parse_l1) { p_id += 1; }
//    std::vector<uint_t> freq_l1(dict_l1_prep.size() + 1, 0);
//    for (auto& p_id : parse_l1) { freq_l1[p_id] += 1; }
//
//
//    // L2, TS: 6-16, 2,9
//    std::size_t w_l2 = 2;
//    uint32_t int_shift = 10;
//    std::vector<std::vector<uint32_t>> dict_l2_prep =
//    {
//    { 0, 6, 16}, // -
//    { 6, 16, 24, 10, 8, 19, 15, 22 , 3, 6, 16 }, // -
//    { 6, 16, 24, 12, 23, 14, 21, 2, 9 }, //
//    { 2, 9, 16, 18, 4, 19, 14, 17, 2, 9 }, //
//    { 2, 9, 11, 4, 19, 13,  6, 16 }, // -
//    { 6, 16, 25, 3, 6, 16 }, //
//    { 6, 16, 24, 16, 18, 8, 23, 15, 22, 2, 9 }, //
//    { 2, 9, 12, 19, 13, 7, 3, 6, 16 }, // -
//    { 6, 16, 24, 16, 20, 8, 23, 13, 7, 3, 5, 11, 4, 19, 15, 22, 1 } //
//    };
//    for (auto& v : dict_l2_prep) { for (auto& e : v) { e += int_shift + 1; } }
//    dict_l2_prep[0].insert(dict_l2_prep[0].begin(), 1, Dollar);
//    dict_l2_prep[8].insert(dict_l2_prep[8].end(), w_l2, Dollar);
//    std::sort(dict_l2_prep.begin(), dict_l2_prep.end());
//
//    std::vector<uint32_t> dict_l2;
//    for (auto& phrase : dict_l2_prep)
//    {
//        for (auto& c : phrase)
//        { dict_l2.push_back(c); }
//        dict_l2.push_back(EndOfWord);
//    }
//    dict_l2.push_back(EndOfDict);
//
//    std::vector<uint32_t> parse_l2 = { 1, 5, 6, 4, 2, 9, 7, 3, 8, 0};
//
//    std::vector<uint_t> freq_l2(dict_l2_prep.size() + 1, 0);
//    for (auto& p_id : parse_l2) { freq_l2[p_id] += 1; }
//    freq_l2[0] = 0;
//
//    std::string bwt = "CCCCC#CCCTTTTC###GGTATAATTACCC#####GGGAGCAAAGGTTTTTGGGTTTAAAAAAAAAAGGTGTTTTTTTTTAAA"
//                      "ATTTTTCCCCCCCCTTTAAAAAAACCCAAAACAAGGGGGTGGGGGCCCCCGGGGGGGGCCCCCCCCAAGCCAAAAACGAAAAA"
//                      "ACCCCCGTCCTTTTCACCCCACGGGTGGTGTGTTTTTTTTGGGGGTGGGGCCCGCGGG";
//
//    rpfbwt::rpfbwt_algo<uint8_t> rpfbwt_algo(dict_l1, freq_l1, w_l1, dict_l2, parse_l2, freq_l2, w_l2);
//    std::vector<uint8_t> easy_easy_and_hard_easy_chars = rpfbwt_algo.l1_bwt(true);
//
//    REQUIRE(not easy_easy_and_hard_easy_chars.empty());
//    bool all_good = true;
//    for (std::size_t bwt_i = 0; bwt_i < bwt.size(); bwt_i++)
//    {
//        all_good = all_good and (easy_easy_and_hard_easy_chars[bwt_i] == bwt[bwt_i]);
//    }
//
//    REQUIRE(all_good);
//}


#include <pfp/sa_support.hpp>
#include <pfp/ra_support.hpp>

TEST_CASE( "RLBWT for yeast", "PFP on yeast.fasta" )
{
    bool all_good = true;
    std::size_t miss_match = 0;

    // Build the rle bwt
    std::string yeast_pfp_path = testfiles_dir + "/yeast.fasta";
    std::size_t blocks = 5;
    std::size_t w_l1 = 10, w_l2 = 5;
    rpfbwt::rpfbwt_algo<char> rpfbwt_algo_1_chunk(yeast_pfp_path, w_l1, w_l2, 1);
    rpfbwt::rpfbwt_algo<char> rpfbwt_algo_n_chunks(yeast_pfp_path, w_l1, w_l2, blocks);
    
    std::vector<char> v_bwt_1c = rpfbwt_algo_1_chunk.l1_bwt(true);
    std::vector<char> v_bwt_nc = rpfbwt_algo_n_chunks.l1_bwt(true);
    
    REQUIRE(v_bwt_1c.size() == v_bwt_nc.size());
    all_good = true; miss_match = 0;
    for (std::size_t i = 0; i < v_bwt_1c.size(); i++)
    {
        all_good = all_good and (v_bwt_1c[i] == v_bwt_nc[i]);
        if (not all_good) { miss_match = i; break; }
    }
    REQUIRE(all_good);
    
    // Check parallel construction of rlebwt on disk
    rpfbwt::rpfbwt_algo<char> rpfbwt_algo_n_chunks_p(yeast_pfp_path, w_l1, w_l2, blocks);
    rpfbwt_algo_n_chunks_p.l1_bwt_parallel();

    // Read in rle bwt
    std::string rle_path = yeast_pfp_path + ".rlebwt";
    rle::RLEString rle_bwt; rle_bwt.load(rle_path);
    
    REQUIRE(v_bwt_1c.size() == rle_bwt.size());
    all_good = true; miss_match = 0;
    for (std::size_t i = 0; i < v_bwt_1c.size(); i++)
    {
        all_good = all_good and (v_bwt_1c[i] == rle_bwt[i]);
        if (not all_good) { miss_match = i; break; }
    }
    REQUIRE(all_good);
    
    // Now check that the rlebwt from disk is correct with sdsl
    std::vector<char> yeast;
    read_fasta_file(std::string(testfiles_dir + "/yeast.fasta").c_str(), yeast, w_l1);
    
    uint8_t num_bytes = 1;
    sdsl::csa_wt<> csa;
    yeast.insert(yeast.begin(), w_l1 - 1, 3);
    yeast.push_back(0);
    sdsl::construct_im(csa, static_cast<const char *>(&yeast[0]), num_bytes);
    std::vector<char> sdsl_bwt;
    for (std::size_t i = 0; i < csa.size(); i++)
    {
        std::size_t adj_i = i + ((yeast.size()) - w_l1 + 1) % (yeast.size());
        sdsl_bwt.push_back(csa.bwt[adj_i]);
    }
    
    all_good = true; miss_match = 0;
    for (std::size_t i = 0; i < rle_bwt.size(); i++)
    {
        std::size_t adj_i = i + (yeast.size() - w_l1 + 1) % (yeast.size());
        all_good = all_good and (sdsl_bwt[adj_i] == rle_bwt[i]);
        //if (not all_good) { miss_match = i; break; }
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