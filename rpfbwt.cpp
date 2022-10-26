//
// Created by marco on 10/21/22.
//

#include <rle/rle_string.hpp>
#include <pfp/pfp.hpp>

#include <CLI/CLI.hpp>
#include <version.hpp>
#include <spdlog/spdlog.h>

int main(int argc, char **argv)
{
    CLI::App app("rpfbwt");

    std::string l1_prefix;
    std::size_t w1;
    std::size_t w2;

    app.add_option("--l1-prefix", l1_prefix, "Level 1 Prefix.")->configurable()->required();
    app.add_option("--w1", w1, "Level 1 window length.")->configurable()->required();
    app.add_option("--w2", w2, "Level 2 window length.")->configurable()->required();
    app.add_flag_callback("--version",rpfbwt::Version::print,"Version number.");
    app.set_config("--configure");
    app.allow_windows_style_options();

    CLI11_PARSE(app, argc, argv);

    // get L1 dictionary, occ and M1
    pfpds::dictionary<uint8_t> l1_d(l1_prefix, w1);

    size_t d1_words; uint_t * occ;
    std::vector<uint_t> l1_freq( 1, 0);
    pfpds::read_file<uint_t> (std::string(l1_prefix + ".occ").c_str(), occ, d1_words); // here occ has the same size as the integers used for gsacak.
    l1_freq.insert(l1_freq.end(),occ, occ + d1_words);
    
    // get L2 PFP
    std::string l2_prefix = l1_prefix + ".parse";
    pfpds::pf_parsing<uint32_t> l2_pfp(l2_prefix, w2);

    // Iterate on the suffixes of l1_d and output easy suffixes to a run length encoded string
    rle::RLEString::RLEncoder easy_suffixes_encoder(l1_prefix + ".easy.rle");

    
}
