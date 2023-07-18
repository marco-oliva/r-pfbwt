//
// Created by marco on 10/21/22.
//

#include <CLI/CLI.hpp>
#include <version.hpp>

#include <rpfbwt_algorithm.hpp>

bool file_exists(std::string path)
{
    struct stat64 stat_buf;
    int rc = stat64(path.c_str(), &stat_buf);
    return rc == 0;
}

int main(int argc, char **argv)
{
    CLI::App app("rpfbwt");

    std::string l1_prefix;
    std::size_t w1;
    std::size_t w2;
    std::size_t threads = 1;
    std::string tmp_dir;
    std::size_t chunks = 50;
    uint32_t int_shift = 10;
    
    bool bwt_only = false;

    app.add_option("--l1-prefix", l1_prefix, "Level 1 Prefix.")->configurable()->required();
    app.add_option("--w1", w1, "Level 1 window length.")->configurable()->required();
    app.add_option("--w2", w2, "Level 2 window length.")->configurable()->required();
    app.add_option("-t, --threads", threads, "Number of threads.")->configurable();
    app.add_option("--chunks", chunks, "Number of chunks.")->configurable()->check(CLI::Range(1,1000));
    app.add_option("--integer-shift", int_shift, "Integer shift used during parsing.")->configurable();
    app.add_option("--tmp-dir", tmp_dir, "Temporary files directory.")->check(CLI::ExistingDirectory)->configurable();
    app.add_flag("--bwt-only", bwt_only, "Only compute the RLBWT. No SA values.")->configurable();
    app.add_flag_callback("--version",rpfbwt::Version::print,"Version number.");
    app.set_config("--configure");
    app.allow_windows_style_options();

    CLI11_PARSE(app, argc, argv);
    
    // Check all files needed
    std::string l1_d_path = l1_prefix + ".dict";
    if (not file_exists(l1_d_path)) { spdlog::error("Could not open required file: {}", l1_d_path); std::exit(EXIT_FAILURE); }
    
    std::string l2_d_path = l1_prefix + ".parse.dict";
    if (not file_exists(l2_d_path)) { spdlog::error("Could not open required file: {}", l2_d_path); std::exit(EXIT_FAILURE); }
    
    std::string l2_p_path = l1_prefix + ".parse.parse";
    if (not file_exists(l2_p_path)) { spdlog::error("Could not open required file: {}", l2_p_path); std::exit(EXIT_FAILURE); }
    
    // Set tmp dir for rle string
    if (not tmp_dir.empty()) { rle::TempFile::setDirectory(tmp_dir); }

    // Build dictionaries and parse
    std::less<uint8_t> uint8_t_comp;
    pfpds::dictionary<uint8_t> l1_d(l1_prefix, w1, uint8_t_comp, true, true, true, true, false, true, false);

    rpfbwt::rpfbwt_algo<uint8_t, uint32_t>::l2_colex_comp l2_comp(l1_d, int_shift);
    pfpds::dictionary<uint32_t, rpfbwt::rpfbwt_algo<uint8_t, uint32_t>::l2_colex_comp> l2_d(l1_prefix + ".parse", w2, l2_comp);
    pfpds::parse l2_p(l1_prefix + ".parse", l2_d.n_phrases() + 1);
    pfpds::pf_parsing<uint32_t, rpfbwt::rpfbwt_algo<uint8_t, uint32_t>::l2_colex_comp, pfpds::pfp_wt_sdsl> l2_pfp(l2_d, l2_p, false, true);
    
    rpfbwt::rpfbwt_algo<uint8_t> rpfbwt_algo(l1_prefix, l1_d, l2_pfp, int_shift, chunks);
    if (bwt_only) {rpfbwt_algo.l1_rlebwt(threads); }
    else { rpfbwt_algo.l1_refined_rindex(threads); }
}
