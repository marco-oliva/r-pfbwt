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
    
    bool bwt_only = false;

    app.add_option("--l1-prefix", l1_prefix, "Level 1 Prefix.")->configurable()->required();
    app.add_option("--w1", w1, "Level 1 window length.")->configurable()->required();
    app.add_option("--w2", w2, "Level 2 window length.")->configurable()->required();
    app.add_option("-t, --threads", threads, "Number of threads.")->configurable();
    app.add_option("--chunks", chunks, "Number of chunks.")->configurable()->check(CLI::Range(1,1000));
    app.add_option("--tmp-dir", tmp_dir, "Temporary files directory.")->check(CLI::ExistingDirectory)->configurable();
    app.add_flag("--bwt-only", bwt_only, "Only compute the RLBWT. No SA values.")->configurable();
    app.add_flag_callback("--version",rpfbwt::Version::print,"Version number.");
    app.set_config("--configure");
    app.allow_windows_style_options();

    CLI11_PARSE(app, argc, argv);
    
    // Check all files needed
    std::string l1_d_path = l1_prefix + ".dict";
    if (not file_exists(l1_d_path)) { spdlog::error("Could not open required file: {}", l1_d_path); std::exit(EXIT_FAILURE); }
    
    std::string l1_occ_path = l1_prefix + ".occ";
    if (not file_exists(l1_occ_path)) { spdlog::error("Could not open required file: {}", l1_occ_path); std::exit(EXIT_FAILURE); }
    
    std::string l2_d_path = l1_prefix + ".parse.dict";
    if (not file_exists(l2_d_path)) { spdlog::error("Could not open required file: {}", l2_d_path); std::exit(EXIT_FAILURE); }
    
    std::string l2_p_path = l1_prefix + ".parse.parse";
    if (not file_exists(l2_p_path)) { spdlog::error("Could not open required file: {}", l2_p_path); std::exit(EXIT_FAILURE); }
    
    // Set tmp dir for rle string
    if (not tmp_dir.empty()) { rle::TempFile::setDirectory(tmp_dir); }
    
    rpfbwt::rpfbwt_algo<uint8_t> rpfbwt_algo(l1_prefix, w1, w2, chunks);
    if (bwt_only) {rpfbwt_algo.l1_rlebwt(threads); }
    else { rpfbwt_algo.l1_refined_rindex(threads); }
    
}
