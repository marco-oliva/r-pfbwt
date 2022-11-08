//
// Created by marco on 10/21/22.
//

#include <CLI/CLI.hpp>
#include <version.hpp>

#include <rpfbwt_algorithm.hpp>

int main(int argc, char **argv)
{
    CLI::App app("rpfbwt");

    std::string l1_prefix;
    std::size_t w1;
    std::size_t w2;
    std::size_t threads = 1;
    std::string tmp_dir;
    std::size_t chunks_per_thread = 5;

    app.add_option("--l1-prefix", l1_prefix, "Level 1 Prefix.")->configurable()->required();
    app.add_option("--w1", w1, "Level 1 window length.")->configurable()->required();
    app.add_option("--w2", w2, "Level 2 window length.")->configurable()->required();
    app.add_option("-t, --threads", threads, "Number of threads.")->configurable();
    app.add_option("--chunks-per-thread", chunks_per_thread, "Number of chunks per thread.")->configurable();
    app.add_option("--tmp-dir", tmp_dir, "Temporary files directory.")->check(CLI::ExistingDirectory)->configurable();
    app.add_flag_callback("--version",rpfbwt::Version::print,"Version number.");
    app.set_config("--configure");
    app.allow_windows_style_options();

    CLI11_PARSE(app, argc, argv);
    
    // Set tmp dir for rle string
    if (not tmp_dir.empty()) { rle::TempFile::setDirectory(tmp_dir); }
    
    // Set threads accordingly to configuration
    omp_set_num_threads(threads);
    
    rpfbwt::rpfbwt_algo<uint8_t> rpfbwt_algo(l1_prefix, w1, w2, threads * chunks_per_thread);
    rpfbwt_algo.l1_bwt_parallel();
}
