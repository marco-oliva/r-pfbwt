//
// Created by marco on 10/21/22.
//

#include <rle/rle_string.hpp>
#include <pfp/pfp.hpp>

#include <CLI/CLI.hpp>
#include <version.hpp>
#include <spdlog/spdlog.h>

#include <rpfbwt_algorithm.hpp>

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

    rpfbwt::rpfbwt_algo<uint8_t> rpfbwt_algo(l1_prefix, w1, w2);
    rpfbwt_algo.l1_bwt(false);
}
