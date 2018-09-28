#ifndef FILE_HPP_
#define FILE_HPP_

#include <filesystem>
#include "simulation.hpp"

inline std::filesystem::path out_dir; // 保存先ディレクトリ

// 保存先ディレクトリのチェック
bool check_outdir(const std::filesystem::path &out_dir);
// バイト数をいいかんじの文字列にする
const std::string byte2str(const uintmax_t &size);
// .profファイルの読み書き
void load_data(const std::string &fname, simulation_t &sim);
void save_data(const std::string &fname, simulation_t &sim);

#endif
