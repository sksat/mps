#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <filesystem>
#include <cmath>
#include <sksat/math/vector.hpp>

using Float = double;

// 3次元ベクトル量: vec_t -> sksat::math::vector

// 粒子
struct particle_t {
	enum type_t {
		liquid,
		wall,
		ghost,
	};
	type_t type;
	sksat::math::vector<Float> pos, vel, acc;
};

std::filesystem::path out_dir; // 保存先ディレクトリ

namespace params {
	// コンパイル時に決定しておく定数
	constexpr Float dt = 0.0001;
	constexpr Float time_max = 1.0;
	constexpr int dim = 3;
	const sksat::math::vector gravity = {0.0, -9.8, 0.0}; //TODO: sksat::math::vectorのconstexprコンストラクタ
	constexpr Float sound_vel = 22.0; // 音速

	// 事前に計算しておく定数
	Float pcl_dst;	// 初期粒子間距離
	Float r, r2;	// 影響半径とその２乗
	Float n0;		// 粒子数密度
	Float lamda;	// ラプラシアンモデルの係数λ

	// 変数
	Float time;
}

// 保存先ディレクトリのチェック
bool check_outdir(const std::filesystem::path &out_dir);
// バイト数をいいかんじの文字列にする
const std::string byte2str(const uintmax_t &size);
// .profファイルの読み書き
void load_data(const std::string &fname, std::vector<particle_t> &particle);
void save_data(const std::string &fname, const std::vector<particle_t> &particle);

// 計算に関わる関数
// 事前に計算できるパラメータを計算
void set_param(const std::vector<particle_t> &particle);
// 重み関数
inline Float weight(const Float dist, const Float re){
	return ((re/dist) - 1.0);
}
// 計算のメインループ
void sim_loop(std::vector<particle_t> &particle);

// arg: init.prof out_dir
int main(int argc, char **argv){
	namespace fs = std::filesystem;
	std::vector<particle_t> particle;

	if(argc != 3) return -1;

	// output directoryの準備
	out_dir = argv[2];
	if(!check_outdir(out_dir)) return -1;

	load_data(argv[1], particle);

	set_param(particle);

	sim_loop(particle);

	return 0;
}

bool check_outdir(const std::filesystem::path &out_dir){
	using namespace std::filesystem;
	if(exists(out_dir)){
		std::cerr << out_dir << " is already exists." << std::endl;
		if(!is_directory(out_dir)){
			std::cerr << out_dir
				<< " is not a directory. please remove or rename." << std::endl;
			return false;
		}
	}else{
		std::cerr << "directory " << out_dir << " does not exists." << std::endl;
		std::cout << "creating " << out_dir << "... ";
		if(!create_directory(out_dir)){
			std::cerr << "failed." << std::endl;
			return false;
		}
		std::cout << "ok." << std::endl;
	}

	auto si = space(out_dir);
	std::cout << "output directory: " << out_dir
		<< " (capacity: " << byte2str(si.capacity)
		<< ", free: " << byte2str(si.free)
		<< ", available: " << byte2str(si.available)
		<< ")" << std::endl;

	//TODO: ディレクトリ容量が足りなくなりそうだったら知らせる

	return true;
}

const std::string byte2str(const uintmax_t &size){
	const auto KB = 1024;
	const auto MB = 1024 * KB;
	const auto GB = 1024 * MB;
	if(size > GB)
		return std::to_string(size / GB) + "GB";
	else if(size > MB)
		return std::to_string(size / MB) + "MB";
	else if(size > KB)
		return std::to_string(size / KB) + "KB";
	else
		return std::to_string(size) + "B";
}

void load_data(const std::string &fname, std::vector<particle_t> &particle){
	std::ifstream ifs(fname);
	int num;
	int liquid = 0, wall =0;
	ifs >> num;
	particle.resize(num);

	std::cout << "loading \"" << fname << "\"";

	for(int i=0;i<num;i++){
		int n, type;
		auto &p = particle[i];
		int t;
		Float prs, pav;
		ifs >> n >> t
			>> p.pos.x >> p.pos.y >> p.pos.z
			>> p.vel.x >> p.vel.y >> p.vel.z
			>> prs >> pav;
		p.type = static_cast<particle_t::type_t>(t);
		if(p.type == particle_t::liquid) liquid++;
		else if(p.type == particle_t::wall) wall++;
	}

	std::cout << std::endl;
	std::cout << "particle num: " << particle.size()
		<< " (liquid: " << liquid << ", wall: " << wall << ")" << std::endl;
}

void save_data(const std::string &fname, const std::vector<particle_t> &particle){
	std::ofstream ofs(fname);
	ofs << particle.size() << std::endl;
	for(int i=0;i<particle.size();i++){
		auto &p = particle[i];
		ofs << i << " " << p.type << " "
			<< p.pos.x << " " << p.pos.y << " " << p.pos.z << " "
			<< p.vel.x << " " << p.vel.y << " " << p.vel.z << " "
			<< 0.0 << " " << 0.0 << std::endl;
	}
}

// set_paramで使う関数
// 液体粒子の最小距離を求める
Float calc_min_dist(const std::vector<particle_t> &particle){
	Float dist2 = 0.0;
	for(int i=0;i<particle.size();i++){
		if(particle[i].type != particle_t::liquid) continue;
		for(int k=0;k<particle.size();k++){
			if(i == k) continue;
			if(particle[k].type != particle_t::liquid) continue;
			auto &pos_i = particle[i].pos;
			auto &pos_k = particle[k].pos;
			auto vd = pos_k - pos_i;
			auto d2 = (vd.x*vd.x) + (vd.y*vd.y) + (vd.z*vd.z);
			if(dist2 == 0.0 || d2 < dist2) dist2 = d2;
		}
	}
	return std::sqrt(dist2);
}

void set_param(const std::vector<particle_t> &particle){
	// 事前に計算できるパラメータを求める
	params::pcl_dst = calc_min_dist(particle); //TODO: とりあえず液体粒子の最小距離にしている．計算を途中で再開する場合これは使えない．

	params::r = params::pcl_dst * 2.1; // 影響半径
	params::r2= params::r * params::r;

	// 粒子を仮想的に格子状に配置
	// n0:		粒子数密度の基準値
	// lamda:	ラプラシアンモデルの係数λ
	Float tn0 = 0.0, tlmd = 0.0;
	for(int ix=-4;ix<5;ix++){
		for(int iy=-4;iy<5;iy++){
			for(int iz=-4;iz<5;iz++){
				sksat::math::vector<Float> d;
				d.x = static_cast<Float>(ix);
				d.y = static_cast<Float>(iy);
				d.z = static_cast<Float>(iz);
				d = params::pcl_dst * d;
				auto d2 = (d.x*d.x) + (d.y*d.y) + (d.z*d.z);
				if(d2 <= params::r2){
					if(d2 == 0.0) continue;
					Float dist = std::sqrt(d2);
					auto w = weight(dist, params::r);
					tn0 += w;
					tlmd+= d2 * w;
				}
			}
		}
	}
	params::n0 = tn0;
	params::lamda = tlmd / params::n0;

	params::time = 0.0;

	// パラメータの表示
	std::cout
		<< "parameters:" << std::endl
		<< "\tparticle distance: " << params::pcl_dst << std::endl
		<< "\tn0: " << params::n0 << std::endl
		<< "\tλ : " << params::lamda << std::endl
		<< "\ttime: " << params::time << std::endl
		<< std::endl;
}

void sim_loop(std::vector<particle_t> &particle){
	size_t iloop = 0; // ループの回数
	size_t ifile = 0; // 保存ファイルの番号 TODO: 途中から計算する場合の初期値の設定
	using clock = std::chrono::high_resolution_clock;
	clock::time_point begin, end;

	std::cout << "start simulation." << std::endl;

	begin = clock::now();
	// メインループ
	while(true){
		// ログ表示
		if(iloop % 100 == 0){
			std::cout << "iloop=" << iloop << ", time=" << params::time << std::endl;
		}

		// ファイル保存
		if(iloop % 100 == 0){
			std::stringstream fname;
			fname << "output"
				<< std::setfill('0') << std::setw(10)
				<< iloop/100
				<< ".prof";
			auto fpath = out_dir / fname.str();
			save_data(fpath, particle);
			// 終了条件
			if(params::time >= params::time_max) break;
		}

		iloop++;
		params::time += params::dt;
	}
	end = clock::now();
	auto diff = end - begin;
	auto mdiff= std::chrono::duration_cast<std::chrono::milliseconds>(diff);
	auto msec = static_cast<double>(mdiff.count());

	std::cout << "end simulation." << std::endl
		<< "simulation time: "
		<< (msec / 1000.0) << "sec" << std::endl;
}
