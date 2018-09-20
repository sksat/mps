#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <cmath>
#include <sksat/math/vector.hpp>

using Float = double;

// 定数
constexpr Float dt = 0.0001;
constexpr Float time_max = 1.0;

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

// arg: init.prof out_dir
int main(int argc, char **argv){
	namespace fs = std::filesystem;
	std::vector<particle_t> particle;

	if(argc != 3) return -1;

	// output directoryの準備
	fs::path out_dir(argv[2]);
	if(!check_outdir(out_dir)) return -1;

	load_data(argv[1], particle);

	set_param(particle);

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
	auto pcl_dst = calc_min_dist(particle); //TODO: とりあえず液体粒子の最小距離にしている．計算を途中で再開する場合これは使えない．

	// パラメータの表示
	std::cout << "particle distance: " << pcl_dst << std::endl;
}
