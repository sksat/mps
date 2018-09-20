#include <iostream>
#include <fstream>
#include <vector>

using Float = double;

// 定数
constexpr Float dt = 0.0001;
constexpr Float time_max = 1.0;

// 3次元ベクトル量
struct vec_t {
	Float x, y, z;
};

// 粒子
struct particle_t {
	vec_t pos, vel, acc;
};

// .profファイルの読み込み
void load_data(const std::string &fname, std::vector<particle_t> &particle);

int main(int argc, char **argv){
	std::vector<particle_t> particle;

	if(argc != 2) return -1;

	load_data(argv[1], particle);
	std::cout << "particle num: " << particle.size() << std::endl;
	return 0;
}

void load_data(const std::string &fname, std::vector<particle_t> &particle){
	std::ifstream ifs(fname);
	int num;
	ifs >> num;
	particle.resize(num);

	for(int i=0;i<num;i++){
		int n, type;
		auto &p = particle[i];
		Float prs, pav;
		ifs >> n >> type
			>> p.pos.x >> p.pos.y >> p.pos.z
			>> p.vel.x >> p.vel.y >> p.vel.z
			>> prs >> pav;
	}
}
