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
	enum type_t {
		liquid,
		wall,
	};
	type_t type;
	vec_t pos, vel, acc;
};

// .profファイルの読み込み
void load_data(const std::string &fname, std::vector<particle_t> &particle);
void save_data(const std::string &fname, std::vector<particle_t> &particle);

int main(int argc, char **argv){
	std::vector<particle_t> particle;

	if(argc != 2) return -1;

	load_data(argv[1], particle);
	save_data("test.prof", particle);
	return 0;
}

void load_data(const std::string &fname, std::vector<particle_t> &particle){
	std::ifstream ifs(fname);
	int num;
	int liquid = 0, wall =0;
	ifs >> num;
	particle.resize(num);

	std::cout << "loading \"" << fname;

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

void save_data(const std::string &fname, std::vector<particle_t> &particle){
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
