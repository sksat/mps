#include <fstream>
#include "file.hpp"

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

void load_data(const std::string &fname, simulation_t &sim){
	std::ifstream ifs(fname);
	int num;
	int fluid = 0, wall =0;
	ifs >> num;
	sim.resize(num);

	std::cout << "loading \"" << fname << "\"";

	for(int i=0;i<num;i++){
		int n, t;
		auto &pos = sim.pos[i];
		auto &vel = sim.vel[i];
		auto &acc = sim.acc[i];
		auto &type= sim.type[i];
		auto &prs = sim.press[i], &pav = sim.pav[i];
		ifs >> n >> t
			>> pos.x >> pos.y >> pos.z
			>> vel.x >> vel.y >> vel.z
			>> prs >> pav;
		type = static_cast<simulation_t::particle_type>(t);
		if(type == simulation_t::fluid) fluid++;
		else if(type == simulation_t::wall) wall++;
	}

	std::cout << std::endl;
	std::cout << "particle num: " << sim.pos.size()
		<< " (fluid: " << fluid << ", wall: " << wall << ")" << std::endl;
}

void save_data(const std::string &fname, simulation_t &sim){
	std::ofstream ofs(fname);
	ofs << sim.pos.size() << std::endl;
	for(int i=0;i<sim.pos.size();i++){
		auto &type= sim.type[i];
		auto &pos = sim.pos[i];
		auto &vel = sim.vel[i];
		auto &prs = sim.press[i];
		auto &pav = sim.pav[i];
		pav /= 100;
		ofs << i << " " << type << " "
			<< pos.x << " " << pos.y << " " << pos.z << " "
			<< vel.x << " " << vel.y << " " << vel.z << " "
			<< prs << " " << pav << std::endl;
		pav = 0.0;
	}
}
