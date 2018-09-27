#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <list>
#include <memory>
#include <filesystem>
#include <cmath>
#include <sksat/math/vector.hpp>

#ifdef OPENMP
#include <omp.h>
#endif

using Float = double;

// 3次元ベクトル量: vec_t -> sksat::math::vector

struct simulation_t {
	using vec3 = sksat::math::vector<Float>;
	enum particle_type {
		fluid,
		wall,
		ghost
	};
	std::vector<vec3> pos, vel, acc;
	std::vector<Float> press, pav;
	std::vector<particle_type> type;

	void resize(const size_t &s){ pos.resize(s); vel.resize(s); acc.resize(s); press.resize(s); pav.resize(s); type.resize(s); }
};

std::filesystem::path out_dir; // 保存先ディレクトリ

namespace params {
	// コンパイル時に決定しておく定数
	constexpr Float dt = 0.0005;				// 時間刻み
	constexpr Float time_max = 1.0;
	constexpr int dim = 3;						// 次元
	constexpr Float kinem_viscous = 0.000001;	// 動粘性係数
	const sksat::math::vector<Float> gravity = {0.0, 0.0, -9.8}; //TODO: sksat::math::vectorのconstexprコンストラクタ
	constexpr Float dens[] = {1000, 1000};		// 密度
	constexpr Float sound_vel = 22.0;			// 音速
	constexpr Float col_rat = 0.2;				// 接近した粒子の反発率
	constexpr Float col = 1.0 + col_rat;
	constexpr Float dst_lmt_rat = 0.9;			// これ以上の接近を許さない距離の係数
	constexpr Float crt_num = 0.1;				// クーラン条件数

	// 事前に計算しておく定数
	Float pcl_dst;			// 初期粒子間距離
	Float r, r2;			// 影響半径とその２乗
	Float rlim, rlim2;		// これ以上の粒子間の接近を許さない距離
	sksat::math::vector<Float> min, max; // 計算領域
	Float n0;				// 粒子数密度
	Float lamda;			// ラプラシアンモデルの係数λ
	Float coeff_viscous;	// 粘性項の係数
	Float coeff_mkpress;	// 圧力計算の係数
	Float coeff_press_grad; // 圧力勾配項の係数

	// 変数
	Float time;

	// バケット
	using bucket_t = std::list<int>;
	std::vector<std::shared_ptr<bucket_t>> bucket;
	Float bsize, bsize2, bsizeinv;
	int nbx, nby, nbz;
	int nbxy, nbxyz;
}

// 保存先ディレクトリのチェック
bool check_outdir(const std::filesystem::path &out_dir);
// バイト数をいいかんじの文字列にする
const std::string byte2str(const uintmax_t &size);
// .profファイルの読み書き
void load_data(const std::string &fname, simulation_t &sim);
void save_data(const std::string &fname, simulation_t &sim);

// 計算に関わる関数
// 事前に計算できるパラメータを計算
void set_param(const simulation_t &sim);
// 重み関数
inline Float weight(const Float dist, const Float re){
	return ((re/dist) - 1.0);
}

void check_particle(simulation_t &sim);

// バケット
inline std::shared_ptr<params::bucket_t> get_bucket(const int ix, const int iy, const int iz){
	int ib = (iz*params::nbxy) + (iy*params::nbx) + ix;
	return params::bucket[ib];
}
std::shared_ptr<params::bucket_t> get_bucket(const sksat::math::vector<Float> &pos);
void make_bucket(const simulation_t &sim);
const std::vector<std::shared_ptr<params::bucket_t>> near_buckets(const sksat::math::vector<Float> &pos);

// 計算のメインループ
void sim_loop(simulation_t &sim);
// 粘性項
void viscous_term(simulation_t &sim);
// 外力項
void external_term(simulation_t &sim);
// 仮の加速度で速度・位置の更新
void update_vp_tmp(simulation_t &sim);
// 衝突応答
void check_collision(simulation_t &sim);
// 仮圧力,圧力の計算
void make_press(simulation_t &sim);
// 圧力勾配項
void press_grad_term(simulation_t &sim);
// 速度・位置の修正
void update_vp(simulation_t &sim);

// arg: init.prof out_dir
int main(int argc, char **argv){
	namespace fs = std::filesystem;
	simulation_t sim;

	if(argc != 3) return -1;

#ifdef OPENMP
	std::cout << "[OpenMP] num_procs: " << omp_get_num_procs() << std::endl;
#endif

	// output directoryの準備
#ifndef BENCH
	out_dir = argv[2];
	if(!check_outdir(out_dir)) return -1;
#endif

	load_data(argv[1], sim);

	set_param(sim);

	sim_loop(sim);

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

// set_paramで使う関数
// 液体粒子の最小距離を求める
Float calc_min_dist(const simulation_t &sim){
	Float dist2 = 0.0;
	for(int i=0;i<sim.pos.size();i++){
		if(sim.type[i] != simulation_t::fluid) continue;
		for(int k=0;k<sim.pos.size();k++){
			if(i == k) continue;
			if(sim.type[k] != simulation_t::fluid) continue;
			auto &pos_i = sim.pos[i];
			auto &pos_k = sim.pos[k];
			auto vd = pos_k - pos_i;
			auto d2 = (vd.x*vd.x) + (vd.y*vd.y) + (vd.z*vd.z);
			if(dist2 == 0.0 || d2 < dist2) dist2 = d2;
		}
	}
	return std::sqrt(dist2);
}

void set_param(const simulation_t &sim){
	// 事前に計算できるパラメータを求める
	params::pcl_dst = calc_min_dist(sim); //TODO: とりあえず液体粒子の最小距離にしている．計算を途中で再開する場合これは使えない．

	params::r = params::pcl_dst * 2.1; // 影響半径
	params::r2= params::r * params::r;

	params::rlim = params::pcl_dst * params::dst_lmt_rat;
	params::rlim2= params::rlim * params::rlim;

	{	// とりあえずハードコード
		auto tmp = params::pcl_dst * 3;
		params::min = {
			0.0-tmp,
			0.0-tmp,
			0.0-tmp
		};
		params::max = {
			1.0+tmp,
			0.2+tmp,
			0.6+params::pcl_dst*30
		};
	}

	// bucket
	params::bsize = params::r * (1.0 + params::crt_num);
	params::bsize2= params::bsize * params::bsize;
	params::bsizeinv= 1.0 / params::bsize;
	{
		auto tmp = params::max - params::min;
		tmp.x /= params::bsize;
		tmp.y /= params::bsize;
		tmp.z /= params::bsize;
		params::nbx = static_cast<int>(tmp.x) + 3;
		params::nby = static_cast<int>(tmp.y) + 3;
		params::nbz = static_cast<int>(tmp.z) + 3;
	}
	params::nbxy = params::nbx * params::nby;
	params::nbxyz= params::nbx * params::nby * params::nbz;
	params::bucket.reserve(params::nbxyz);
	for(int i=0;i<params::nbxyz;i++)
		params::bucket.push_back(std::make_shared<params::bucket_t>());

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

	// 粘性項の係数
	params::coeff_viscous = 2.0 * params::kinem_viscous * params::dim / (params::n0 * params::lamda);

	// 圧力計算の係数
	params::coeff_mkpress = params::sound_vel*params::sound_vel / params::n0;

	// 圧力勾配項の係数
	params::coeff_press_grad = -1.0 * params::dim / params::n0;

	params::time = 0.0;

	// パラメータの表示
	std::cout
		<< "parameters:" << std::endl
		<< "\tparticle distance: " << params::pcl_dst << std::endl
		<< "\tr: " << params::r << std::endl
		<< "\trlim: " << params::rlim << std::endl
		<< "\tmin: (" << params::min.x << ", " << params::min.y << ", " << params::min.z << ")" << std::endl
		<< "\tmax: (" << params::max.x << ", " << params::max.y << ", " << params::max.z << ")" << std::endl
		<< "\tbucket size: " << params::bsize << std::endl
		<< "\tbucket num: " << params::nbxyz
		<< " (" << params::nbx << ", " << params::nby << ", " << params::nbz << ")" << std::endl
		<< "\tn0: " << params::n0 << std::endl
		<< "\tλ : " << params::lamda << std::endl
		<< "\tviscous coeff: " << params::coeff_viscous << std::endl
		<< "\tmkpress coeff: " << params::coeff_mkpress << std::endl
		<< "\ttime: " << params::time << std::endl
		<< std::endl;
}

void check_particle(simulation_t &sim){
	for(int i=0;i<sim.pos.size();i++){
		auto &pos = sim.pos[i];
		auto &vel = sim.vel[i];
		auto &type= sim.type[i];
		if(	pos.x < params::min.x || params::max.x < pos.x ||
			pos.y < params::min.y || params::max.y < pos.y ||
			pos.z < params::min.z || params::max.z < pos.z){
			type = simulation_t::ghost;
			vel = {0.0, 0.0, 0.0};
		}
	}
}

std::shared_ptr<params::bucket_t> get_bucket(const sksat::math::vector<Float> &pos){
	auto p = pos - params::min;
	p = p * params::bsizeinv;
	int ix = static_cast<int>(p.x) + 1;
	int iy = static_cast<int>(p.y) + 1;
	int iz = static_cast<int>(p.z) + 1;
	return get_bucket(ix, iy, iz);
}

void make_bucket(const simulation_t &sim){
	for(auto &b : params::bucket) b->clear();

	for(int i=0;i<sim.pos.size();i++){
		if(sim.type[i] == simulation_t::ghost) continue;
		auto b = get_bucket(sim.pos[i]);
		b->push_back(i);
	}
}

const std::vector<std::shared_ptr<params::bucket_t>> near_buckets(const sksat::math::vector<Float> &pos){
	std::vector<std::shared_ptr<params::bucket_t>> neigh(3*3*3);
	neigh.clear();
	auto p = (pos - params::min) * params::bsizeinv;
	int ix = static_cast<int>(p.x) + 1;
	int iy = static_cast<int>(p.y) + 1;
	int iz = static_cast<int>(p.z) + 1;
	for(int dz=-1;dz<=1;dz++){
		for(int dy=-1;dy<=1;dy++){
			for(int dx=-1;dx<=1;dx++){
				neigh.push_back(get_bucket(ix+dx, iy+dy, iz+dz));
			}
		}
	}
	return neigh;
}

void sim_loop(simulation_t &sim){
	size_t iloop = 0; // ループの回数
	size_t ifile = 0; // 保存ファイルの番号 TODO: 途中から計算する場合の初期値の設定
	using clock = std::chrono::high_resolution_clock;
	clock::time_point begin, end;

	std::cout << "start simulation." << std::endl;

	begin = clock::now();
	// メインループ
	while(true){
#ifndef BENCH
		// ログ表示
		if(iloop % 10 == 0){
			std::cout << "iloop=" << iloop << ", time=" << params::time << std::endl;
		}
#endif // BENCH

		// ファイル保存
		if(iloop % 10 == 0){
#ifndef BENCH
			std::stringstream fname;
			fname << "output"
				<< std::setfill('0') << std::setw(10)
				<< iloop/100
				<< ".prof";
			auto fpath = out_dir / fname.str();
			save_data(fpath, sim);
#endif // BENCH
			// 終了条件
			if(params::time >= params::time_max) break;
		}

		// 計算部分

		// バケット生成
		make_bucket(sim);

		// 仮の加速度 <- 粘性項
		viscous_term(sim);

		// 仮の加速度 <- 外力(重力)項
		external_term(sim);

		// 仮の速度,仮の位置を更新
		update_vp_tmp(sim);
		check_particle(sim);

		// 衝突判定
		check_collision(sim);

		// 仮の圧力
		make_press(sim);

		// 加速度の修正量 <- 仮の圧力
		press_grad_term(sim);

		// 速度,位置を修正
		update_vp(sim);
		check_particle(sim);

		// 圧力の修正
		make_press(sim);

		for(int i=0;i<sim.pos.size();i++)
			sim.pav[i] += sim.press[i];

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

// 粘性項
void viscous_term(simulation_t &sim){
	const sksat::math::vector<Float> vec_zero = {0.0, 0.0, 0.0};
#ifdef OPENMP
	#pragma omp parallel for
#endif
	for(int i=0;i<sim.pos.size();i++){
		if(sim.type[i] != simulation_t::fluid) continue;
		auto acc = vec_zero;
		// 周囲の粒子との計算 TODO: バケット構造で近傍粒子探索
		auto near = near_buckets(sim.pos[i]);
		for(auto &b : near){ // 周囲のバケット
			for(auto begin = b->begin(); begin!=b->end(); begin++){
				auto &k = *begin;
//			for(auto &k : *b){
				if(i == k) continue;
				if(sim.type[k] == simulation_t::ghost) continue; // ゴースト粒子は無視
				auto vd = sim.pos[k] - sim.pos[i];
				auto dist2 = (vd.x*vd.x) + (vd.y*vd.y) + (vd.z*vd.z);
				if(params::r2 <= dist2) continue; // 影響半径外の粒子は無視
				auto dist = std::sqrt(dist2);
				auto w = weight(dist, params::r); // 重み
				auto vel_diff= sim.vel[k] - sim.vel[i];
				acc += vel_diff * w;
			}
		}

		sim.acc[i] = acc * params::coeff_viscous;
	}
}

// 外力項
void external_term(simulation_t &sim){
#ifdef OPENMP
	#pragma omp parallel for
#endif
	for(int i=0;i<sim.pos.size();i++){
		if(sim.type[i] != simulation_t::fluid) continue;
		sim.acc[i] += params::gravity; // 重力ベクトル
	}
}

// 仮の加速度を使って速度と位置を更新する
void update_vp_tmp(simulation_t &sim){
#ifdef OPENMP
	#pragma omp parallel for
#endif
	for(int i=0;i<sim.pos.size();i++){
		if(sim.type[i] != simulation_t::fluid) continue;
		sim.vel[i] += sim.acc[i] * params::dt;
		sim.pos[i] += sim.vel[i] * params::dt;
		sim.acc[i] = {0.0, 0.0, 0.0};
	}
}

// 剛体衝突
void check_collision(simulation_t &sim){
#ifdef OPENMP
#pragma omp parallel
#endif
{
#ifdef OPENMP
	#pragma omp for
#endif
	for(int i=0;i<sim.pos.size();i++){
		if(sim.type[i] != simulation_t::fluid) continue;
		sksat::math::vector<Float> v = sim.vel[i];
		auto near = near_buckets(sim.pos[i]);
		for(auto &b : near){ // 周囲のバケット
			for(auto &k : *b){
				if(i == k) continue;
				if(sim.type[k] == simulation_t::ghost) continue;
				auto pd = sim.pos[k] - sim.pos[i];
				auto dist2 = (pd.x*pd.x) + (pd.y*pd.y) + (pd.z*pd.z);
				if(params::rlim2 <= dist2) continue;
				auto vd = sim.vel[i] - sim.vel[k];
				auto fDT = vd.x*pd.x + vd.y*pd.y + vd.z*pd.z;
				if(fDT > 0.0){
					fDT *= params::col * params::dens[sim.type[k]]
						/ ((params::dens[sim.type[i]]+params::dens[sim.type[k]])*dist2);
					v -= pd*fDT;
				}
			}
		}
		sim.acc[i] = v;
	}
#ifdef OPENMP
	#pragma omp for
#endif
	for(int i=0;i<sim.pos.size();i++){
		sim.vel[i] = sim.acc[i];
	}
} // omp parallel
}

// 粒子数密度から仮の圧力を求める
void make_press(simulation_t &sim){
#ifdef OPENMP
	#pragma omp parallel for
#endif
	for(int i=0;i<sim.pos.size();i++){
		Float ni = 0.0; // 粒子数密度
		if(sim.type[i] == simulation_t::ghost) continue;
		auto near = near_buckets(sim.pos[i]);
		for(auto b_ = near.begin(); b_!=near.end(); b_++){
			auto& b = *b_;
//		for(auto &b : near){ // 周囲のバケット
			for(auto &k : *b){
				if(i == k) continue;
				if(sim.type[k] == simulation_t::ghost) continue;
				auto pd = sim.pos[k] - sim.pos[i];
				auto dist2 = (pd.x*pd.x) + (pd.y*pd.y) + (pd.z*pd.z);
				if(params::r2 <= dist2) continue;
				auto dist = std::sqrt(dist2);
				ni += weight(dist, params::r);
			}
		}
		sim.press[i] = (ni > params::n0)		// ni>n0なら内部粒子，そうでなければ自由表面(圧力0)
			* (ni - params::n0)
			* params::coeff_mkpress
			* params::dens[sim.type[i]];
	}
}

// 圧力勾配項
void press_grad_term(simulation_t &sim){
#ifdef OPENMP
	#pragma omp parallel for
#endif
	for(int i=0;i<sim.pos.size();i++){
		if(sim.type[i] != simulation_t::fluid) continue;
		Float press_min = 0.0;
		// 影響半径内の粒子の圧力の最小値を求める
		for(auto &b : near_buckets(sim.pos[i])){ // 周囲のバケット
			for(auto &k : *b){
				if(i == k) continue;
				if(sim.type[k] == simulation_t::ghost) continue;
				auto pd = sim.pos[k] - sim.pos[i];
				auto dist2 = (pd.x*pd.x) + (pd.y*pd.y) + (pd.z*pd.z);
				if(params::r2 <= dist2) continue;
				if(press_min > sim.press[i]) press_min = sim.press[i];
			}
		}
		sksat::math::vector<Float> acc = {0.0, 0.0, 0.0};

		for(auto &b : near_buckets(sim.pos[i])){ // 周囲のバケット
			for(auto &k : *b){
				if(i == k) continue;
				if(sim.type[k] == simulation_t::ghost) continue;
				auto pd = sim.pos[k] - sim.pos[i];
				auto dist2 = (pd.x*pd.x) + (pd.y*pd.y) + (pd.z*pd.z);
				if(params::r2 <= dist2) continue;
				auto dist = std::sqrt(dist2);
				auto w = weight(dist, params::r);
				acc += pd * (w * (sim.press[k] - press_min) / dist2);
			}
		}
		sim.acc[i] = acc * (params::coeff_press_grad / params::dens[simulation_t::fluid]);
	}
}

void update_vp(simulation_t &sim){
#ifdef OPENMP
	#pragma omp parallel for
#endif
	for(int i=0;i<sim.pos.size();i++){
		if(sim.type[i] != simulation_t::fluid) continue;
		sim.vel[i] += sim.acc[i] * params::dt;
		sim.pos[i] += sim.acc[i] * params::dt * params::dt;
		sim.acc[i] = {0.0, 0.0, 0.0};
	}
}
