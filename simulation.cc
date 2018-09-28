#include <chrono>
#include <iomanip>
#include <sstream>
#include "simulation.hpp"
#include "file.hpp"

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
