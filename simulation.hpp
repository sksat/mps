#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#include <sksat/math/vector.hpp>
#include "common.hpp"
#include "params.hpp"

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

// 計算に関わる関数
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

#endif
