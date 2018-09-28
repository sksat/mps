#ifndef PARAMS_HPP_
#define PARAMS_HPP_

#include <iostream>
#include <list>
#include <vector>
#include <memory>
#include "common.hpp"

namespace params {
	// コンパイル時に決定しておく定数
	constexpr Float dt = 0.0005;				// 時間刻み
	constexpr Float time_max = 0.1;
	constexpr int dim = 3;						// 次元
	constexpr Float kinem_viscous = 0.000001;	// 動粘性係数
	inline const sksat::math::vector<Float> gravity = {0.0, 0.0, -9.8}; //TODO: sksat::math::vectorのconstexprコンストラクタ
	constexpr Float dens[] = {1000, 1000};		// 密度
	constexpr Float sound_vel = 22.0;			// 音速
	constexpr Float col_rat = 0.2;				// 接近した粒子の反発率
	constexpr Float col = 1.0 + col_rat;
	constexpr Float dst_lmt_rat = 0.9;			// これ以上の接近を許さない距離の係数
	constexpr Float crt_num = 0.1;				// クーラン条件数

	// 事前に計算しておく定数
	inline Float pcl_dst;			// 初期粒子間距離
	inline Float r, r2;			// 影響半径とその２乗
	inline Float rlim, rlim2;		// これ以上の粒子間の接近を許さない距離
	inline sksat::math::vector<Float> min, max; // 計算領域
	inline Float n0;				// 粒子数密度
	inline Float lamda;			// ラプラシアンモデルの係数λ
	inline Float coeff_viscous;	// 粘性項の係数
	inline Float coeff_mkpress;	// 圧力計算の係数
	inline Float coeff_press_grad; // 圧力勾配項の係数

	// 変数
	inline Float time;

	// バケット
	using bucket_t = std::list<int>;
	inline std::vector<std::shared_ptr<bucket_t>> bucket;
	inline Float bsize, bsize2, bsizeinv;
	inline int nbx, nby, nbz;
	inline int nbxy, nbxyz;
}

// 計算に関わる関数
// 事前に計算できるパラメータを計算
struct simulation_t;
void set_param(const simulation_t &sim);

#endif
