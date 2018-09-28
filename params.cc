#include "params.hpp"
#include "simulation.hpp"

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
