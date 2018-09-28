#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <cmath>
#include <sksat/math/vector.hpp>

// OpenMP
#ifdef OPENMP
	#include <omp.h>
#endif

// 浮動小数点数
using Float = double;

// 3次元ベクトル量: vec_t -> sksat::math::vector

#endif
