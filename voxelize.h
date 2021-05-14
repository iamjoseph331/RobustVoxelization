#ifndef VOXELIZE_H
#define VOXELIZE_H

#define INF 99999999.0

struct tetrahedra
{
	Eigen::Vector3f a, b, c, d;
	Eigen::Vector3f bound_min, bound_max;
};

#endif