/*
*  implicit.h
*  flip3D
*/
//------------------------------------------------------------------------------
//reference to https://code.google.com/archive/p/flip3d/ ando`s filp code.
//I just use this for surface reconstruction.
//------------------------------------------------------------------------------
#pragma once
#ifndef IMPLICIT_H
#define IMPLICIT_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include <vector>
#include "Sorter.h"
namespace implicit {
	double implicit_func(sorter *sort, double p[3], double density);
};


using namespace std;

static double implicit_func(vector<Particle *> &neighbors, double p[3], double density, int gn) {
	double phi = 8.0*density / gn;
	for (int m = 0; m<neighbors.size(); m++) 
	{
		Particle &np = *neighbors[m];
		//if (np.type == WALL) 
		//{
		//	if (length(np.p, p) < density / gn) return 4.5*density / gn;
		//	continue;
		//}
		double p2[3];
		p2[0] = np._p.X;
		p2[1] = np._p.Y;
		p2[2] = np._p.Z;
		double d = length(p2, p);
		if (d < phi) {
			phi = d;
		}
	}
	return phi - density / gn;
}

static double implicit::implicit_func(sorter *sort, double p[3], double density) {
	int gn = sort->getCellSizeX();
	vector<Particle *> neighbors = sort->getNeigboringParticles_cell(
		fmax(0, fmin(gn - 1, gn*p[0])),
		fmax(0, fmin(gn - 1, gn*p[1])),
		fmax(0, fmin(gn - 1, gn*p[2])), 2, 2, 2
		);
	return implicit_func(neighbors, p, density, gn);
}
#endif