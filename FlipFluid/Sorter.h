/*
*  sorter.cpp
*  flip3D
*/
//------------------------------------------------------------------------------
//reference to https://code.google.com/archive/p/flip3d/ ando`s filp code.
//I just use this for surface reconstruction.
//------------------------------------------------------------------------------
#pragma once

#ifndef _SORTER_H
#define _SORTER_H

#include <vector>
#include "util.h"
using namespace std;
extern ofstream debug;
class sorter {
public:
	sorter(int _gnx, int _gny, int _gnz);
	~sorter();

	void sort(std::vector<Particle> &Particles);
	std::vector<Particle *> getNeigboringParticles_wall(int i, int j, int k, int w, int h, int d);
	std::vector<Particle *> getNeigboringParticles_cell(int i, int j, int k, int w, int h, int d);
	double levelset(int i, int j, int k, double ***halfwall, double density);

	int	 getCellSizeX(){ return gnx; }
	int	 getCellSizeY(){ return gnx; }
	int	 getCellSizeZ(){ return gnx; }

	int	 getNumParticleAt(int i, int j, int k);
	void markWater(char ***A, double ***halfwall, double density);
	void deleteAllParticles();

protected:
	std::vector<Particle *> ***cells;
	int gnx;
	int gny;
	int gnz;
};

#endif