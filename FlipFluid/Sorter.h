/*
*  sorter.cpp
*  flip3D
*/


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


sorter::sorter(int _gnx, int _gny, int _gnz) {
	gnx = _gnx;
	gny = _gny;
	gnz = _gnz;
	cells = memAlloc3D<vector<Particle *> >(gnx, gny, gnz);
}

sorter::~sorter() {
}

void sorter::sort(std::vector<Particle> &_Particles) {
	// Clear All Cells
	FOR_EACH_CELL(gnx, gny, gnz) 
	{
		cells[i][j][k].clear();
	}

		// Store Into The Cells
		for (int n = 0; n<_Particles.size(); n++)
		{
			Particle *p = &(_Particles[n]);
			double pos[3];

			pos[0] = p->_p.X;
			pos[1] = p->_p.Y;
			pos[2] = p->_p.Z;

			int i = fmax(0, fmin(gnx - 1, gnx*pos[0]));
			int j = fmax(0, fmin(gny - 1, gny*pos[1]));
			int k = fmax(0, fmin(gnz - 1, gnz*pos[2]));
			cells[i][j][k].push_back(p);
		}
#if 0
		FOR_EACH_CELL(gnx, gny, gnz)
		{
			for (int l = 0; l < cells[i][j][k].size(); l++)
			{
				debug << i << "  " << j << "  " << k << "  " << cells[i][j][k][l] << endl;
			}
			
		}
		debug << endl;
#endif
}

std::vector<Particle *> sorter::getNeigboringParticles_wall(int i, int j, int k, int w, int h, int d) {
	std::vector<Particle *> res;
	for (int si = i - w; si <= i + w - 1; si++) for (int sj = j - h; sj <= j + h - 1; sj++) for (int sk = k - d; sk <= k + d - 1; sk++) {
		if (si < 0 || si > gnx - 1 || sj < 0 || sj > gny - 1 || sk < 0 || sk > gnz - 1) continue;
		for (int a = 0; a<cells[si][sj][sk].size(); a++) {
			Particle *p = cells[si][sj][sk][a];
			res.push_back(p);
		}
	}
	return res;
}

std::vector<Particle *> sorter::getNeigboringParticles_cell(int i, int j, int k, int w, int h, int d) {
	std::vector<Particle *> res;
	for (int si = i - w; si <= i + w; si++) for (int sj = j - h; sj <= j + h; sj++) for (int sk = k - d; sk <= k + d; sk++) {
		if (si < 0 || si > gnx - 1 || sj < 0 || sj > gny - 1 || sk < 0 || sk > gnz - 1) continue;
		for (int a = 0; a<cells[si][sj][sk].size(); a++) {
			Particle *p = cells[si][sj][sk][a];
			res.push_back(p);
		}
	}
#if 0
	debug << i << "  " << j << "  " << k << "  " << res.size() <<endl;
#endif
	return res;
}

int	 sorter::getNumParticleAt(int i, int j, int k) {
	return cells[i][j][k].size();
}

double sorter::levelset(int i, int j, int k, double ***halfwall, double density) {
	double accm = 0.0;
	for (int a = 0; a<cells[i][j][k].size(); a++)
	{
		//if (cells[i][j][k][a]->type == FLUID) {
			accm += cells[i][j][k][a]->density;

		//else {
			//return 1.0;
		//}
	}
	double n0 = 1.0 / (density*density*density);
	return 0.2*n0 - accm;
}

void sorter::markWater(char ***A, double ***halfwall, double density) {
	FOR_EACH_CELL(gnx, gny, gnz) 
	{
		A[i][j][k] = AIR;
		//for (int a = 0; a<cells[i][j][k].size(); a++) {
		//	if (cells[i][j][k][a]->type == WALL) {
		//		A[i][j][k] = WALL;
		//	}
		//}
		if (A[i][j][k] != SOLID) A[i][j][k] = levelset(i, j, k, halfwall, density) < 0.0 ? WATER : AIR;
	}
}

void sorter::deleteAllParticles() {
	FOR_EACH_CELL(gnx, gny, gnz) 
	{
		for (int a = 0; a<cells[i][j][k].size(); a++) {
			delete cells[i][j][k][a];
		}
	}
}
#endif