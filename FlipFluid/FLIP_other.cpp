#include "FLIP.h"

void FlipFluid::resample(Vec3 &p, Vec3 &u, double re) {
	// Variables for Neighboring Particles
	double wsum = 0.0;
	double save[3] = { u.X, u.Y, u.Z };
	u.X = u.Y = u.Z = 0.0;

	// Gather Neighboring Particles
	vector<Particle*> v_Of_p;
	getNeighbourPar(p.X * _grid->xNum, p.Y * _grid->yNum, p.Z * _grid->zNum, v_Of_p);

	for (int i = 0; i < v_Of_p.size(); i++)
	{
		Particle *np = v_Of_p[i];

		double dist2 = length2(p, np->_p);

		double w = PMASS * sharp_kernel(dist2, re);

		u.X += w * np->_v.X;
		u.Y += w * np->_v.Y;
		u.Z += w * np->_v.Z;
		wsum += w;

	}
	if (wsum)
	{
		u.X /= wsum;
		u.Y /= wsum;
		u.Z /= wsum;
	}
	else
	{
		u.X = save[0];
		u.Y = save[1];
		u.Z = save[2];
	}
}
//------------------------------------------------------------------------------
//reference to https://code.google.com/archive/p/flip3d/ ando`s filp code.
//velocity correct
//------------------------------------------------------------------------------
void FlipFluid::correct()
{
	double re = DENSITY / 32;
	ParCellIdentification();

	for (int i = 0; i < _particle.size(); i++)
	{
		double spring[3] = { 0.0, 0.0, 0.0 };

		int x = fmax(0, fmin(_grid->xNum - 1, _particle[i]._p.X * _grid->xNum));
		int y = fmax(0, fmin(_grid->yNum - 1, _particle[i]._p.Y * _grid->yNum));
		int z = fmax(0, fmin(_grid->zNum - 1, _particle[i]._p.Z * _grid->zNum));

		vector<Particle*> v_Of_p;
		getNeighbourPar(x, y, z, v_Of_p);

		for (int j = 0; j < v_Of_p.size(); j++)
		{
			Particle *np = v_Of_p[j];

			double dist = length(_particle[i]._p, np->_p);
			double sq_dist = dist * dist;

			double w = SPRING * PMASS * smooth_kernel(sq_dist, re);
			if (sq_dist > 0.1* re)
			{
				spring[0] += w * (_particle[i]._p.X - np->_p.X) / dist * re;
				spring[1] += w * (_particle[i]._p.Y - np->_p.Y) / dist * re;
				spring[2] += w * (_particle[i]._p.Z - np->_p.Z) / dist * re;
			}
			else
			{
				spring[0] += 0.01*re / dt *(rand() % 101) / 100.0;
				spring[1] += 0.01*re / dt *(rand() % 101) / 100.0;
				spring[2] += 0.01*re / dt *(rand() % 101) / 100.0;
			}
		}
		_particle[i]._tempp.X = _particle[i]._p.X + dt*spring[0];
		_particle[i]._tempp.Y = _particle[i]._p.Y + dt*spring[1];
		_particle[i]._tempp.Z = _particle[i]._p.Z + dt*spring[2];
	}
	// Resample New Velocity
	for (int i = 0; i < _particle.size(); i++)
	{
		Particle *p = &_particle[i];
		p->_tempv.X = p->_v.X;
		p->_tempv.Y = p->_v.Y;
		p->_tempv.Z = p->_v.Z;
		resample(p->_tempp, p->_tempv, re);
	}

	// Update
	for (int i = 0; i < _particle.size(); i++)
	{
		Particle *p = &_particle[i];
		p->_p.X = p->_tempp.X;
		p->_p.Y = p->_tempp.Y;
		p->_p.Z = p->_tempp.Z;
		p->_v.X = p->_tempv.X;
		p->_v.Y = p->_tempv.Y;
		p->_v.Z = p->_tempv.Z;
	}
}
//------------------------------------------------------------------------------
//if a particle move into solid cell
//then reposition this particle to grid(0,0,0) and sample a velocity
//------------------------------------------------------------------------------
void FlipFluid::reposition()
{
	vector<uint> repoList;
	for (int i = 0; i < _particle.size(); i++)
	{
		int x = fmax(0, fmin(_grid->xNum - 1, _particle[i]._p.X * _grid->xNum));
		int y = fmax(0, fmin(_grid->yNum - 1, _particle[i]._p.Y * _grid->yNum));
		int z = fmax(0, fmin(_grid->zNum - 1, _particle[i]._p.Z * _grid->zNum));
		if (_grid->_gMark[x][y][z] == SOLID)
		{
#if 0
			debug << x << "  " << y << "  " << z << endl;
#endif
			repoList.push_back(i);
		}

	}
	if (repoList.empty()) return;

	// First Search for Deep Water
	vector<index> waters;
	while (waters.size() < repoList.size()) {
		FOR_EACH_CELL(_grid->xNum, _grid->yNum, _grid->zNum)
		{
#if 1
			if (i > 0 && _grid->_gMark[i - 1][j][k] != WATER) continue;
			if (i < _grid->xNum - 1 && _grid->_gMark[i + 1][j][k] != WATER) continue;
			if (j > 0 && _grid->_gMark[i][j - 1][k] != WATER) continue;
			if (j < _grid->yNum - 1 && _grid->_gMark[i][j + 1][k] != WATER) continue;
			if (k > 0 && _grid->_gMark[i][j][k - 1] != WATER) continue;
			if (k < _grid->zNum - 1 && _grid->_gMark[i][j][k + 1] != WATER) continue;
#endif
			if (_grid->_gMark[i][j][k] != WATER) continue;

			index aPos = { i, j, k };
			waters.push_back(aPos);
			if (waters.size() >= repoList.size())
			{
				i = _grid->xNum; j = _grid->yNum; k = _grid->zNum;
			}
		}
		if (waters.empty()) return;
	}

	// Shuffle
	random_shuffle(waters.begin(), waters.end());

	double h = 1.0 / GRIDN;
	for (int i = 0; i < repoList.size(); i++)
	{
		Particle &p = _particle[repoList[i]];
		p._p.X = h*(waters[i].X + 0.5 + 0.5*(rand() % 101) / 100);
		p._p.Y = h*(waters[i].Y + 0.5 + 0.5*(rand() % 101) / 100);
		p._p.Z = h*(waters[i].Z + 0.5 + 0.5*(rand() % 101) / 100);
#if 0
		debug << waters[i].X << "  " << waters[i].Y << "  " << waters[i].Z << endl;
#endif
	}

	ParCellIdentification();

	for (int i = 0; i < repoList.size(); i++)
	{
		Particle &p = _particle[repoList[i]];
		Velocity u = { 0.0, 0.0, 0.0 };
		resample(p._p, u, h);
		p._v = u;
	}
}