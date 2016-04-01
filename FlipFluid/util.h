#ifndef UTIL_H
#define UTIL_H
//------------------------------------------------------------------------------
//Inluding macro defination, structs,utility functions, grid defination.
//------------------------------------------------------------------------------
#include <algorithm>
#include <vector>
#define FOR_EACH_CELL(I, J, K)  for(int i = 0; i < I; i++)for(int j = 0; j < J; j++)for(int k = 0; k < K; k++)
#define FOR_EACH_CELL_FLIP(I, J, K)  for(int i = I-1; i >= 0; i--)for(int j = J-1; j >= 0; j--)for(int k = K-1; k >= 0; k--)

#define GRIDN  32
#define DENSITY 0.5
#define dt     0.01666667
#define PMASS   1.0
#define SOLID   0x0001
#define AIR   0x0002
#define WATER   0x0003
#define VCOEFFICIENT 0.05
#define SPRING 50

//#define USING_GV2ADVECT

typedef unsigned int uint;

// struct for index of grid.
typedef struct
{
	uint X;
	uint Y;
	uint Z;
}index;
// struct for index of vector3, position and velocity of particle.
typedef struct
{
	double X;
	double Y;
	double Z;
}Vec3, Position, Velocity;

//square of length of two vector3
double length2(Vec3& a, Vec3& b)
{
	return (a.X - b.X)*(a.X - b.X) + (a.Y - b.Y) *(a.Y - b.Y) + (a.Z - b.Z)*(a.Z - b.Z);
}
//length of two vector3
double length(Vec3& a, Vec3& b)
{
	return sqrt(length2(a, b));
}
//length of two double[3]
double length(double p0[3], double p1[3]) {
	return hypotf(hypotf(p0[0] - p1[0], p0[1] - p1[1]), p0[2] - p1[2]);
}

double hypot2(double a, double b, double c) {
	return a*a + b*b + c*c;
}
//square of length of two double[3]
double length2(double p0[3], double p1[3]) {
	return hypot2(p0[0] - p1[0], p0[1] - p1[1], p0[2] - p1[2]);
}

//Smooth Kernel Functions
double smooth_kernel(double d, double h)
{
	return fmax(1.0 - d / (h*h), 0.0);
}
double sharp_kernel(double d, double h)
{
	return fmax((h*h) / fmax(d, 1.0e-5) - 1.0, 0.0);
}
//copy a to b
void mycopy(double ***src, double ***des, int x, int y, int z)
{
	FOR_EACH_CELL(x, y, z)
	{
		des[i][j][k] = src[i][j][k];
	}
}
//allocate memory
template<class T>
T*** memAlloc3D(int x, int y, int z)
{
	T ***mem = new T**[x];
	for (int i = 0; i < x; i++)
	{
		mem[i] = new T*[y];
		for (int j = 0; j < y; j++)
		{
			mem[i][j] = new T[z];
		}
	}
	return mem;
}
//grid defination
struct grid
{
	int xNum;
	int yNum;
	int zNum;
	double deltaH;

	double ***_vx;				//grid x velocity
	double ***_vy;				
	double ***_vz;				
	double ***_vx_save;			//grid temporary x velocity
	double ***_vy_save;
	double ***_vz_save;

	double ***_pressure;		//pressure between grids
	double ***_div;				//divergence between grids
	char ***_gMark;
	std::vector<int> ***_particleInCell;    //indicate which particle is in _particleInCell[*] 
	grid()
	{
		xNum = GRIDN;
		yNum = GRIDN;
		zNum = GRIDN;
		deltaH = 1.0 / GRIDN;

		_vx = memAlloc3D<double>(xNum + 1, yNum, zNum);
		_vy = memAlloc3D<double>(xNum, yNum + 1, zNum);
		_vz = memAlloc3D<double>(xNum, yNum, zNum + 1);
		_vx_save = memAlloc3D<double>(xNum + 1, yNum, zNum);
		_vy_save = memAlloc3D<double>(xNum, yNum + 1, zNum);
		_vz_save = memAlloc3D<double>(xNum, yNum, zNum + 1);

		_pressure = memAlloc3D<double>(xNum, yNum, zNum);
		_div = memAlloc3D<double>(xNum, yNum, zNum);
		_gMark = memAlloc3D<char>(xNum, yNum, zNum);

		_particleInCell = memAlloc3D<std::vector<int>>(xNum, yNum, zNum);
	}
	void reset()
	{
		FOR_EACH_CELL(xNum, yNum, zNum)
		{
			_pressure[i][j][k] = 0.0;
			_div[i][j][k] = 0.0;
			_gMark[i][j][k] = AIR;
		}
		FOR_EACH_CELL(xNum + 1, yNum, zNum)
		{
			_vx[i][j][k] = _vx_save[i][j][k] = 0.0;
		}
		FOR_EACH_CELL(xNum, yNum + 1, zNum)
		{
			_vy[i][j][k] = _vy_save[i][j][k] = 0.0;
		}
		FOR_EACH_CELL(xNum, yNum, zNum + 1)
		{
			_vz[i][j][k] = _vz_save[i][j][k] = 0.0;
		}
	}
};
//Particle defination
struct Particle
{
	Position _p;				//particle temporary position
	Position _tempp;			//particle temporary velocity
	Velocity _v;				//particle temporary position
	Velocity _tempv;			//particle temporary velocity
	double density;
	Particle(double px, double py, double pz, double vx, double vy, double vz)
	{
		_p.X = px;
		_p.Y = py;
		_p.Z = pz;
		_v.X = vx;
		_v.Y = vy;
		_v.Z = vz;
		density = DENSITY;
	}
};
#endif
