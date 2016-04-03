//--------------------------------------------------------------------------------------------------------
//# FLIP_FLUID
//This is a flip fluid implement acooding to article <<3D Particle in Cell / Fluid Implicit Particle Fluid
//Solver using OpenMP directives>>.
//--------------------------------------------------------------------------------------------------------
#pragma once
#ifndef FLIP_H_
#define FLIP_H_

#include <windows.h>
#include <iostream>
#include <vector>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include "App.h"
#include "util.h"
#include "solver.h"
#include "Timer.h"
#include "Sorter.h"
#include "Mesher.h"

static ofstream debug("debug.txt");


static std::vector<double> vertices;
static std::vector<double> normals;
static std::vector<int> faces;
static double MaxDensity;

#define DRAW_PARTICLES
//#define OB							//you can decomment this line to see result of obstacle
//#define FACES							//whether do surface reconstruction
//#define USING_GV2ADVECT				//use grid velocity to advect particle (maybe this is wrong?)
//#define CORRECT						//post correct particle
//#define COMPUTE_DENSITY				//whether compute density per frame
class FlipFluid : public App
{
public:
	FlipFluid();
	FlipFluid(int w, int h);
	~FlipFluid();

	bool Init();
	void UpdateScene();
	void Rendering();
	void onResize(GLFWwindow* window, int w, int h);

	void onMouseWheel(GLFWwindow* window, double x, double y);
	void onMouseMove(GLFWwindow* window, double xd, double yd);
	void onMouseButton(GLFWwindow* window, int button, int action, int mods);
	void onKey(GLFWwindow* window, int key, int scancode, int action, int mods);

	void reset();
	void reset2();
	void DrawSolid();
	void DrawFaces();
	void pourWater();
	void DrawCube();
	void DrawGrid();
	void DrawParticle();
	void DrawVOfLocalP();

	//------------------------------------------------------------------------------
	//Traversal all particles and identify which cell include them.
	//------------------------------------------------------------------------------ 
	void ParCellIdentification();

	//------------------------------------------------------------------------------
	//Mark the cell to AIR or WATER or SOLID
	//------------------------------------------------------------------------------ 
	void FluidSurfaceMark();

	void addGravityF();

	//------------------------------------------------------------------------------
	//get neighbour particles of one cell (total 26 cells(except solid)) 
	//------------------------------------------------------------------------------ 
	void getNeighbourPar(int x, int y, int z, vector<Particle*>& v_Of_p, int parnum = -1);
	//------------------------------------------------------------------------------
	//get neighbour particles of one grid surface in X, Y, X direction.(total 18 cells(except solid))
	//------------------------------------------------------------------------------ 
	void getNeighbourX(int x, int y, int z, vector<Particle*>& v_Of_p);
	void getNeighbourY(int x, int y, int z, vector<Particle*>& v_Of_p);
	void getNeighbourZ(int x, int y, int z, vector<Particle*>& v_Of_p);
	//------------------------------------------------------------------------------
	//compute desity
	//------------------------------------------------------------------------------ 
	void computeDensity();
	//------------------------------------------------------------------------------
	//Set all boundary velocity to zero
	//------------------------------------------------------------------------------ 
	void setBoundaryVe();
	//------------------------------------------------------------------------------
	//map particle velocity to grid
	//------------------------------------------------------------------------------ 
	void P2G();
	//------------------------------------------------------------------------------
	//Tri-linear interpolation
	//------------------------------------------------------------------------------ 
	double trilinInterpolation(double ***q, double x, double y, double z, int w, int h, int d);
	//------------------------------------------------------------------------------
	//transform grid velocity to particles
	//------------------------------------------------------------------------------ 
	void G2P(double ***vx, double ***vy, double ***vz);
	void G2P_particle(Particle *p, double *v, double ***vx, double ***vy, double ***vz);
	//------------------------------------------------------------------------------
	//calculate the divergence of velocity in grid
	//¡°how much¡± the grid velocity is changing between two continuous grid cells
	//------------------------------------------------------------------------------ 
	double DivOfVelocity(int x, int y, int z);
	//------------------------------------------------------------------------------
	//sample a velocity for the grid surface near water surface.
	//when a surface is not water, not near water and is air or near solid.
	//------------------------------------------------------------------------------ 
	void velocityExtrapolation();
	//------------------------------------------------------------------------------
	//solve pressure
	//and advect grid velocity
	//------------------------------------------------------------------------------ 

	void solvePressureAndNewV();
	//------------------------------------------------------------------------------
	//reference to https://code.google.com/archive/p/flip3d/ ando`s filp code.
	//velocity resample
	//------------------------------------------------------------------------------
	void resample(Vec3 &p, Vec3 &u, double re);
	//------------------------------------------------------------------------------
	//reference to https://code.google.com/archive/p/flip3d/ ando`s filp code.
	//velocity correct
	//------------------------------------------------------------------------------
	void correct();
	//------------------------------------------------------------------------------
	//if a particle move into solid cell
	//then reposition this particle to grid(0,0,0) and sample a velocity
	//------------------------------------------------------------------------------
	void reposition();

	//------------------------------------------------------------------------------
	//simulate function.
	//------------------------------------------------------------------------------
	void simulate();
private:
	grid *_grid;
	std::vector<Particle> _particle;

	Solver *_Solver;
	sorter *_sorter;

};
#endif