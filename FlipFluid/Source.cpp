//--------------------------------------------------------------------------------------------------------
//# FLIP_FLUID
//This is a flip fluid implement acooding to article <<3D Particle in Cell / Fluid Implicit Particle Fluid
//Solver using OpenMP directives>>.
//--------------------------------------------------------------------------------------------------------
#include "FLIP.h"

int main(void)
{
	FlipFluid *theApp = new FlipFluid(800, 600);

	if (!theApp->Init())
		return 0;
	theApp->Run();
	delete theApp;
	return 0;
}