#include "FLIP.h"

int oldX = 0, oldY = 0;
float rX = 0.0, rY = 0.0;
bool isMouseButton = false;

int drawOB = 1;


FlipFluid::FlipFluid() :App()
{
}
FlipFluid::FlipFluid(int w, int h) : App(w, h)
{

}
FlipFluid::~FlipFluid()
{

}

void FlipFluid::reset()
{
	int x = _grid->xNum / 3;
	int y = _grid->yNum;
	int z = _grid->zNum / 2;
	double dH = _grid->deltaH;
	double dH_4 = dH / 4;

	_particle.clear();
	for (int i = 0 + 8; i < x + 10; i++)
		for (int j = 8; j < y - 5; j++)
			for (int k = 2; k < z; k++)
			{
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));


				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
			}
	for (int i = 0; i < _grid->xNum; i++)
		for (int j = 0; j < _grid->yNum; j++)
			for (int k = 0; k < 2; k++)
			{
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));


				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
			}
	_grid->reset();
}
void FlipFluid::reset2()
{
	int x = _grid->xNum / 3;
	int y = _grid->yNum;
	int z = _grid->zNum / 2;
	double dH = _grid->deltaH;
	double dH_4 = dH / 4;

	_particle.clear();
	for (int i = 0; i < x; i++)
		for (int j = y - _grid->yNum / 3; j < y; j++)
			for (int k = 2; k < z; k++)
			{
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));


				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
			}
	for (int i = _grid->xNum - x; i < _grid->xNum; i++)
		for (int j = 0; j < _grid->yNum / 3; j++)
			for (int k = 2; k < z; k++)
			{
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));


				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
			}
	for (int i = 0; i < _grid->xNum; i++)
		for (int j = 0; j < _grid->yNum; j++)
			for (int k = 0; k < 2; k++)
			{
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));


				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
				_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
			}
	_grid->reset();
}
void FlipFluid::DrawSolid()
{
	int x = 26;
	int y = 12;
	int z = 0;
	int xn = 1;
	int yn = 8;
	int zn = 10;
	double dH = _grid->deltaH;
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	glBegin(GL_QUADS);
	//front
	glColor3f(0.5f, 0.0f, 0.0f);
	glVertex3d(x * dH, y * dH, 0.0f);
	glVertex3d(x * dH + dH, y * dH, 0.0f);
	glVertex3d(x * dH + dH, y * dH, zn * dH + dH);
	glVertex3d(x * dH, y * dH, zn * dH + dH);
	//left
	glColor3f(0.0f, 0.5f, 0.0f);
	glVertex3d(x * dH, y * dH, 0.0f);
	glVertex3d(x * dH, y * dH, zn * dH + dH);
	glVertex3d(x * dH, (y + yn + 1) * dH, zn * dH + dH);
	glVertex3d(x * dH, (y + yn + 1) * dH, 0.0);
	//right
	glColor3f(0.0f, 0.0f, 0.5f);
	glVertex3d(x * dH + dH, y * dH, zn * dH + dH);
	glVertex3d(x * dH + dH, y * dH, 0.0f);
	glVertex3d(x * dH + dH, (y + yn + 1) * dH, 0.0f);
	glVertex3d(x * dH + dH, (y + yn + 1) * dH, zn * dH + dH);
	//top
	glColor3f(0.5f, 0.5f, 0.5f);
	glVertex3d(x * dH, y * dH, zn * dH + dH);
	glVertex3d(x * dH + dH, y * dH, zn * dH + dH);
	glVertex3d(x * dH + dH, (y + yn + 1) * dH, zn * dH + dH);
	glVertex3d(x * dH, (y + yn + 1) * dH, zn * dH + dH);
	//back
	glColor3f(0.5f, 0.5f, 0.0f);
	glVertex3d(x * dH + dH, (y + yn + 1) * dH, zn * dH + dH);
	glVertex3d(x * dH + dH, (y + yn + 1) * dH, 0.0f);
	glVertex3d(x * dH, (y + yn + 1) * dH, 0.0f);
	glVertex3d(x * dH, (y + yn + 1) * dH, zn * dH + dH);
	glEnd();
	glDisable(GL_DEPTH_TEST);
	//FOR_EACH_CELL(_grid->xNum, _grid->yNum, _grid->zNum)
	//{
	//	if (_grid->_gMark[i][j][k] = SOLID;)
	//}
}
void FlipFluid::DrawFaces()
{
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glColor3f(0.8f, 0.8f, 0.8f);
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < faces.size() / 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			uint index = faces[i * 3 + j];
			double posX = vertices[index * 3];
			double posY = vertices[index * 3 + 1];
			double posZ = vertices[index * 3 + 2];
			glNormal3d(normals[index * 3], normals[index * 3 + 1], normals[index * 3 + 2]);
			glVertex3d(posX, posY, posZ);
#if 0
			debug << posX << "  " << posY << "  " << posZ << endl;
#endif
		}
	}
	//debug << "=======================" << endl;
	glEnd();
	glDisable(GL_DEPTH_TEST);
}
void FlipFluid::pourWater()
{
	int i = 0;
	double dH = _grid->deltaH;
	double dH_4 = dH / 4;
	for (int j = 4; j < _grid->yNum / 2; j++)
		for (int k = _grid->zNum - 10; k < _grid->zNum - 4; k++)
		{
			_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4, k*dH + dH_4, 10, 0, -2));
			_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4, 10, 0, -2));
			_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4, k*dH + dH_4 + dH_4 + dH_4, 10, 0, -2));
			_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4 + dH_4 + dH_4, 10, 0, -2));
		}
}
void FlipFluid::DrawCube()
{
	glColor3f(0.5f, 0.5f, 0.5f);
	//glBegin(GL_TRIANGLES);
	//	glVertex3f(0.0f, 0.0f, 0.0f);
	//	glVertex3f(1.f, 0.f, 0.0f);
	//	glVertex3f(0.0f, 1.f, 0.0f);
	//glEnd();
	glBegin(GL_LINES);
	//cout << _grid->deltaX << " " <<_grid->width <<endl;
	//x axis
	glColor3f(0.5f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(1.f, 0.f, 0.0f);
	//y axis
	glColor3f(0.0f, 0.5f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.f, 1.f, 0.0f);
	//zaxis
	glColor3f(0.0f, 0.0f, 0.5f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.f, 0.f, 1.0f);
	//set to gray
	glColor3f(0.5f, 0.5f, 0.5f);

	glVertex3f(1.0f, 1.0f, 0.0f);
	glVertex3f(1.f, 0.f, 0.0f);

	glVertex3f(1.0f, 1.0f, 0.0f);
	glVertex3f(0.f, 1.f, 0.0f);

	glVertex3f(1.0f, 1.0f, 0.0f);
	glVertex3f(1.f, 1.f, 1.0f);

	glVertex3f(1.0f, 0.0f, 1.0f);
	glVertex3f(1.0f, 0.0f, 0.0f);

	glVertex3f(1.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 1.0f);

	glVertex3f(1.0f, 0.0f, 1.0f);
	glVertex3f(1.0f, 1.0f, 1.0f);

	glVertex3f(0.0f, 1.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 1.0f);

	glVertex3f(0.0f, 1.0f, 1.0f);
	glVertex3f(0.0f, 1.0f, 0.0f);

	glVertex3f(0.0f, 1.0f, 1.0f);
	glVertex3f(1.0f, 1.0f, 1.0f);
	glEnd();
}
void FlipFluid::DrawGrid()
{
	glColor3f(0.5f, 0.5f, 0.5f);
	glBegin(GL_LINES);
	double dH = _grid->deltaH;
	FOR_EACH_CELL(_grid->xNum + 1, _grid->yNum + 1, _grid->zNum + 1)
	{
		//x
		glVertex3d(0, j * dH, k * dH);
		glVertex3d(1, j * dH, k * dH);
		//y		
		glVertex3d(i * dH, 0, k * dH);
		glVertex3d(i * dH, 1, k * dH);
		//z		 
		glVertex3d(i * dH, j * dH, 0);
		glVertex3d(i * dH, j * dH, 1);
	}
	glEnd();
}
void FlipFluid::DrawVOfLocalP()
{
	glColor3f(1.0f, 1.0f, 0.0f);

	glBegin(GL_LINES);
	FOR_EACH_CELL(_grid->xNum, _grid->yNum, _grid->zNum)
	{
		//if (_grid->_gMark[i][j][k] == WATER &&_grid->_gMark[i + 1][j][k] == AIR)
		if ((i == 11 || i == 10 || i == 12) && 10 < j&&j < 16 && k == 3 && _grid->_particleInCell[i][j][k].size())
		{
			Particle* p = &(_particle[_grid->_particleInCell[i][j][k][0]]);
			glVertex3d(p->_p.X, p->_p.Y, p->_p.Z);
			glVertex3d(p->_p.X + p->_v.X, p->_p.Y + p->_v.Y, p->_p.Z + p->_v.Z);
#if 0
			debug << p->_v.X << "  " << p->_v.Y << "  " << p->_v.Z << endl;
#endif
		}
	}
	//debug << "==================================" << endl;
	glEnd();
}
void FlipFluid::DrawParticle()
{
	glColor3f(1.0f, 1.0f, 1.0f);
	glPointSize(1.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < _particle.size(); i++)
	{
		//cout << _particle[i]._p.X << " " << _particle[i]._p.Y << " "<<_particle[i]._p.Z << endl;
		glVertex3d(_particle[i]._p.X, _particle[i]._p.Y, _particle[i]._p.Z);
	}
	glEnd();
}


//------------------------------------------------------------------------------
//Traversal all particles and identify which cell include them.
//------------------------------------------------------------------------------ 
void FlipFluid::ParCellIdentification()
{
	for (int i = 0; i < _grid->xNum; i++)
		for (int j = 0; j < _grid->yNum; j++)
			for (int k = 0; k < _grid->zNum; k++)
			{
				_grid->_particleInCell[i][j][k].clear();
			}
	for (int i = 0; i < _particle.size(); i++)
	{
		int x = fmax(0, fmin(_grid->xNum - 1, _particle[i]._p.X * _grid->xNum));			   //position * grid number
		int y = fmax(0, fmin(_grid->yNum - 1, _particle[i]._p.Y * _grid->yNum));			   //
		int z = fmax(0, fmin(_grid->zNum - 1, _particle[i]._p.Z * _grid->zNum));			   //
#if 0
		debug << x << " " << y << " " << z << endl;
#endif
		_grid->_particleInCell[x][y][z].push_back(i);
	}
}

//------------------------------------------------------------------------------
//Mark the cell to AIR or WATER or SOLID
//------------------------------------------------------------------------------ 
void FlipFluid::FluidSurfaceMark()
{
	double alphaNLo = 0.2 * 8 * DENSITY;
	//set all mark to air
	FOR_EACH_CELL(_grid->xNum, _grid->yNum, _grid->zNum)
	{
		_grid->_gMark[i][j][k] = AIR;
	}
	for (int i = 0; i < _grid->xNum; i++)
		for (int j = 0; j < _grid->yNum; j++)
			for (int k = 0; k < _grid->zNum; k++)
			{
				//have a particle mark to water
				//if (_grid->_particleInCell[i][j][k].size())
				//	_grid->_gMark[i][j][k] = WATER;

				double densitySum = 0.0;
				for (int l = 0; l < _grid->_particleInCell[i][j][k].size(); l++)
				{
					densitySum += _particle[_grid->_particleInCell[i][j][k][l]].density;
				}
				//density in a cell is bigger than alphaNlo makt to WATER
				if (alphaNLo < densitySum)
				{
					_grid->_gMark[i][j][k] = WATER;
				}

			}
	//mark obstacle cell to SOLID.
#ifdef OB
	int x = 26;
	for (int y = 12; y < 20; y++)
		for (int z = 0; z < 10; z++)
			_grid->_gMark[x][y][z] = SOLID;
#endif
}

void FlipFluid::addGravityF()
{
	for (int i = 0; i < _particle.size(); i++)
	{
		_particle[i]._v.Z -= 9.81 * dt;
	}
}

//------------------------------------------------------------------------------
//get neighbour particles of one cell (total 26 cells(except solid)) 
//------------------------------------------------------------------------------ 
void FlipFluid::getNeighbourPar(int x, int y, int z, vector<Particle*>& v_Of_p, int parnum)
{
	v_Of_p.clear();
	for (int i = x - 1; i <= x + 1; i++)
		for (int j = y - 1; j <= y + 1; j++)
			for (int k = z - 1; k <= z + 1; k++)
			{
				if (i < 0 || i > _grid->xNum - 1 || j < 0 || j > _grid->yNum - 1 || k < 0 || k > _grid->zNum - 1)
					continue;

				for (int l = 0; l < _grid->_particleInCell[i][j][k].size(); l++)
				{

					if (l != parnum)
						v_Of_p.push_back(&(_particle[_grid->_particleInCell[i][j][k][l]]));
				}
			}
}
//------------------------------------------------------------------------------
//get neighbour particles of one grid surface in X, Y, X direction.(total 18 cells(except solid))
//------------------------------------------------------------------------------ 
void FlipFluid::getNeighbourX(int x, int y, int z, vector<Particle*>& v_Of_p)
{
	v_Of_p.clear();
	for (int i = x - 1; i <= x; i++)
		for (int j = y - 1; j <= y + 1; j++)
			for (int k = z - 1; k <= z + 1; k++)
			{
				if (i < 0 || i > _grid->xNum - 2 || j < 0 || j > _grid->yNum - 1 || k < 0 || k > _grid->zNum - 1)
					continue;

				for (int l = 0; l < _grid->_particleInCell[i][j][k].size(); l++)
				{
					v_Of_p.push_back(&(_particle[_grid->_particleInCell[i][j][k][l]]));
				}
			}
}
void FlipFluid::getNeighbourY(int x, int y, int z, vector<Particle*>& v_Of_p)
{
	v_Of_p.clear();
	for (int i = x - 1; i <= x + 1; i++)
		for (int j = y - 1; j <= y; j++)
			for (int k = z - 1; k <= z + 1; k++)
			{
				if (i < 0 || i > _grid->xNum - 1 || j < 0 || j > _grid->yNum - 2 || k < 0 || k > _grid->zNum - 1)
					continue;

				for (int l = 0; l < _grid->_particleInCell[i][j][k].size(); l++)
				{
					v_Of_p.push_back(&(_particle[_grid->_particleInCell[i][j][k][l]]));
				}
			}
}
void FlipFluid::getNeighbourZ(int x, int y, int z, vector<Particle*>& v_Of_p)
{
	v_Of_p.clear();
	for (int i = x - 1; i <= x + 1; i++)
		for (int j = y - 1; j <= y + 1; j++)
			for (int k = z - 1; k <= z; k++)
			{
				if (i < 0 || i > _grid->xNum - 1 || j < 0 || j > _grid->yNum - 1 || k < 0 || k > _grid->zNum - 2)
					continue;

				for (int l = 0; l < _grid->_particleInCell[i][j][k].size(); l++)
				{
					v_Of_p.push_back(&(_particle[_grid->_particleInCell[i][j][k][l]]));
				}
			}
}
//------------------------------------------------------------------------------
//compute desity
//------------------------------------------------------------------------------ 
void FlipFluid::computeDensity()
{
	for (int i = 0; i < _particle.size(); i++)
	{
		vector<Particle*> v_Of_p;
		getNeighbourPar(_particle[i]._p.X * _grid->xNum, _particle[i]._p.Y * _grid->yNum, _particle[i]._p.Z * _grid->zNum, v_Of_p, i);

		double sum = 0.0;
		for (int j = 0; j < v_Of_p.size(); j++)
		{
			double w = smooth_kernel(length2(_particle[i]._p, v_Of_p[j]->_p), _grid->deltaH * 2.0);
			sum += PMASS * w;
		}
		_particle[i].density = sum / MaxDensity;

	}
}
//------------------------------------------------------------------------------
//Set all boundary velocity to zero
//------------------------------------------------------------------------------ 
void FlipFluid::setBoundaryVe()
{
	//set x
	int i = 0;
	for (int j = 0; j < _grid->yNum; j++)
		for (int k = 0; k < _grid->zNum; k++)
			_grid->_vx[i][j][k] = 0.0;

	i = _grid->xNum;
	for (int j = 0; j < _grid->yNum; j++)
		for (int k = 0; k < _grid->zNum; k++)
			_grid->_vx[i][j][k] = 0.0;
	//set y
	i = 0;
	for (int j = 0; j < _grid->xNum; j++)
		for (int k = 0; k < _grid->zNum; k++)
			_grid->_vy[j][i][k] = 0.0;

	i = _grid->yNum;
	for (int j = 0; j < _grid->xNum; j++)
		for (int k = 0; k < _grid->zNum; k++)
			_grid->_vy[j][i][k] = 0.0;

	//set z
	i = 0;
	for (int j = 0; j < _grid->xNum; j++)
		for (int k = 0; k < _grid->yNum; k++)
			_grid->_vz[j][k][i] = 0.0;

	i = _grid->zNum;
	for (int j = 0; j < _grid->xNum; j++)
		for (int k = 0; k < _grid->yNum; k++)
			_grid->_vz[j][k][i] = 0.0;
#ifdef OB
	int x = 26;
	for (int y = 12; y < 20; y++)
		for (int z = 0; z < 10; z++)
		{
			_grid->_vx[x][y][z] = 0.0;
			_grid->_vx[x + 1][y][z] = 0.0;
			_grid->_vy[x][y][z] = 0.0;
			_grid->_vy[x][y + 1][z] = 0.0;
			_grid->_vz[x][y][z] = 0.0;
			_grid->_vz[x][y][z + 1] = 0.0;
		}
#endif
}
//------------------------------------------------------------------------------
//map particle velocity to grid
//------------------------------------------------------------------------------ 
void FlipFluid::P2G()
{
	vector<Particle*> neighbours;
	double dH = _grid->deltaH;
	FOR_EACH_CELL(_grid->xNum + 1, _grid->yNum + 1, _grid->zNum + 1)
	{
		//map x
		if (j < _grid->yNum && k < _grid->zNum)
		{
			Vec3 pos = { dH * i, dH * j + 0.5 * dH, dH * k + 0.5 * dH };
			getNeighbourX(i, j, k, neighbours);													//get neighbour particles
			double sumW = 0;
			double sumX = 0;
			for (int l = 0; l < neighbours.size(); l++)
			{
				double w = sharp_kernel(length2(pos, neighbours[l]->_p), 1.4);					//use stiff kernel weight to sample
				sumW += w;
				sumX += w * neighbours[l]->_v.X;
			}
			_grid->_vx[i][j][k] = sumW ? sumX / sumW : 0.0;
		}
		//map y
		if (i < _grid->xNum && k < _grid->zNum)
		{
			Vec3 pos = { dH * i + 0.5*dH, dH * j, dH * k + 0.5 * dH };
			getNeighbourY(i, j, k, neighbours);
			double sumW = 0;
			double sumX = 0;
			for (int l = 0; l < neighbours.size(); l++)
			{
				double w = sharp_kernel(length2(pos, neighbours[l]->_p), 1.4);
				sumW += w;
				sumX += w * neighbours[l]->_v.Y;
			}
			_grid->_vy[i][j][k] = sumW ? sumX / sumW : 0.0;
		}
		//map z
		if (i < _grid->xNum && j < _grid->yNum)
		{
			Vec3 pos = { dH * i + 0.5 * dH, dH * j + 0.5 * dH, dH * k };
			getNeighbourZ(i, j, k, neighbours);
			double sumW = 0;
			double sumX = 0;
			for (int l = 0; l < neighbours.size(); l++)
			{
				double w = sharp_kernel(length2(pos, neighbours[l]->_p), 1.4);
				sumW += w;
				sumX += w * neighbours[l]->_v.Z;
			}
			_grid->_vz[i][j][k] = sumW ? sumX / sumW : 0.0;
		}
	}
}
//------------------------------------------------------------------------------
//Tri-linear interpolation
//------------------------------------------------------------------------------ 
double FlipFluid::trilinInterpolation(double ***q, double x, double y, double z, int w, int h, int d)
{
	x = fmax(0.0, fmin(w, x));
	y = fmax(0.0, fmin(h, y));
	z = fmax(0.0, fmin(d, z));
	int i = min(x, w - 2);
	int j = min(y, h - 2);
	int k = min(z, d - 2);

	return	(k + 1 - z)*(((i + 1 - x)*q[i][j][k] + (x - i)*q[i + 1][j][k])*(j + 1 - y) + ((i + 1 - x)*q[i][j + 1][k] + (x - i)*q[i + 1][j + 1][k])*(y - j)) +
		(z - k)*(((i + 1 - x)*q[i][j][k + 1] + (x - i)*q[i + 1][j][k + 1])*(j + 1 - y) + ((i + 1 - x)*q[i][j + 1][k + 1] + (x - i)*q[i + 1][j + 1][k + 1])*(y - j));
}
//------------------------------------------------------------------------------
//transform grid velocity to particles
//------------------------------------------------------------------------------ 
void FlipFluid::G2P(double ***vx, double ***vy, double ***vz)
{
	for (int i = 0; i < _particle.size(); i++)
	{
		double x = _particle[i]._p.X;
		double y = _particle[i]._p.Y;
		double z = _particle[i]._p.Z;

		_particle[i]._v.X = trilinInterpolation(vx, x * _grid->xNum, y * _grid->yNum - 0.5, z * _grid->zNum - 0.5, _grid->xNum + 1, _grid->yNum, _grid->xNum);
		_particle[i]._v.Y = trilinInterpolation(vy, x * _grid->xNum - 0.5, y * _grid->yNum, z * _grid->zNum - 0.5, _grid->xNum, _grid->yNum + 1, _grid->xNum);
		_particle[i]._v.Z = trilinInterpolation(vz, x * _grid->xNum - 0.5, y * _grid->yNum - 0.5, z * _grid->zNum, _grid->xNum, _grid->yNum, _grid->xNum + 1);
	}
}
void FlipFluid::G2P_particle(Particle *p, double *v, double ***vx, double ***vy, double ***vz)
{
	double x = p->_p.X;
	double y = p->_p.Y;
	double z = p->_p.Z;

	v[0] = trilinInterpolation(vx, x * _grid->xNum, y * _grid->yNum - 0.5, z * _grid->zNum - 0.5, _grid->xNum + 1, _grid->yNum, _grid->xNum);
	v[1] = trilinInterpolation(vy, x * _grid->xNum - 0.5, y * _grid->yNum, z * _grid->zNum - 0.5, _grid->xNum, _grid->yNum + 1, _grid->xNum);
	v[2] = trilinInterpolation(vz, x * _grid->xNum - 0.5, y * _grid->yNum - 0.5, z * _grid->zNum, _grid->xNum, _grid->yNum, _grid->xNum + 1);
}
//------------------------------------------------------------------------------
//calculate the divergence of velocity in grid
//¡°how much¡± the grid velocity is changing between two continuous grid cells
//------------------------------------------------------------------------------ 
double FlipFluid::DivOfVelocity(int x, int y, int z)
{
	double uxleft;
	double uxright;
	double uyleft;
	double uyright;
	double uzleft;
	double uzright;
	double dH = _grid->deltaH;
	uxleft = _grid->_vx[x][y][z];
	uxright = _grid->_vx[x + 1][y][z];

	uyleft = _grid->_vy[x][y][z];
	uyright = _grid->_vy[x][y + 1][z];

	uzleft = _grid->_vz[x][y][z];
	uzright = _grid->_vz[x][y][z + 1];
#if 0
	debug << "Grid: " << x << " " << y << " " << z << endl;
	debug << uxright << " " << uxleft << " " << uyright << " " << uyleft << " " << uzright << " " << uzleft << endl;
	debug << (uxright - uxleft) / dH + (uyright - uyleft) / dH + (uzright - uzleft) / dH << endl;
#endif
	return ((uxright - uxleft) + (uyright - uyleft) + (uzright - uzleft)) / dH;
}
//------------------------------------------------------------------------------
//sample a velocity for the grid surface near water surface.
//when a surface is not water, not near water and is air or near solid.
//------------------------------------------------------------------------------ 
void FlipFluid::velocityExtrapolation()
{//somthing wrong
	//mark face
	static	char ***faceMarkX[2] = { memAlloc3D<char>(_grid->xNum + 1, _grid->yNum, _grid->zNum), memAlloc3D<char>(_grid->xNum + 1, _grid->yNum, _grid->zNum) };
	static	char ***faceMarkY[2] = { memAlloc3D<char>(_grid->xNum, _grid->yNum + 1, _grid->zNum), memAlloc3D<char>(_grid->xNum, _grid->yNum + 1, _grid->zNum) };
	static	char ***faceMarkZ[2] = { memAlloc3D<char>(_grid->xNum, _grid->yNum, _grid->zNum + 1), memAlloc3D<char>(_grid->xNum, _grid->yNum, _grid->zNum + 1) };
	//x
	for (int i = 0; i < _grid->xNum + 1; i++)
		for (int j = 0; j < _grid->yNum; j++)
			for (int k = 0; k < _grid->zNum; k++)
			{
				//fluid or near fluid
				faceMarkX[0][i][j][k] = (i > 0 && _grid->_gMark[i - 1][j][k] == WATER) || (i < _grid->xNum && _grid->_gMark[i][j][k] == WATER);
				//near solid or air
				faceMarkX[1][i][j][k] = (i <= 0 || i >= _grid->xNum) || _grid->_gMark[i - 1][j][k] == AIR || _grid->_gMark[i][j][k] == SOLID;
			}
	//y
	for (int i = 0; i < _grid->xNum; i++)
		for (int j = 0; j < _grid->yNum + 1; j++)
			for (int k = 0; k < _grid->zNum; k++)
			{
				//fluid or near fluid
				faceMarkY[0][i][j][k] = (j > 0 && _grid->_gMark[i][j - 1][k] == WATER) || (j < _grid->yNum && _grid->_gMark[i][j][k] == WATER);
				//near solid or air
				faceMarkY[1][i][j][k] = (j <= 0 || j >= _grid->yNum) || _grid->_gMark[i][j - 1][k] == AIR || _grid->_gMark[i][j][k] == SOLID;
			}
	//z
	for (int i = 0; i < _grid->xNum; i++)
		for (int j = 0; j < _grid->yNum; j++)
			for (int k = 0; k < _grid->zNum + 1; k++)
			{
				//fluid or near fluid
				faceMarkZ[0][i][j][k] = (k > 0 && _grid->_gMark[i][j][k - 1] == WATER) || (k < _grid->zNum && _grid->_gMark[i][j][k] == WATER);
				//near solid or air
				faceMarkZ[1][i][j][k] = (k <= 0 || k >= _grid->zNum) || _grid->_gMark[i][j][k - 1] == AIR || _grid->_gMark[i][j][k] == SOLID;
			}

	//extropote
	//x
	for (int i = 1; i < _grid->xNum; i++)
		for (int j = 0; j < _grid->yNum; j++)
			for (int k = 0; k < _grid->zNum; k++)
			{
				if (!faceMarkX[0][i][j][k] && faceMarkX[1][i][j][k])
				{
					int wsum = 0;
					double sum = 0.0;
					int q[][3] = { { i - 1, j, k }, { i + 1, j, k }, { i, j - 1, k }, { i, j + 1, k }, { i, j, k - 1 }, { i, j, k + 1 } };
					for (int qk = 0; qk < 6; qk++) {
						if (q[qk][0] >= 0 && q[qk][0] < _grid->xNum + 1 && q[qk][1] >= 0 && q[qk][1] < _grid->yNum  && q[qk][2] >= 0 && q[qk][2] < _grid->zNum)
						{
							if (faceMarkX[0][q[qk][0]][q[qk][1]][q[qk][2]])
							{
								wsum++;
								sum += _grid->_vx[q[qk][0]][q[qk][1]][q[qk][2]];
							}
						}
					}
					if (wsum) _grid->_vx[i][j][k] = sum / wsum;
				}
			}
	//y
	for (int i = 0; i < _grid->xNum; i++)
		for (int j = 1; j < _grid->yNum; j++)
			for (int k = 0; k < _grid->zNum; k++)
			{
				if (!faceMarkY[0][i][j][k] && faceMarkY[1][i][j][k])
				{
					int wsum = 0;
					double sum = 0.0;
					int q[][3] = { { i - 1, j, k }, { i + 1, j, k }, { i, j - 1, k }, { i, j + 1, k }, { i, j, k - 1 }, { i, j, k + 1 } };
					for (int qk = 0; qk < 6; qk++) {
						if (q[qk][0] >= 0 && q[qk][0] < _grid->xNum && q[qk][1] >= 0 && q[qk][1] < _grid->yNum + 1 && q[qk][2] >= 0 && q[qk][2] < _grid->zNum)
						{
							if (faceMarkY[0][q[qk][0]][q[qk][1]][q[qk][2]])
							{
								wsum++;
								sum += _grid->_vy[q[qk][0]][q[qk][1]][q[qk][2]];
							}
						}
					}
					if (wsum) _grid->_vy[i][j][k] = sum / wsum;
				}
			}
	//z
	for (int i = 0; i < _grid->xNum; i++)
		for (int j = 0; j < _grid->yNum; j++)
			for (int k = 1; k < _grid->zNum; k++)
			{
				if (!faceMarkZ[0][i][j][k] && faceMarkZ[1][i][j][k])
				{
					int wsum = 0;
					double sum = 0.0;
					int q[][3] = { { i - 1, j, k }, { i + 1, j, k }, { i, j - 1, k }, { i, j + 1, k }, { i, j, k - 1 }, { i, j, k + 1 } };
					for (int qk = 0; qk < 6; qk++) {
						if (q[qk][0] >= 0 && q[qk][0] < _grid->xNum && q[qk][1] >= 0 && q[qk][1] < _grid->yNum  && q[qk][2] >= 0 && q[qk][2] < _grid->zNum + 1)
						{
							if (faceMarkZ[0][q[qk][0]][q[qk][1]][q[qk][2]])
							{
								wsum++;
								sum += _grid->_vz[q[qk][0]][q[qk][1]][q[qk][2]];
							}
						}
					}
					if (wsum) _grid->_vz[i][j][k] = sum / wsum;
				}
			}

}
//------------------------------------------------------------------------------
//solve pressure
//and advect grid velocity
//------------------------------------------------------------------------------ 

void FlipFluid::solvePressureAndNewV()
{
	//get divgence
	FOR_EACH_CELL(_grid->xNum, _grid->yNum, _grid->zNum)
	{
		if (_grid->_gMark[i][j][k] == WATER)
			_grid->_div[i][j][k] = -DivOfVelocity(i, j, k);
	}

	_Solver->solve(_grid->_pressure, _grid->_div);

	//x
	FOR_EACH_CELL(_grid->xNum + 1, _grid->yNum, _grid->zNum)
	{
		if (0 < i && i < _grid->xNum)
		{
			double pf = _grid->_pressure[i][j][k];
			double pb = _grid->_pressure[i - 1][j][k];
			_grid->_vx[i][j][k] -= (pf - pb) / _grid->deltaH;

		}
	}

	//y
	FOR_EACH_CELL(_grid->xNum, _grid->yNum + 1, _grid->zNum)
	{
		if (0 < j && j < _grid->yNum)
		{
			double pf = _grid->_pressure[i][j][k];
			double pb = _grid->_pressure[i][j - 1][k];
			_grid->_vy[i][j][k] -= (pf - pb) / _grid->deltaH;
		}
	}
	//z
	FOR_EACH_CELL(_grid->xNum, _grid->yNum, _grid->zNum + 1)
	{
		if (0 < k && k < _grid->zNum)
		{
			double pf = _grid->_pressure[i][j][k];
			double pb = _grid->_pressure[i][j][k - 1];
			_grid->_vz[i][j][k] -= (pf - pb) / _grid->deltaH;
		}
	}
}

//------------------------------------------------------------------------------
//simulate function.
//------------------------------------------------------------------------------
void FlipFluid::simulate()
{
	//match particle to corresponding cell
	Timer::tic();
	ParCellIdentification();
	Timer::toc("ParCellIdentification");
	//grid marking
	Timer::tic();
	FluidSurfaceMark();
	Timer::toc("FluidSurfaceMark");
	//calculate density
#ifdef COMPUTE_DENSITY
	Timer::tic();
	computeDensity();
	Timer::toc("computeDensity");
#endif
	//calculate gravity force
	Timer::tic();
	addGravityF();
	Timer::toc("addGravity");
	//store particle velocties to grid
	Timer::tic();
	P2G();
	Timer::toc("P2G");
	//temp store current grid velocities
	mycopy(_grid->_vx, _grid->_vx_save, _grid->xNum + 1, _grid->yNum, _grid->zNum);
	mycopy(_grid->_vy, _grid->_vy_save, _grid->xNum, _grid->yNum + 1, _grid->zNum);
	mycopy(_grid->_vz, _grid->_vz_save, _grid->xNum, _grid->yNum, _grid->zNum + 1);
	//solve pressure to grid //calculate new grid velocities
	setBoundaryVe();
	Timer::tic();
	solvePressureAndNewV();
	//set boundary velocities to zero
	setBoundaryVe();
	Timer::toc("solvePressureAndNewV");
	//re evaluate grid velocities in regions that mightnot be defined
	velocityExtrapolation();
	//flip
	//calculate new flip velocity
	//subtract current grid from saved grid velocities and store 'em
	Timer::tic();
	FOR_EACH_CELL(_grid->xNum + 1, _grid->yNum, _grid->zNum)
	{
		_grid->_vx_save[i][j][k] = _grid->_vx[i][j][k] - _grid->_vx_save[i][j][k];
	}
	FOR_EACH_CELL(_grid->xNum, _grid->yNum + 1, _grid->zNum)
	{
		_grid->_vy_save[i][j][k] = _grid->_vy[i][j][k] - _grid->_vy_save[i][j][k];
	}
	FOR_EACH_CELL(_grid->xNum, _grid->yNum, _grid->zNum + 1)
	{
		_grid->_vz_save[i][j][k] = _grid->_vz[i][j][k] - _grid->_vz_save[i][j][k];
	}
	for (int i = 0; i < _particle.size(); i++)
	{
		_particle[i]._tempv.X = _particle[i]._v.X;
		_particle[i]._tempv.Y = _particle[i]._v.Y;
		_particle[i]._tempv.Z = _particle[i]._v.Z;
	}

	G2P(_grid->_vx_save, _grid->_vy_save, _grid->_vz_save);
	//save temp as flip v
	for (int i = 0; i < _particle.size(); i++)
	{
		_particle[i]._tempv.X = _particle[i]._v.X + _particle[i]._tempv.X;
		_particle[i]._tempv.Y = _particle[i]._v.Y + _particle[i]._tempv.Y;
		_particle[i]._tempv.Z = _particle[i]._v.Z + _particle[i]._tempv.Z;
	}
	//pic
	//calculate new pic velocity
	G2P(_grid->_vx, _grid->_vy, _grid->_vz);
	//interpolate pic flip velocities
	for (int i = 0; i < _particle.size(); i++)
	{
		_particle[i]._v.X = _particle[i]._v.X * VCOEFFICIENT + (1 - VCOEFFICIENT) * _particle[i]._tempv.X;
		_particle[i]._v.Y = _particle[i]._v.Y * VCOEFFICIENT + (1 - VCOEFFICIENT) * _particle[i]._tempv.Y;
		_particle[i]._v.Z = _particle[i]._v.Z * VCOEFFICIENT + (1 - VCOEFFICIENT) * _particle[i]._tempv.Z;
	}
	Timer::toc("calculate new velocity");
	//advect particles
	Timer::tic();

#ifdef USING_GV2ADVECT
	for (int i = 0; i < _particle.size(); i++)
	{
		double v[3];
		G2P_particle(&_particle[i], v, _grid->_vx, _grid->_vy, _grid->_vz);
		_particle[i]._p.X += v[0] * dt;
		_particle[i]._p.Y += v[1] * dt;
		_particle[i]._p.Z += v[2] * dt;
	}
#else
	for (int i = 0; i < _particle.size(); i++)
	{

		_particle[i]._p.X += _particle[i]._v.X * dt;
		_particle[i]._p.Y += _particle[i]._v.Y * dt;
		_particle[i]._p.Z += _particle[i]._v.Z * dt;

	}
#endif

	Timer::toc("Adevect");
	//debug << "====================================" << endl;
	//reposition particles
	reposition();
	//correct
#ifdef CORRECT
	Timer::tic();
	correct();
	Timer::toc("correct");
#endif
	//detect collisions amongst particles
	for (int i = 0; i < _particle.size(); i++)
	{
		//_particle[i]._p.X = fmax(0.0, fmin(_particle[i]._p.X, 1.0));
		//_particle[i]._p.Y = fmax(0.0, fmin(_particle[i]._p.Y, 1.0));
		//_particle[i]._p.Z = fmax(0.0, fmin(_particle[i]._p.Z, 1.0));
		if (_particle[i]._p.X < 0)
		{
			_particle[i]._p.X = 0;
			_particle[i]._v.X = -_particle[i]._v.X;
		}
		if (_particle[i]._p.X > 1.0)
		{
			_particle[i]._p.X = 1.0;
			_particle[i]._v.X = -_particle[i]._v.X;
		}
		if (_particle[i]._p.Y < 0)
		{
			_particle[i]._p.Y = 0;
			_particle[i]._v.Y = -_particle[i]._v.Y;
		}
		if (_particle[i]._p.Y > 1.0)
		{
			_particle[i]._p.Y = 1.0;
			_particle[i]._v.Y = -_particle[i]._v.Y;
		}
		if (_particle[i]._p.Z < 0)
		{
			_particle[i]._p.Z = 0;
			_particle[i]._v.Z = -_particle[i]._v.Z;
		}
		if (_particle[i]._p.Z > 1.0)
		{
			_particle[i]._p.Z = 1.0;
			_particle[i]._v.Z = -_particle[i]._v.Z;
		}
	}
}

bool FlipFluid::Init()
{
	if (!App::Init())
		return false;

	_grid = new grid;
	_Solver = new Solver(_grid);
	_sorter = new sorter(_grid->xNum, _grid->yNum, _grid->zNum);
	Timer::start();
	//calculate max density
	double dH = _grid->deltaH;
	double dH_4 = dH * 0.25;
	//initiate particle position and velocity
	FOR_EACH_CELL(10, 10, 10)
	{
		_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4, k*dH + dH_4, 0, 0, 0));
		_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4, 0, 0, 0));
		_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
		_particle.push_back(Particle(i*dH + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));


		_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4, k*dH + dH_4, 0, 0, 0));
		_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4, 0, 0, 0));
		_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
		_particle.push_back(Particle(i*dH + dH_4 + dH_4 + dH_4, j*dH + dH_4 + dH_4 + dH_4, k*dH + dH_4 + dH_4 + dH_4, 0, 0, 0));
	}
	ParCellIdentification();
	MaxDensity = 1.0;
	computeDensity();
	MaxDensity = 0.0;
	for (int i = 0; i<_particle.size(); i++)
	{
		Particle &p = _particle[i];
		MaxDensity = fmax(MaxDensity, p.density);
	}
	_particle.clear();

	reset();

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	GLfloat lightPos[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat ambientLight[] = { 0.0, 0.7, 0.8, 1.0 };
	//glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, 1);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientLight);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, ambientLight);
	//glLightfv(GL_LIGHT0, GL_SPECULAR, specular);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);
	//glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glShadeModel(GL_SMOOTH);

	return true;
}

void FlipFluid::onResize(GLFWwindow* window, int w, int h)
{
	App::onResize(window, w, h);
}

void FlipFluid::UpdateScene()
{
	float ratio;
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	ratio = width / (float)height;
	glViewport(0, 0, width, height);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0f, ratio, 0.1f, 10.0f);
	//glOrtho(0, ratio, 0.f, 1.f, 1.f, -1.f);

	gluLookAt(0.5, -1.5, 1.5,  /* eye is at (0,0,5) */
		0.5, 0.5, 0.5,      /* center is at (0,0,0) */
		0.0, 1.0, 0.);      /* up is in postivie Y direction */

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(rX, 1, 0, 0);
	glRotatef(rY, 0, 0, 1);
	//Timer::tic();
	simulate();
	//Timer::toc("simulation");
#ifdef FACES
	_sorter->sort(_particle);
	mesher::generateMesh(_grid->_gMark, _sorter, _particle, DENSITY, GRIDN + 1, vertices, normals, faces);
#endif
}
void FlipFluid::Rendering()
{
	static GLfloat mat_shininess[] = { 80.0 };

	// static GLfloat mat_diffuse[] = { 0.5, 0.2, 0.3, 1.0 };
	static GLfloat mat_diffuse[] = { 0.23, 0.23, 0.86, 1.0 };
	static GLfloat mat_specular[] = { 0.23, 0.23, 1.0, 1.0 };

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

	DrawCube();
	//DrawGrid();
	//DrawVOfLocalP();
#ifdef	DRAW_PARTICLES
	DrawParticle();
#endif
#ifdef FACES
	DrawFaces();
#endif
#ifdef OB
	if (drawOB)
	{
		DrawSolid();
	}

#endif
	Timer::AddFrame();
	string ss = "FPS: " + to_string(Timer::fps());
	glfwSetWindowTitle(window, ss.c_str());
}
void FlipFluid::onMouseWheel(GLFWwindow* window, double x, double y)
{

}
void FlipFluid::onMouseMove(GLFWwindow* window, double x, double y)
{
	if (isMouseButton)
	{
		rY += (x - oldX) / 5.0f;
		rX += (y - oldY) / 5.0f;
	}

	oldX = x;
	oldY = y;
}
void FlipFluid::onMouseButton(GLFWwindow* window, int button, int action, int mods)
{
	double xd, yd;

	if ((button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_PRESS))
	{
		glfwGetCursorPos(window, &xd, &yd);
		isMouseButton = !isMouseButton;
		oldX = xd;
		oldY = yd;
	}
	else if ((button == GLFW_MOUSE_BUTTON_1) && (action == GLFW_RELEASE))
	{
		isMouseButton = !isMouseButton;
	}
}
void FlipFluid::onKey(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	App::onKey(window, key, scancode, action, mods);
	if ((key == GLFW_KEY_Q) && (action == GLFW_PRESS))
	{
		exit(0);
	}
	if ((key == GLFW_KEY_T) && (action == GLFW_PRESS))
	{
		Timer::printTime();
	}
	if ((key == GLFW_KEY_1) && (action == GLFW_PRESS))
	{
		reset();
	}
	if ((key == GLFW_KEY_2) && (action == GLFW_PRESS))
	{
		reset2();
	}
	if ((key == GLFW_KEY_P))
	{
		pourWater();
	}
	if ((key == GLFW_KEY_A))
	{
		drawOB = !drawOB;
	}
}