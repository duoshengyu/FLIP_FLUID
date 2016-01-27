/*
*  mesher.cpp
*  flip3D
*/
#ifndef MESHER_H
#define MESHER_H
#include "sorter.h"
#include <vector>

namespace mesher {
	void generateMesh(char ***A, sorter *sort, std::vector<Particle> &particles, double density, int mg,
		std::vector<double> &vertices, std::vector<double> &normals, std::vector<int> &faces);
}
#include "util.h"
#include "implicit.h"
#include "CIsoSurface.h"
#include <math.h>
using namespace std;

static void smoothMesh(sorter *sort, vector<double> &vertices, vector<double> &normals, vector<vector<int> > &connections,
	int numit, int do_normal = 1, int do_vertices = 1)
{
	vector<double> copy_vetices;
	vector<double> copy_normals;

	copy_vetices.resize(vertices.size());
	copy_normals.resize(normals.size());

	for (int k = 0; k < numit; k++) {
		for (int vi = 0; vi < connections.size(); vi++) {
			double wsum = 0.0;
			double v[3] = { 0.0, 0.0, 0.0 };
			double vn[3] = { 0.0, 0.0, 0.0 };

			// Add Itself
			{
				bool isOK = true;
				for (int n = 0; n < 3; n++) {
					if (!isnormal(vertices[3 * vi + n])) isOK = false;
					if (!isnormal(normals[3 * vi + n])) isOK = false;
				}

				if (isOK) {
					double w = 1.0;
					for (int n = 0; n < 3; n++) {
						v[n] += w*vertices[3 * vi + n];
						vn[n] += w*normals[3 * vi + n];
					}
					wsum += w;
				}
			}

			// Add Neighbors
			for (int i = 0; i < connections[vi].size(); i++) {
				int ni = connections[vi][i];
				double w = 1.0;

				bool isOK = true;
				for (int n = 0; n < 3; n++) {
					if (!isnormal(vertices[3 * ni + n])) isOK = false;
					if (!isnormal(normals[3 * ni + n])) isOK = false;
				}

				if (isOK) {
					for (int n = 0; n < 3; n++) {
						v[n] += w*vertices[3 * ni + n];
						vn[n] += w*normals[3 * ni + n];
					}
					wsum += w;
				}
			}

			double len = hypotf(vn[0], hypotf(vn[1], vn[2]));
			for (int n = 0; n < 3; n++) {
				if (wsum) copy_vetices[3 * vi + n] = v[n] / (double)wsum;
				if (len) copy_normals[3 * vi + n] = vn[n] / len;
			}
		}


		for (int vi = 0; vi < connections.size(); vi++) {
			double a = 1.0;
			if (isnan(a)) {
				a = 1.0;
				for (int n = 0; n < 3; n++) vertices[3 * vi + n] = normals[3 * vi + n] = 0.0;
			}

			for (int n = 0; n < 3; n++) {
				if (do_vertices) vertices[3 * vi + n] = copy_vetices[3 * vi + n];
				if (do_normal) normals[3 * vi + n] = copy_normals[3 * vi + n];
			}
		}
	}
}

static vector<vector<int> > buildConnection(std::vector<double> &vertices, std::vector<double> &normals, std::vector<int> &faces)
{
	vector<vector<int> > connections;
	connections.resize(vertices.size() / 3);
	for (int i = 0; i < faces.size(); i += 3) {
		connections[faces[i]].push_back(faces[i + 1]);
		connections[faces[i]].push_back(faces[i + 2]);
		connections[faces[i + 1]].push_back(faces[i]);
		connections[faces[i + 1]].push_back(faces[i + 2]);
		connections[faces[i + 2]].push_back(faces[i]);
		connections[faces[i + 2]].push_back(faces[i + 1]);
	}
	return connections;
}
/*
static void alignNormalToWalls(sorter *sort, vector<double> &vertices, vector<double> &normals)
{
int gn = sort->getCellSizeX();
for (int i = 0; i<vertices.size() / 3; i++) {
double p[3] = { vertices[3 * i + 0], vertices[3 * i + 1], vertices[3 * i + 2] };
double n[3] = { normals[3 * i + 0], normals[3 * i + 1], normals[3 * i + 2] };

vector<Particle *> neighbors = sort->getNeigboringParticles_cell(fmax(0, fmin(gn - 1, gn*p[0])),
fmax(0, fmin(gn - 1, gn*p[1])),
fmax(0, fmin(gn - 1, gn*p[2])), 1, 1, 1);
double min_len = 99999.0;
Particle *min_p = NULL;
//for (int m = 0; m<neighbors.size(); m++)
//{
//	Particle *np = neighbors[m];
//	if (np->type == WALL)
//	{
//		double len = length(np->p, p);
//		if (len < min_len) {
//			min_p = np;
//			min_len = len;
//		}
//	}
//}

if (min_len < 2.0 / gn && n[0] * min_p->n[0] + n[1] * min_p->n[1] + n[2] * min_p->n[2] < 0) {
n[0] = -min_p->n[0];
n[1] = -min_p->n[1];
n[2] = -min_p->n[2];

normals[3 * i + 0] = n[0];
normals[3 * i + 1] = n[1];
normals[3 * i + 2] = n[2];
}
}
}
*/
void mesher::generateMesh(char ***A, sorter *sort, std::vector<Particle> &particles, double density, int mg,
	std::vector<double> &vertices, std::vector<double> &normals, std::vector<int> &faces) {
	static double *dense_grid = new double[mg*mg*mg];

	// Create Density Field
	FOR_EACH_CELL(mg,mg,mg) {
		double h = 1.0 / (double)(mg - 1);
		double x = i*h;
		double y = j*h;
		double z = k*h;
		double p[3] = { x, y, z };
		double value = implicit::implicit_func(sort, p, density);
		if (i == 0 || i == mg - 1 || j == 0 || j == mg - 1 || k == 0 || k == mg - 1) {
			value = fmax(value, 0.01);
		}
		dense_grid[(mg*mg)*k + mg*j + i] = -value;
#if 0
		debug << i << "  " << j << "  " << k << "  "<< value <<endl;
#endif
	}

	// Extract mesh
	vertices.clear();
	normals.clear();
	faces.clear();
	CIsoSurface<double> surface;
	surface.GenerateSurface(dense_grid, 0, mg - 1, mg - 1, mg - 1, 1.0 / (mg - 1), 1.0 / (mg - 1), 1.0 / (mg - 1));

	for (int i = 0; i<surface.m_nVertices; i++) 
	{
		for (int n = 0; n<3; n++) 
			vertices.push_back(surface.m_ppt3dVertices[i][n]);
	}


	for (unsigned int i = 0; i < surface.m_nTriangles; i++) 
	{
		unsigned int id0, id1, id2;
		id0 = surface.m_piTriangleIndices[i * 3];
		id1 = surface.m_piTriangleIndices[i * 3 + 1];
		id2 = surface.m_piTriangleIndices[i * 3 + 2];

		double v0[3] = { vertices[3 * id0 + 0], vertices[3 * id0 + 1], vertices[3 * id0 + 2] };
		double v1[3] = { vertices[3 * id1 + 0], vertices[3 * id1 + 1], vertices[3 * id1 + 2] };
		double v2[3] = { vertices[3 * id2 + 0], vertices[3 * id2 + 1], vertices[3 * id2 + 2] };

		double maxlen = 0.1;
		double maxlen2 = maxlen*maxlen;
		if (length2(v0, v1) < maxlen2 && length2(v1, v2) < maxlen2 && length2(v0, v2) < maxlen2) {
			faces.push_back(id0);
			faces.push_back(id1);
			faces.push_back(id2);
		}
	}

	for (int i = 0; i<surface.m_nNormals; i++) {
		for (int n = 0; n<3; n++) normals.push_back(surface.m_pvec3dNormals[i][n]);
	}

	// Smooth Vertices
	vector<vector<int> > connections = buildConnection(vertices, normals, faces);
	smoothMesh(sort, vertices, normals, connections, 2, 1, 1);

	// Align Normal Directions
	// alignNormalToWalls( sort, vertices, normals );
}

#endif