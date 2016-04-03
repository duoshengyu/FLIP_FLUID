#ifndef SOLVER_H
#define SOLVER_H
//------------------------------------------------------------------------------
//PCG solver. Use to solve divergence free pressure.
//Input is divergence of grid, output is pressure.
//------------------------------------------------------------------------------
#include "util.h"

class Solver
{
public:
	Solver(grid *grid);
	~Solver();
	void CG(double ***x, double ***b);
	void MICprecon();
	void applyPrecon(double ***z, double ***r);

	//copy a to b
	void copy(double ***a, double ***b)
	{
		FOR_EACH_CELL(_gridX, _gridY, _gridZ)
		{
			b[i][j][k] = a[i][j][k];
		}
	}

	//dot product between two vector (grid based), have not assert check.
	double dot_product(double ***a, double ***b)
	{
		double sum = 0.0;
		FOR_EACH_CELL(_gridX, _gridY, _gridZ)
		{
			sum += (a[i][j][k] * b[i][j][k]);
		}
		return sum;
	}

	//Calculate A * x.
	void Ax(double ***result, double ***x)
	{
		FOR_EACH_CELL(_gridX, _gridY, _gridZ)
		{
			if (_mark[i][j][k] == WATER)
			{
				result[i][j][k] = (6.0 * x[i][j][k] - x_ref(x, i - 1, j, k, i, j, k) 
					- x_ref(x, i + 1, j, k, i, j, k) - x_ref(x, i, j - 1, k, i, j, k)
					- x_ref(x, i, j + 1, k, i, j, k) - x_ref(x, i, j, k - 1, i, j, k) 
					- x_ref(x, i, j, k + 1, i, j, k)) / _deltaX2;
			}
			else
				result[i][j][k] = 0.0;
		}
	}

	//result = b - Ax.
	void b_Ax(double ***result, double ***b, double ***Ax)
	{
		FOR_EACH_CELL(_gridX, _gridY, _gridZ)
		{
			if (_mark[i][j][k] == WATER)
			{
				result[i][j][k] = b[i][j][k] - Ax[i][j][k];
			}
			else
				result[i][j][k] = 0.0;
		}
	}

	//x = s * alpha.
	void xplusAlphaS(double ***x, double ***s, double alpha)
	{
		FOR_EACH_CELL(_gridX, _gridY, _gridZ)
		{
			if (_mark[i][j][k] == WATER)
			{
				x[i][j][k] += s[i][j][k] * alpha;
			}
			else
				x[i][j][k] = 0.0;
		}
	}
	//s = z + beta*s
	void splusBetaSZ(double ***x, double ***s, double alpha)
	{
		FOR_EACH_CELL(_gridX, _gridY, _gridZ)
		{
			if (_mark[i][j][k] == WATER)
			{
				x[i][j][k] = s[i][j][k] + x[i][j][k] * alpha;
			}
			else
				x[i][j][k] = 0.0;
		}
	}
	//get corresponding divergence between grid (i j k) and gird(I J K).
	double x_ref(double ***x, int i, int j, int k, int I, int J, int K)
	{
		i = min(max(0, i), _gridX - 1);
		j = min(max(0, j), _gridY - 1);
		k = min(max(0, k), _gridZ - 1);
		if (_mark[i][j][k] == WATER)
			return x[i][j][k];
		if (_mark[i][j][k] == SOLID)
			return x[I][J][K];
		return 0.0;
	}
	double Aplus(int i, int j, int k, int I, int J, int K)
	{
		if (0 > i || i > _gridX - 1 || 0 > j || j > _gridY - 1 || 0 > k || k > _gridZ - 1 || _mark[i][j][k] != WATER) return 0.0;
		if (0 > I || I > _gridX - 1 || 0 > J || J > _gridY - 1 || 0 > K || K > _gridZ - 1 || _mark[I][J][K] != WATER) return 0.0;
		return -1.0;
	}
	double Adiag(int i, int j, int k)
	{
		double diag = 6.0;
		if (_mark[i][j][k] != WATER) return diag;
		int index[6][3] = { { i - 1, j, k }, { i + 1, j, k }, { i, j - 1, k }, { i, j + 1, k }, { i, j, k - 1 }, { i, j, k + 1 } };
		for (int l = 0; l < 6; l++)
		{
			int I = index[l][0];
			int J = index[l][1];
			int K = index[l][2];
			if (I < 0 || I > _gridX - 1 || J < 0 || J > _gridY - 1 || K < 0 || K > _gridZ - 1 || _mark[I][J][K] == SOLID) diag -= 1.0;
		}
		return diag;
	}
	double getDataInVector(double ***data, int i, int j, int k)
	{
		if (i < 0 || i > _gridX - 1 || j < 0 || j > _gridY - 1 || k < 0 || k > _gridZ - 1 || _mark[i][j][k] != WATER) return 0.0;
		return data[i][j][k];
	}
	void solve(double ***x, double ***b)
	{
		//build preconditioner
		MICprecon();
		//cg solve
		CG(x, b);
	}
private:
	int _gridX;
	int _gridY;
	int _gridZ;
	int _gridReso;
	double _deltaX2;
	char ***_mark;
	double ***pre;



	double ***z;
	double ***r;
	double ***s;
	double ***q;

};

#endif