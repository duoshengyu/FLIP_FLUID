#include "solver.h"

Solver::Solver(grid *grid)
{
	_gridX = grid->xNum;
	_gridY = grid->yNum;
	_gridZ = grid->zNum;
	_gridReso = _gridX * _gridY * _gridZ;
	_deltaX2 = grid->deltaH * grid->deltaH;
	_mark = grid->_gMark;
	pre = memAlloc3D<double>(_gridX, _gridY, _gridZ);

	z = memAlloc3D<double>(_gridX, _gridY, _gridZ);
	r = memAlloc3D<double>(_gridX, _gridY, _gridZ);
	s = memAlloc3D<double>(_gridX, _gridY, _gridZ);
	q = memAlloc3D<double>(_gridX, _gridY, _gridZ);
}
void Solver::MICprecon()
{
	double tau = 0.97;

	FOR_EACH_CELL(_gridX, _gridY, _gridZ)
	{
		if (_mark[i][j][k] != WATER)
			continue;
		double diag = Adiag(i, j, k);
		double aplusi_1jk = Aplus(i - 1, j, k, i, j, k);
		double aplusij_1k = Aplus(i, j - 1, k, i, j, k);
		double aplusijk_1 = Aplus(i, j, k - 1, i, j, k);
		double preconi_1jk = getDataInVector(pre, i - 1, j, k);
		double preconij_1k = getDataInVector(pre, i, j - 1, k);
		double preconijk_1 = getDataInVector(pre, i, j, k - 1);
#if 0
		debug << diag << "  " << aplusi_1jk << "  " << aplusij_1k << "  " << aplusijk_1 << "  " << preconi_1jk << "  " << preconij_1k << "  " << preconijk_1 << endl;
#endif
		double e = diag
			- (aplusi_1jk * preconi_1jk) * (aplusi_1jk * preconi_1jk)
			- (aplusij_1k * preconij_1k) * (aplusij_1k * preconij_1k)
			- (aplusijk_1 * preconijk_1) * (aplusijk_1 * preconijk_1);

#ifdef MIC
		double tauEntry1 = aplusi_1jk * (aplusi_1jk + aplusi_1jk) * preconi_1jk * preconi_1jk;
		double tauEntry2 = aplusij_1k * (aplusij_1k + aplusij_1k) * preconij_1k * preconij_1k;
		double tauEntry3 = aplusijk_1 * (aplusijk_1 + aplusijk_1) * preconijk_1 * preconijk_1;

		e -= tau * (tauEntry1 + tauEntry2 + tauEntry3);
#endif
		pre[i][j][k] = 1.0 / sqrt(e + 10e-30);
	}
}
void Solver::applyPrecon(double ***z, double ***r)
{
	//First solve Lq = r
	FOR_EACH_CELL(_gridX, _gridY, _gridZ)
	{
		if (_mark[i][j][k] != WATER)
			continue;

		double left = Aplus(i - 1, j, k, i, j, k) * getDataInVector(pre, i - 1, j, k) * getDataInVector(q, i - 1, j, k);
		double back = Aplus(i, j - 1, k, i, j, k) * getDataInVector(pre, i, j - 1, k) * getDataInVector(q, i, j - 1, k);
		double bottom = Aplus(i, j, k - 1, i, j, k) * getDataInVector(pre, i, j, k - 1) * getDataInVector(q, i, j, k - 1);

		double t = r[i][j][k] - left - back - bottom;
		q[i][j][k] = t * getDataInVector(pre, i, j, k);
	}
	//next solve L^Tz = q
	FOR_EACH_CELL_FLIP(_gridX, _gridY, _gridZ)
	{
		if (_mark[i][j][k] != WATER)
			continue;

		double right = Aplus(i, j, k, i, j, k) * getDataInVector(pre, i, j, k) * getDataInVector(z, i + 1, j, k);
		double front = Aplus(i, j, k, i, j, k) * getDataInVector(pre, i, j, k) * getDataInVector(z, i, j + 1, k);
		double top = Aplus(i, j, k, i, j, k) * getDataInVector(pre, i, j, k) * getDataInVector(z, i, j, k + 1);

		double t = q[i][j][k] - right - front - top;
		z[i][j][k] = t * getDataInVector(pre, i, j, k);
	}
}
void Solver::CG(double ***x, double ***b)
{
#if 0
	FOR_EACH_CELL(_gridX, _gridY, _gridZ)
	{
		debug << x[i][j][k] << endl;
	}
#endif
	Ax(z, x);//z = Ax
#if 0
	FOR_EACH_CELL(_gridX, _gridY, _gridZ)
	{
		debug << z[i][j][k] << endl;
	}
#endif

	b_Ax(r, b, z);//r = b-Ax;
#if 0
	FOR_EACH_CELL(_gridX, _gridY, _gridZ)
	{
		debug << b[i][j][k] << endl;
	}
#endif
	double error2 = dot_product(r, r);

	applyPrecon(z, r);//z = apply(r)
#if 0
	FOR_EACH_CELL(_gridX, _gridY, _gridZ)
	{
		debug << z[i][j][k] << endl;
	}
#endif
	copy(z, s);//s = z

	double a = dot_product(z, r);
	double eps = 1.0e-6 * _gridReso;

	for (int i = 0; i < _gridReso; i++)
	{
		Ax(z, s);// z = Ax
		double alpha = a / dot_product(z, s);// alpha = a/(z`s)
		xplusAlphaS(x, s, alpha);//x += alpha * s;
		xplusAlphaS(r, z, -alpha);//r-= alpha * z;

		double error = dot_product(r, r);
		error2 = fmin(error, error2);
#if 0
		debug << "error:" << error << endl;
#endif
		if (error2 <= eps)
		{
#if 0
			debug << i << endl;
#endif
			break;
		}

		applyPrecon(z, r);
		double anew = dot_product(z, r);
		double beta = anew / a;
		splusBetaSZ(s, z, beta);//s = z + beta*s
		a = anew;
	}

}