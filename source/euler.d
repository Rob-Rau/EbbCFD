/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.euler;

import std.math;

import numd.linearalgebra.matrix;

import ebb.mesh;

alias Vec = Vector!4;
alias Mat = Matrix!(4, 4);

T noop(T)(T x)
{
	return x;
}

Mat L(int kx, int ky)(Vec q)
{
	immutable double rho = q[0];
	immutable double u = q[1]/rho;
	immutable double v = q[2]/rho;
	immutable double E = q[$-1];
	immutable double p = ((gamma - 1)*(E - 0.5*rho*(u^^2.0 + v^^2.0)));
	immutable double a = sqrt(gamma*(p/rho));
	return L!(kx, ky)(u, v, a);
}

Mat L(int kx, int ky)(double u, double v, double a)
{	
	immutable double phi2 = 0.5*(gamma - 1)*(u^^2.0 + v^^2.0);
	immutable double beta = 1.0/(2.0*a^^2.0);
	immutable double kxt = kx/sqrt(kx^^2.0 + ky^^2.0);
	immutable double kyt = ky/sqrt(kx^^2.0 + ky^^2.0);
	immutable double theta = kxt*v + kyt*v;

	Mat mat = Mat([1 - phi2/a^^2.0, (gamma - 1)*u/a^^2.0, (gamma - 1)*v/a^^2.0, -(gamma - 1)/a^^2.0, 
					-(kyt*u - kxt*v), kyt, -kxt, 0,
					beta*(phi2 - a*theta), beta*(kxt*a - (gamma - 1)*u), beta*(kyt*a - (gamma - 1)*v), beta*(gamma - 1),
					beta*(phi2 + a*theta), -beta*(kxt*a + (gamma - 1)*u), -beta*(kyt*a + (gamma - 1)*v), beta*(gamma - 1)]);
	return mat;
}

Mat R(int kx, int ky)(Vec q)
{
	immutable double rho = q[0];
	immutable double u = q[1]/rho;
	immutable double v = q[2]/rho;
	immutable double E = q[$-1];
	immutable double p = ((gamma - 1)*(E - 0.5*rho*(u^^2.0 + v^^2.0)));
	immutable double a = sqrt(gamma*(p/rho));
	
	return R!(kx, ky)(u, v, a);
}

Mat R(int kx, int ky)(double u, double v, double a)
{
	immutable double phi2 = 0.5*(gamma - 1)*(u^^2.0 + v^^2.0);
	immutable double kxt = kx/sqrt(kx^^2.0 + ky^^2.0);
	immutable double kyt = ky/sqrt(kx^^2.0 + ky^^2.0);
	immutable double theta = kxt*v + kyt*v;

	Mat mat = Mat([1, 0, 1, 1,
					u, kyt, u + kxt*a, u - kxt*a,
					v, -kxt, v + kyt*a, v - kyt*a,
					phi2/(gamma - 1), kyt*u - kxt*v, (phi2 + a^^2.0)/(gamma - 1) + a*theta, (phi2 + a^^2.0)/(gamma - 1) - a*theta]);
	return mat;
}

Mat Lam(int kx, int ky, alias func = noop)(Vec q)
{
	immutable double rho = q[0];
	immutable double u = q[1]/rho;
	immutable double v = q[2]/rho;
	immutable double E = q[$-1];
	immutable double p = ((gamma - 1)*(E - 0.5*rho*(u^^2.0 + v^^2.0)));
	immutable double a = sqrt(gamma*(p/rho));
	
	return Lam!(kx, ky, func)(u, v, a);
}
	
Mat Lam(int kx, int ky, alias func = noop)(double u, double v, double a)
{
	immutable double U = kx*u + ky*v;
	Mat mat = Mat([func(U), 0, 0, 0, 0, func(U), 0, 0, 0, 0, func(U + a), 0, 0, 0, 0, func(U - a)]);
	return mat;
}

immutable double gamma = 1.4;

template getVelocity(int dim)
	if((dim == 0) || (dim == 1))
{
	double getVelocity(ref Cell cell)
	{
		double vel;
		static if(dim == 0)
		{
			vel = cell.q[1]/cell.q[0];
		}
		else static if(dim == 1)
		{
			vel = cell.q[2]/cell.q[0];
		}

		return vel;
	}

	double[][] getVelocity(ref Mesh mesh)
	{
		double[][] vel = new double[][](mesh.M, mesh.N);
		for(int i = 0; i < mesh.N; i++)
		{
			for(int j = 0; j < mesh.M; j++)
			{
				static if(dim == 0)
				{
					vel[j][i] = mesh[i,j].q[1]/mesh[i,j].q[0];
				}
				else static if(dim == 1)
				{
					vel[j][i] = mesh[i,j].q[2]/mesh[i,j].q[0];
				}
			}
		}
		return vel;
	}
}

double[][] getDensity(ref Mesh mesh)
{
	double[][] rho = new double[][](mesh.M, mesh.N);
	for(int i = 0; i < mesh.N; i++)
	{
		for(int j = 0; j < mesh.M; j++)
		{
			rho[j][i] = mesh[i,j].q[0];
		}
	}
	return rho;
}

@nogc double getPressure(ref Vector!4 q)
{
	return (gamma - 1)*(q[3] - 0.5*q[0]*((q[1]/q[0])^^2.0 + (q[2]/q[0])^^2.0));
}

double getPressure(ref Cell cell)
{
	return (gamma - 1)*(cell.q[3] - 0.5*cell.q[0]*((cell.q[1]/cell.q[0])^^2.0 + (cell.q[2]/cell.q[0])^^2.0));
}

double[][] getPressure(ref Mesh mesh)
{
	double[][] p = new double[][](mesh.M, mesh.N);
	for(int i = 0; i < mesh.N; i++)
	{
		for(int j = 0; j < mesh.M; j++)
		{
			if(mesh[i,j].cellType == CellType.Normal)
			{
				p[j][i] = (gamma - 1)*(mesh[i,j].q[3] - 0.5*mesh[i,j].q[0]*((mesh[i,j].q[1]/mesh[i,j].q[0])^^2.0 + (mesh[i,j].q[2]/mesh[i,j].q[0])^^2.0));
			}
			else
			{
				p[j][i] = double.nan;
			}
		}
	}
	return p;
}

double[][] getMach(ref Mesh mesh)
{
	double[][] M = new double[][](mesh.M, mesh.N);
	for(int i = 0; i < mesh.N; i++)
	{
		for(int j = 0; j < mesh.M; j++)
		{
			if(mesh[i,j].cellType == CellType.Normal)
			{
				double p = (gamma - 1)*(mesh[i,j].q[3] - 0.5*mesh[i,j].q[0]*((mesh[i,j].q[1]/mesh[i,j].q[0])^^2.0 + (mesh[i,j].q[2]/mesh[i,j].q[0])^^2.0));
				double rho = mesh[i,j].q[0];
				double a = sqrt(gamma*p/rho);
				double u = mesh[i,j].q[1]/rho;
				double v = mesh[i,j].q[2]/rho;
				double speed = sqrt(u^^2.0 + v^^2.0);
				M[j][i] = speed/a;
			}
			else
			{
				M[j][i] = double.nan;
			}
		}
	}
	return M;
}

double[][] getSpeed(ref Mesh mesh)
{
	double[][] speed = new double[][](mesh.M, mesh.N);
	for(int i = 0; i < mesh.N; i++)
	{
		for(int j = 0; j < mesh.M; j++)
		{
			if(mesh[i,j].cellType == CellType.Normal)
			{
				double rho = mesh[i,j].q[0];
				double u = mesh[i,j].q[1]/rho;
				double v = mesh[i,j].q[2]/rho;
				speed[j][i] = sqrt(u^^2.0 + v^^2.0);
			}
			else
			{
				speed[j][i] = double.nan;
			}
		}
	}
	return speed;
}

double[] getSpeed2(ref Mesh mesh)
{
	double[] speed = new double[mesh.M*mesh.N];
	for(int i = 0; i < mesh.N*mesh.M; i++)
	{
		if(mesh.cells[i].cellType == CellType.Normal)
		{
			double rho = mesh.cells[i].q[0];
			double u = mesh.cells[i].q[1]/rho;
			double v = mesh.cells[i].q[2]/rho;
			speed[i] = sqrt(u^^2.0 + v^^2.0);
		}
		else
		{
			speed[i] = 0;//double.nan;
		}
	}
	return speed;
}

double getSoundSpeed(ref Cell cell)
{
	double p = getPressure(cell);
	double rho = cell.q[0];
	return sqrt(gamma*p/rho);
}

double[] getSoundSpeed(ref Mesh mesh)
{
	double[] a = new double[mesh.M*mesh.N];
	for(int i = 0; i < mesh.N*mesh.M; i++)
	{
		if(mesh.cells[i].cellType == CellType.Normal)
		{
			double p = (gamma - 1)*(mesh.cells[i].q[3] - 0.5*mesh.cells[i].q[0]*((mesh.cells[i].q[1]/mesh.cells[i].q[0])^^2.0 + (mesh.cells[i].q[2]/mesh.cells[i].q[0])^^2.0));
			double rho = mesh.cells[i].q[0];
			a[i] = sqrt(gamma*p/rho);
		}
		else
		{
			a[i] = 0;//double.nan;
		}
	}
	return a;
}

double getEnergy(ref Cell cell)
{
	return cell.q[3];
}

double[][] getEnergy(ref Mesh mesh)
{
	double[][] e = new double[][](mesh.M, mesh.N);
	for(int i = 0; i < mesh.N; i++)
	{
		for(int j = 0; j < mesh.M; j++)
		{
			e[j][i] = mesh[i,j].q[3];
		}
	}
	return e;
}