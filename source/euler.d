/+ Copyright (c) 2016 Robert F. Rau II +/
import std.math;

import numd.linearalgebra.matrix;

import mesh;

alias Vec = Vector!4;
alias Mat = Matrix!(4, 4);

T noop(T)(T x)
{
	return x;
}

Mat L(int kx, int ky)(Vec q)
{
	immutable double ρ = q[0];
	immutable double u = q[1]/ρ;
	immutable double v = q[2]/ρ;
	immutable double E = q[$-1];
	immutable double p = ((gamma - 1)*(E - 0.5*ρ*(u^^2.0 + v^^2.0)));
	immutable double a = sqrt(gamma*(p/ρ));
	return L!(kx, ky)(u, v, a);
}

Mat L(int kx, int ky)(double u, double v, double a)
{	
	immutable double Φ2 = 0.5*(gamma - 1)*(u^^2.0 + v^^2.0);
	immutable double β = 1.0/(2.0*a^^2.0);
	immutable double kxt = kx/sqrt(kx^^2.0 + ky^^2.0);
	immutable double kyt = ky/sqrt(kx^^2.0 + ky^^2.0);
	immutable double θ = kxt*v + kyt*v;

	Mat mat = Mat([1 - Φ2/a^^2.0, (gamma - 1)*u/a^^2.0, (gamma - 1)*v/a^^2.0, -(gamma - 1)/a^^2.0, 
					-(kyt*u - kxt*v), kyt, -kxt, 0,
					β*(Φ2 - a*θ), β*(kxt*a - (gamma - 1)*u), β*(kyt*a - (gamma - 1)*v), β*(gamma - 1),
					β*(Φ2 + a*θ), -β*(kxt*a + (gamma - 1)*u), -β*(kyt*a + (gamma - 1)*v), β*(gamma - 1)]);
	return mat;
}

Mat R(int kx, int ky)(Vec q)
{
	immutable double ρ = q[0];
	immutable double u = q[1]/ρ;
	immutable double v = q[2]/ρ;
	immutable double E = q[$-1];
	immutable double p = ((gamma - 1)*(E - 0.5*ρ*(u^^2.0 + v^^2.0)));
	immutable double a = sqrt(gamma*(p/ρ));
	
	return R!(kx, ky)(u, v, a);
}

Mat R(int kx, int ky)(double u, double v, double a)
{
	immutable double Φ2 = 0.5*(gamma - 1)*(u^^2.0 + v^^2.0);
	immutable double β = 1.0/(2.0*a^^2.0);
	immutable double kxt = kx/sqrt(kx^^2.0 + ky^^2.0);
	immutable double kyt = ky/sqrt(kx^^2.0 + ky^^2.0);
	immutable double θ = kxt*v + kyt*v;

	Mat mat = Mat([1, 0, 1, 1,
					u, kyt, u + kxt*a, u - kxt*a,
					v, -kxt, v + kyt*a, v - kyt*a,
					Φ2/(gamma - 1), kyt*u - kxt*v, (Φ2 + a^^2.0)/(gamma - 1) + a*θ, (Φ2 + a^^2.0)/(gamma - 1) - a*θ]);
	return mat;
}

Mat Lam(int kx, int ky, alias func = noop)(Vec q)
{
	immutable double ρ = q[0];
	immutable double u = q[1]/ρ;
	immutable double v = q[2]/ρ;
	immutable double E = q[$-1];
	immutable double p = ((gamma - 1)*(E - 0.5*ρ*(u^^2.0 + v^^2.0)));
	immutable double a = sqrt(gamma*(p/ρ));
	
	return Lam!(kx, ky, func)(u, v, a);
}
	
Mat Lam(int kx, int ky, alias func = noop)(double u, double v, double a)
{
	immutable double U = kx*u + ky*v;
	Mat mat = Mat([func(U), 0, 0, 0, 0, func(U), 0, 0, 0, 0, func(U + a), 0, 0, 0, 0, func(U - a)]);
	return mat;
}

immutable double gamma = 1.4;

double[][] getVelocity(int dim)(ref Mesh mesh)
{
	double[][] vel = new double[][](mesh.M, mesh.N);
	//double[][] vel = new double[][](mesh.N, mesh.M);
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