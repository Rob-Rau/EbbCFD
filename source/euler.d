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

@nogc double getPressure(ref Vector!4 q)
{
	return (gamma - 1)*(q[3] - 0.5*q[0]*((q[1]/q[0])^^2.0 + (q[2]/q[0])^^2.0));
}
