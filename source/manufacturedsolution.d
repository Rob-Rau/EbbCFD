/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.manufacturedsolution;

import core.atomic;
import core.stdc.stdio : fopen, fwrite, fopen, printf, snprintf;
import core.sys.posix.signal;

import std.algorithm : canFind, countUntil, min, map, max, reduce, sum;
import std.conv;
import std.file;
import std.math;
import std.stdio;

import numd.utility;
import numd.linearalgebra.matrix;

import ebb.config;
import ebb.euler;
import ebb.exception;
import ebb.flux;
import ebb.integrators;
import ebb.limiters;
import ebb.mesh;
import ebb.io;
import ebb.mpid;

immutable double ar = 0.9;
immutable double br = 0.04;
immutable double cr = -2;
immutable double dr = 1;
immutable double au = 0.1;
immutable double bu = 0.02;
immutable double cu = 1;
immutable double du = 1;
immutable double av = 0.05;
immutable double bv = -0.1;
immutable double cv = 0.7;
immutable double dv = 1.3;
immutable double ap = 1;
immutable double bp = 0.05;
immutable double cp = 2;
immutable double dp = -1;

@nogc double rho(double x, double y)
{
	return ar + br*sin(cr*x + dr*y);
}

@nogc double drdx(double x, double y)
{
	return br*cr*cos(cr*x + dr*y);
}

@nogc double d2rdx2(double x, double y)
{
	return -br*cr^^2.0*sin(cr*x + dr*y);
}

@nogc double drdy(double x, double y)
{
	return br*dr*cos(cr*x + dr*y);
}

@nogc double d2rdy2(double x, double y)
{
	return -br*dr^^2.0*sin(cr*x + dr*y);
}

@nogc double u(double x, double y)
{
	return au + bu*cos(cu*x + du*y);
}

@nogc double dudx(double x, double y)
{
	return -bu*cu*sin(cu*x + du*y);
}

@nogc double d2udxy(double x, double y)
{
	return -bu*cu*du*cos(cu*x + du*y);
}

@nogc double d2udx2(double x, double y)
{
	return -bu*cu^^2.0*cos(cu*x + du*y);
}

@nogc double dudy(double x, double y)
{
	return -bu*du*sin(cu*x + du*y);
}

@nogc double d2udy2(double x, double y)
{
	return -bu*du^^2.0*cos(cu*x + du*y);
}

@nogc double v(double x, double y)
{
	return av + bv*cos(cv*x + dv*y);
}

@nogc double dvdx(double x, double y)
{
	return -bv*cv*sin(cv*x + dv*y);
}

@nogc double d2vdx2(double x, double y)
{
	return -bv*cv^^2.0*cos(cv*x + dv*y);
}

@nogc double d2vdxy(double x, double y)
{
	return -bv*cv*dv*cos(cv*x + dv*y);
}

@nogc double dvdy(double x, double y)
{
	return -bv*dv*sin(cv*x + dv*y);
}

@nogc double d2vdy2(double x, double y)
{
	return -bv*dv^^2.0*cos(cv*x + dv*y);
}

@nogc double p(double x, double y)
{
	return ap + bp*sin(cp*x + dp*y);
}

@nogc double dpdx(double x, double y)
{
	return bp*cp*cos(cp*x + dp*y);
}

@nogc double d2pdx2(double x, double y)
{
	return -bp*cp^^2.0*sin(cp*x + dp*y);
}

@nogc double dpdy(double x, double y)
{
	return bp*dp*cos(cp*x + dp*y);
}

@nogc double d2pdy2(double x, double y)
{
	return -bp*dp^^2.0*sin(cp*x + dp*y);
}

Vector!4 sourceTerm(double x, double y, Config config)
{
	auto S = Vector!4(0);

	auto q = solution(x, y);
	auto dq = solutionGradient(x, y);

	auto dFxdx = Vector!4(0);
	auto dFydy = Vector!4(0);

	auto u = u(x, y);
	auto v = v(x, y);
	auto p = p(x, y);
	auto dudx = dudx(x, y);
	auto dvdy = dvdy(x, y);
	auto dpdx = dpdx(x, y);
	auto dpdy = dpdy(x, y);

	dFxdx[0] = dq[1,0];
	dFxdx[1] = 2*q[1]*dudx + u^^2*dq[0,0] + dpdx;
	dFxdx[2] = v*dFxdx[0] + q[1]*dvdx(x, y);
	dFxdx[3] = u*dq[3,0] + q[3]*dudx + u*dpdx + p*dudx;

	dFydy[0] = dq[2,1];
	dFydy[1] = v*dq[1,1] + q[1]*dvdy;
	dFydy[2] = 2*q[2]*dvdy + v^^2*dq[0,1] + dpdy;
	dFydy[3] = v*dq[3,1] + q[3]*dvdy + v*dpdy + p*dvdy;

	S = dFxdx + dFydy;
	
	if(config.viscosity)
	{
		auto dGxdx = Vector!4(0);
		auto dGydy = Vector!4(0);
		immutable double R = config.physicalConfig.R;
		immutable double mu = config.physicalConfig.mu;
		immutable double Pr = config.physicalConfig.Pr;
		immutable double k = (gamma*R*mu)/((gamma - 1)*Pr);

		immutable auto dudy = dudy(x, y);
		immutable auto d2udx2 = d2udx2(x, y);
		immutable auto d2udy2 = d2udy2(x, y);
		immutable auto d2udxy = d2udxy(x, y);

		immutable auto dvdx = dvdx(x, y);
		immutable auto d2vdx2 = d2vdx2(x, y);
		immutable auto d2vdy2 = d2vdy2(x, y);
		immutable auto d2vdxy = d2vdxy(x, y);

		immutable auto d2pdx2 = d2pdx2(x, y);
		immutable auto d2pdy2 = d2pdy2(x, y);
		
		immutable auto r = q[0];
		immutable auto drdx = dq[0,0];
		immutable auto drdy = dq[0,1];
		immutable auto d2rdx2 = d2rdx2(x, y);
		immutable auto d2rdy2 = d2rdy2(x, y);

		dGxdx[1] = mu*(2.0/3.0)*(2*d2udx2 - d2vdxy);
		dGxdx[2] = mu*(d2udxy + d2vdx2);

		dGxdx[3] = mu*(2.0/3.0)*(2*(u*d2udx2 + dudx^^2.0) - (u*d2vdxy + dudx*dvdy));
		dGxdx[3] += mu*((v*d2udxy + dvdx*dudy) + (v*d2vdx2 + dvdx^^2.0));
		immutable auto d2Tdx2 = (1.0/R)*((1.0/r)*d2pdx2 + (1.0/r^^2.0)*drdx*dpdx - (1.0/r^^3.0)*(r*dpdx - 2.0*p*drdx)*drdx - (p/r^^2.0)*d2rdx2);
		dGxdx[3] += k*d2Tdx2;

		dGydy[1] = mu*(d2udy2 + d2vdxy);
		dGydy[2] = mu*(2.0/3.0)*(2*d2vdy2 - d2udxy);
		dGydy[3] = mu*(2.0/3.0)*(2.0*(v*d2vdy2 + dvdy^^2.0) - (v*d2udxy + dvdy*dudx));
		dGydy[3] += mu*( u*d2udy2 + dudy^^2.0 + u*d2vdxy + dudy*dvdx);
		immutable auto d2Tdy2 = (1.0/R)*((1.0/r)*d2pdy2 + (1.0/r^^2.0)*drdy*dpdy - (1.0/r^^3.0)*(r*dpdy - 2.0*p*drdy)*drdy - (p/r^^2.0)*d2rdy2);
		dGydy[3] += k*d2Tdy2;
		
		S -= (dGxdx + dGydy);
	}

	return S;
}

@nogc Vector!4 solution(double x, double y)
{
	auto q = Vector!4(0);

	q[0] = rho(x, y);
	q[1] = rho(x, y)*u(x, y);
	q[2] = rho(x, y)*v(x, y);
	q[3] = p(x, y)/(gamma - 1) + 0.5*rho(x, y)*(u(x, y)^^2 + v(x, y)^^2);

	return q;
}

@nogc Matrix!(4,2) solutionGradient(double x, double y)
{
	auto grad = Matrix!(4,2)(0);

	// d(r)
	grad[0,0] = drdx(x, y);
	grad[0,1] = drdy(x, y);

	// d(ru)
	grad[1,0] = u(x, y)*drdx(x, y) + rho(x, y)*dudx(x, y);
	grad[1,1] = u(x, y)*drdy(x, y) + rho(x, y)*dudy(x, y);

	// d(rv)
	grad[2,0] = v(x, y)*drdx(x, y) + rho(x, y)*dvdx(x, y);
	grad[2,1] = v(x, y)*drdy(x, y) + rho(x, y)*dvdy(x, y);

	// d(rE)
	grad[3,0] = 1/(gamma - 1)*dpdx(x, y) + 0.5*(2*rho(x, y)*u(x, y)*dudx(x, y) + u(x, y)^^2*drdx(x, y) + 2*rho(x, y)*v(x, y)*dvdx(x, y) + v(x, y)^^2*drdx(x, y));
	grad[3,1] = 1/(gamma - 1)*dpdy(x, y) + 0.5*(2*rho(x, y)*u(x, y)*dudy(x, y) + u(x, y)^^2*drdy(x, y) + 2*rho(x, y)*v(x, y)*dvdy(x, y) + v(x, y)^^2*drdy(x, y));
	
	return grad;
}