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

Vector!4 sourceTerm(double x, double y, Config config)
{
	auto S = Vector!4(0);
	
	double mu = config.physicalConfig.mu;
	double Pr = config.physicalConfig.Pr;

	S[0] = br*cr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) + br*dr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) - bu*cu*sin(cu*x + du*y)*(ar + br*sin(cr*x + dr*y)) - bv*dv*sin(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y));

	S[1] = (mu*(br*cr*dr*sin(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) + bv*cv*dv*cos(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y)) + br*bv*cr*dv*cos(cr*x + dr*y)*sin(cv*x + dv*y) + br*bv*cv*dr*cos(cr*x + dr*y)*sin(cv*x + dv*y)))/(3*(ar + br*sin(cr*x + dr*y))) + bp*cp*cos(cp*x + dp*y) + (4*mu*(br*cr^^2*sin(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) + bu*cu^^2*cos(cu*x + du*y)*(ar + br*sin(cr*x + dr*y)) + 2*br*bu*cr*cu*cos(cr*x + dr*y)*sin(cu*x + du*y)))/(3*(ar + br*sin(cr*x + dr*y))) + (mu*(br*dr^^2*sin(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) + bu*du^^2*cos(cu*x + du*y)*(ar + br*sin(cr*x + dr*y)) + 2*br*bu*dr*du*cos(cr*x + dr*y)*sin(cu*x + du*y)))/(ar + br*sin(cr*x + dr*y)) + br*cr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y))^^2 + br*dr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y))*(av + bv*cos(cv*x + dv*y)) - 2*bu*cu*sin(cu*x + du*y)*(au + bu*cos(cu*x + du*y))*(ar + br*sin(cr*x + dr*y)) - bu*du*sin(cu*x + du*y)*(av + bv*cos(cv*x + dv*y))*(ar + br*sin(cr*x + dr*y)) - bv*dv*sin(cv*x + dv*y)*(au + bu*cos(cu*x + du*y))*(ar + br*sin(cr*x + dr*y)) - (br*dr^^2*mu*sin(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)))/(ar + br*sin(cr*x + dr*y)) - (br^^2*cr^^2*mu*cos(cr*x + dr*y)^^2*((4*au)/3 + (4*bu*cos(cu*x + du*y))/3))/(ar + br*sin(cr*x + dr*y))^^2 - (br*cr^^2*mu*sin(cr*x + dr*y)*((4*au)/3 + (4*bu*cos(cu*x + du*y))/3))/(ar + br*sin(cr*x + dr*y)) + (4*br*cr*mu*cos(cr*x + dr*y)*(br*cr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) - bu*cu*sin(cu*x + du*y)*(ar + br*sin(cr*x + dr*y))))/(3*(ar + br*sin(cr*x + dr*y))^^2) + (br*dr*mu*cos(cr*x + dr*y)*(br*cr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) - bv*cv*sin(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y))))/(ar + br*sin(cr*x + dr*y))^^2 - (2*br*cr*mu*cos(cr*x + dr*y)*(br*dr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) - bv*dv*sin(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y))))/(3*(ar + br*sin(cr*x + dr*y))^^2) + (br*dr*mu*cos(cr*x + dr*y)*(br*dr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) - bu*du*sin(cu*x + du*y)*(ar + br*sin(cr*x + dr*y))))/(ar + br*sin(cr*x + dr*y))^^2 - (br^^2*dr^^2*mu*cos(cr*x + dr*y)^^2*(au + bu*cos(cu*x + du*y)))/(ar + br*sin(cr*x + dr*y))^^2 - (br*cr*dr*mu*sin(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)))/(ar + br*sin(cr*x + dr*y)) + (br^^2*cr*dr*mu*cos(cr*x + dr*y)^^2*((2*av)/3 + (2*bv*cos(cv*x + dv*y))/3))/(ar + br*sin(cr*x + dr*y))^^2 + (br*cr*dr*mu*sin(cr*x + dr*y)*((2*av)/3 + (2*bv*cos(cv*x + dv*y))/3))/(ar + br*sin(cr*x + dr*y)) - (br^^2*cr*dr*mu*cos(cr*x + dr*y)^^2*(av + bv*cos(cv*x + dv*y)))/(ar + br*sin(cr*x + dr*y))^^2 - (4*br*bu*cr*cu*mu*cos(cr*x + dr*y)*sin(cu*x + du*y))/(3*(ar + br*sin(cr*x + dr*y))) - (br*bv*cr*dv*mu*cos(cr*x + dr*y)*sin(cv*x + dv*y))/(ar + br*sin(cr*x + dr*y)) + (2*br*bv*cv*dr*mu*cos(cr*x + dr*y)*sin(cv*x + dv*y))/(3*(ar + br*sin(cr*x + dr*y))) - (br*bu*dr*du*mu*cos(cr*x + dr*y)*sin(cu*x + du*y))/(ar + br*sin(cr*x + dr*y));

	S[2] = (mu*(br*cr*dr*sin(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) + bu*cu*du*cos(cu*x + du*y)*(ar + br*sin(cr*x + dr*y)) + br*bu*cr*du*cos(cr*x + dr*y)*sin(cu*x + du*y) + br*bu*cu*dr*cos(cr*x + dr*y)*sin(cu*x + du*y)))/(ar + br*sin(cr*x + dr*y)) - (2*mu*(br*cr*dr*sin(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) + bv*cv*dv*cos(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y)) + br*bv*cr*dv*cos(cr*x + dr*y)*sin(cv*x + dv*y) + br*bv*cv*dr*cos(cr*x + dr*y)*sin(cv*x + dv*y)))/(3*(ar + br*sin(cr*x + dr*y))) + bp*dp*cos(cp*x + dp*y) + (mu*(br*cr^^2*sin(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) + bv*cv^^2*cos(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y)) + 2*br*bv*cr*cv*cos(cr*x + dr*y)*sin(cv*x + dv*y)))/(ar + br*sin(cr*x + dr*y)) + (4*mu*(br*dr^^2*sin(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) + bv*dv^^2*cos(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y)) + 2*br*bv*dr*dv*cos(cr*x + dr*y)*sin(cv*x + dv*y)))/(3*(ar + br*sin(cr*x + dr*y))) + br*dr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y))^^2 + br*cr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y))*(av + bv*cos(cv*x + dv*y)) - bu*cu*sin(cu*x + du*y)*(av + bv*cos(cv*x + dv*y))*(ar + br*sin(cr*x + dr*y)) - bv*cv*sin(cv*x + dv*y)*(au + bu*cos(cu*x + du*y))*(ar + br*sin(cr*x + dr*y)) - 2*bv*dv*sin(cv*x + dv*y)*(av + bv*cos(cv*x + dv*y))*(ar + br*sin(cr*x + dr*y)) - (br*cr^^2*mu*sin(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)))/(ar + br*sin(cr*x + dr*y)) - (br^^2*dr^^2*mu*cos(cr*x + dr*y)^^2*((4*av)/3 + (4*bv*cos(cv*x + dv*y))/3))/(ar + br*sin(cr*x + dr*y))^^2 - (br*dr^^2*mu*sin(cr*x + dr*y)*((4*av)/3 + (4*bv*cos(cv*x + dv*y))/3))/(ar + br*sin(cr*x + dr*y)) + (br*cr*mu*cos(cr*x + dr*y)*(br*cr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) - bv*cv*sin(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y))))/(ar + br*sin(cr*x + dr*y))^^2 - (2*br*dr*mu*cos(cr*x + dr*y)*(br*cr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) - bv*cv*sin(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y))))/(3*(ar + br*sin(cr*x + dr*y))^^2) + (br*cr*mu*cos(cr*x + dr*y)*(br*dr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) - bu*du*sin(cu*x + du*y)*(ar + br*sin(cr*x + dr*y))))/(ar + br*sin(cr*x + dr*y))^^2 + (4*br*dr*mu*cos(cr*x + dr*y)*(br*dr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) - bv*dv*sin(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y))))/(3*(ar + br*sin(cr*x + dr*y))^^2) - (br^^2*cr^^2*mu*cos(cr*x + dr*y)^^2*(av + bv*cos(cv*x + dv*y)))/(ar + br*sin(cr*x + dr*y))^^2 - (br*cr*dr*mu*sin(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)))/(ar + br*sin(cr*x + dr*y)) + (br^^2*cr*dr*mu*cos(cr*x + dr*y)^^2*((2*au)/3 + (2*bu*cos(cu*x + du*y))/3))/(ar + br*sin(cr*x + dr*y))^^2 + (br*cr*dr*mu*sin(cr*x + dr*y)*((2*au)/3 + (2*bu*cos(cu*x + du*y))/3))/(ar + br*sin(cr*x + dr*y)) - (br^^2*cr*dr*mu*cos(cr*x + dr*y)^^2*(au + bu*cos(cu*x + du*y)))/(ar + br*sin(cr*x + dr*y))^^2 - (br*bv*cr*cv*mu*cos(cr*x + dr*y)*sin(cv*x + dv*y))/(ar + br*sin(cr*x + dr*y)) + (2*br*bu*cr*du*mu*cos(cr*x + dr*y)*sin(cu*x + du*y))/(3*(ar + br*sin(cr*x + dr*y))) - (br*bu*cu*dr*mu*cos(cr*x + dr*y)*sin(cu*x + du*y))/(ar + br*sin(cr*x + dr*y)) - (4*br*bv*dr*dv*mu*cos(cr*x + dr*y)*sin(cv*x + dv*y))/(3*(ar + br*sin(cr*x + dr*y)));

	S[3] = (au + bu*cos(cu*x + du*y))*
			(bp*cp*cos(cp*x + dp*y) - (ar/2 + (br*sin(cr*x + dr*y))/2)*
			(2*bu*cu*sin(cu*x + du*y)*(au + bu*cos(cu*x + du*y)) + 2*bv*cv*sin(cv*x + dv*y)*(av + bv*cos(cv*x + dv*y))) + 
			(bp*cp*cos(cp*x + dp*y))/(gamma - 1) + (br*cr*cos(cr*x + dr*y)*((au + bu*cos(cu*x + du*y))^^2 + (av + bv*cos(cv*x + dv*y))^^2))/2) + 
			(av + bv*cos(cv*x + dv*y))*(bp*dp*cos(cp*x + dp*y) - (ar/2 + (br*sin(cr*x + dr*y))/2)*(2*bu*du*sin(cu*x + du*y)*(au + bu*cos(cu*x + du*y)) + 
			2*bv*dv*sin(cv*x + dv*y)*(av + bv*cos(cv*x + dv*y))) + (bp*dp*cos(cp*x + dp*y))/(gamma - 1) + (br*dr*cos(cr*x + dr*y)*
			((au + bu*cos(cu*x + du*y))^^2 + (av + bv*cos(cv*x + dv*y))^^2))/2) + (mu*(au + bu*cos(cu*x + du*y))*(br*cr*dr*sin(cr*x + dr*y)*
			(av + bv*cos(cv*x + dv*y)) + bv*cv*dv*cos(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y)) + br*bv*cr*dv*cos(cr*x + dr*y)*sin(cv*x + dv*y) + 
			br*bv*cv*dr*cos(cr*x + dr*y)*sin(cv*x + dv*y)))/(ar + br*sin(cr*x + dr*y)) - (mu*((2*av)/3 + (2*bv*cos(cv*x + dv*y))/3)*(br*cr*dr*sin(cr*x + dr*y)*
			(au + bu*cos(cu*x + du*y)) + bu*cu*du*cos(cu*x + du*y)*(ar + br*sin(cr*x + dr*y)) + br*bu*cr*du*cos(cr*x + dr*y)*sin(cu*x + du*y) + 
			br*bu*cu*dr*cos(cr*x + dr*y)*sin(cu*x + du*y)))/(ar + br*sin(cr*x + dr*y)) - (mu*((2*au)/3 + (2*bu*cos(cu*x + du*y))/3)*
			(br*cr*dr*sin(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) + bv*cv*dv*cos(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y)) + 
			br*bv*cr*dv*cos(cr*x + dr*y)*sin(cv*x + dv*y) + br*bv*cv*dr*cos(cr*x + dr*y)*sin(cv*x + dv*y)))/(ar + br*sin(cr*x + dr*y)) - 
			bu*cu*sin(cu*x + du*y)*(ap + bp*sin(cp*x + dp*y) + (ar/2 + (br*sin(cr*x + dr*y))/2)*((au + bu*cos(cu*x + du*y))^^2 + 
			(av + bv*cos(cv*x + dv*y))^^2) + (ap + bp*sin(cp*x + dp*y))/(gamma - 1)) - bv*dv*sin(cv*x + dv*y)*(ap + bp*sin(cp*x + dp*y) + 
			(ar/2 + (br*sin(cr*x + dr*y))/2)*((au + bu*cos(cu*x + du*y))^^2 + (av + bv*cos(cv*x + dv*y))^^2) + (ap + bp*sin(cp*x + dp*y))/(gamma - 1)) - 
			(mu*(au + bu*cos(cu*x + du*y))*(gamma/Pr - 4/3)*(br*cr^^2*sin(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) + bu*cu^^2*cos(cu*x + du*y)*
			(ar + br*sin(cr*x + dr*y)) + 2*br*bu*cr*cu*cos(cr*x + dr*y)*sin(cu*x + du*y)))/(ar + br*sin(cr*x + dr*y)) - 
			(mu*(av + bv*cos(cv*x + dv*y))*(gamma/Pr - 1)*(br*cr^^2*sin(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) + bv*cv^^2*cos(cv*x + dv*y)*
			(ar + br*sin(cr*x + dr*y)) + 2*br*bv*cr*cv*cos(cr*x + dr*y)*sin(cv*x + dv*y)))/(ar + br*sin(cr*x + dr*y)) - 
			(mu*(au + bu*cos(cu*x + du*y))*(gamma/Pr - 1)*(br*dr^^2*sin(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) + bu*du^^2*cos(cu*x + du*y)*
			(ar + br*sin(cr*x + dr*y)) + 2*br*bu*dr*du*cos(cr*x + dr*y)*sin(cu*x + du*y)))/(ar + br*sin(cr*x + dr*y)) - 
			(mu*(av + bv*cos(cv*x + dv*y))*(gamma/Pr - 4/3)*(br*dr^^2*sin(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) + bv*dv^^2*cos(cv*x + dv*y)*
			(ar + br*sin(cr*x + dr*y)) + 2*br*bv*dr*dv*cos(cr*x + dr*y)*sin(cv*x + dv*y)))/(ar + br*sin(cr*x + dr*y)) + (mu*(au + bu*cos(cu*x + du*y))*
			(av + bv*cos(cv*x + dv*y))*(br*cr*dr*sin(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) + bu*cu*du*cos(cu*x + du*y)*(ar + br*sin(cr*x + dr*y)) + 
			br*bu*cr*du*cos(cr*x + dr*y)*sin(cu*x + du*y) + br*bu*cu*dr*cos(cr*x + dr*y)*sin(cu*x + du*y)))/(ar + br*sin(cr*x + dr*y)) + 
			(gamma*mu*(br*cr*cos(cr*x + dr*y)*(2*bu*cu*sin(cu*x + du*y)*(au + bu*cos(cu*x + du*y)) + 2*bv*cv*sin(cv*x + dv*y)*(av + bv*cos(cv*x + dv*y))) - 
			(ar/2 + (br*sin(cr*x + dr*y))/2)*(2*bu^^2*cu^^2*sin(cu*x + du*y)^^2 + 2*bv^^2*cv^^2*sin(cv*x + dv*y)^^2 - 2*bu*cu^^2*cos(cu*x + du*y)*
			(au + bu*cos(cu*x + du*y)) - 2*bv*cv^^2*cos(cv*x + dv*y)*(av + bv*cos(cv*x + dv*y))) + (br*cr^^2*sin(cr*x + dr*y)*((au + bu*cos(cu*x + du*y))^^2 + 
			(av + bv*cos(cv*x + dv*y))^^2))/2 + (bp*cp^^2*sin(cp*x + dp*y))/(gamma - 1)))/(Pr*(ar + br*sin(cr*x + dr*y))) + (gamma*mu*(br*dr*cos(cr*x + dr*y)*
			(2*bu*du*sin(cu*x + du*y)*(au + bu*cos(cu*x + du*y)) + 2*bv*dv*sin(cv*x + dv*y)*(av + bv*cos(cv*x + dv*y))) - (ar/2 + (br*sin(cr*x + dr*y))/2)*(2*bu^^2*du^^2*sin(cu*x + du*y)^^2 + 2*bv^^2*dv^^2*sin(cv*x + dv*y)^^2 - 2*bu*du^^2*cos(cu*x + du*y)*(au + bu*cos(cu*x + du*y)) - 2*bv*dv^^2*cos(cv*x + dv*y)*(av + bv*cos(cv*x + dv*y))) + (br*dr^^2*sin(cr*x + dr*y)*((au + bu*cos(cu*x + du*y))^^2 + (av + bv*cos(cv*x + dv*y))^^2))/2 + (bp*dp^^2*sin(cp*x + dp*y))/(gamma - 1)))/(Pr*(ar + br*sin(cr*x + dr*y))) + (br^^2*cr^^2*mu*cos(cr*x + dr*y)^^2*((au + bu*cos(cu*x + du*y))^^2*(gamma/Pr - 4/3) + (av + bv*cos(cv*x + dv*y))^^2*(gamma/Pr - 1) - (gamma*((au + bu*cos(cu*x + du*y))^^2/2 + (av + bv*cos(cv*x + dv*y))^^2/2 + (ap + bp*sin(cp*x + dp*y))/((gamma - 1)*(ar + br*sin(cr*x + dr*y)))))/Pr))/(ar + br*sin(cr*x + dr*y))^^2 + (br^^2*dr^^2*mu*cos(cr*x + dr*y)^^2*((au + bu*cos(cu*x + du*y))^^2*(gamma/Pr - 1) + (av + bv*cos(cv*x + dv*y))^^2*(gamma/Pr - 4/3) - (gamma*((au + bu*cos(cu*x + du*y))^^2/2 + (av + bv*cos(cv*x + dv*y))^^2/2 + (ap + bp*sin(cp*x + dp*y))/((gamma - 1)*(ar + br*sin(cr*x + dr*y)))))/Pr))/(ar + br*sin(cr*x + dr*y))^^2 + (br*cr^^2*mu*sin(cr*x + dr*y)*((au + bu*cos(cu*x + du*y))^^2*(gamma/Pr - 4/3) + (av + bv*cos(cv*x + dv*y))^^2*(gamma/Pr - 1) - (gamma*((au + bu*cos(cu*x + du*y))^^2/2 + (av + bv*cos(cv*x + dv*y))^^2/2 + (ap + bp*sin(cp*x + dp*y))/((gamma - 1)*(ar + br*sin(cr*x + dr*y)))))/Pr))/(ar + br*sin(cr*x + dr*y)) + (br*dr^^2*mu*sin(cr*x + dr*y)*((au + bu*cos(cu*x + du*y))^^2*(gamma/Pr - 1) + (av + bv*cos(cv*x + dv*y))^^2*(gamma/Pr - 4/3) - (gamma*((au + bu*cos(cu*x + du*y))^^2/2 + (av + bv*cos(cv*x + dv*y))^^2/2 + (ap + bp*sin(cp*x + dp*y))/((gamma - 1)*(ar + br*sin(cr*x + dr*y)))))/Pr))/(ar + br*sin(cr*x + dr*y)) + (br*cr*mu*cos(cr*x + dr*y)*(2*bu*cu*sin(cu*x + du*y)*(au + bu*cos(cu*x + du*y))*(gamma/Pr - 4/3) - (gamma*(bu*cu*sin(cu*x + du*y)*(au + bu*cos(cu*x + du*y)) + bv*cv*sin(cv*x + dv*y)*(av + bv*cos(cv*x + dv*y)) - (bp*cp*cos(cp*x + dp*y))/((gamma - 1)*(ar + br*sin(cr*x + dr*y))) + (br*cr*cos(cr*x + dr*y)*(ap + bp*sin(cp*x + dp*y)))/((gamma - 1)*(ar + br*sin(cr*x + dr*y))^^2)))/Pr + 2*bv*cv*sin(cv*x + dv*y)*(av + bv*cos(cv*x + dv*y))*(gamma/Pr - 1)))/(ar + br*sin(cr*x + dr*y)) - (2*bv*dv*mu*sin(cv*x + dv*y)*(br*cr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) - bu*cu*sin(cu*x + du*y)*(ar + br*sin(cr*x + dr*y))))/(3*(ar + br*sin(cr*x + dr*y))) + (bu*du*mu*sin(cu*x + du*y)*(br*cr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) - bv*cv*sin(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y))))/(ar + br*sin(cr*x + dr*y)) - (2*bu*cu*mu*sin(cu*x + du*y)*(br*dr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) - bv*dv*sin(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y))))/(3*(ar + br*sin(cr*x + dr*y))) + (br*dr*mu*cos(cr*x + dr*y)*(2*bu*du*sin(cu*x + du*y)*(au + bu*cos(cu*x + du*y))*(gamma/Pr - 1) - (gamma*(bu*du*sin(cu*x + du*y)*(au + bu*cos(cu*x + du*y)) + bv*dv*sin(cv*x + dv*y)*(av + bv*cos(cv*x + dv*y)) - (bp*dp*cos(cp*x + dp*y))/((gamma - 1)*(ar + br*sin(cr*x + dr*y))) + (br*dr*cos(cr*x + dr*y)*(ap + bp*sin(cp*x + dp*y)))/((gamma - 1)*(ar + br*sin(cr*x + dr*y))^^2)))/Pr + 2*bv*dv*sin(cv*x + dv*y)*(av + bv*cos(cv*x + dv*y))*(gamma/Pr - 4/3)))/(ar + br*sin(cr*x + dr*y)) + (br*dr*mu*cos(cr*x + dr*y)*(br*cr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) - bv*cv*sin(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y)))*(au + bu*cos(cu*x + du*y)))/(ar + br*sin(cr*x + dr*y))^^2 + (bu*cu*mu*sin(cu*x + du*y)*(br*dr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) - bu*du*sin(cu*x + du*y)*(ar + br*sin(cr*x + dr*y)))*(av + bv*cos(cv*x + dv*y)))/(ar + br*sin(cr*x + dr*y)) + (bv*cv*mu*sin(cv*x + dv*y)*(br*dr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) - bu*du*sin(cu*x + du*y)*(ar + br*sin(cr*x + dr*y)))*(au + bu*cos(cu*x + du*y)))/(ar + br*sin(cr*x + dr*y)) - (bu*cu*mu*sin(cu*x + du*y)*(br*cr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) - bu*cu*sin(cu*x + du*y)*(ar + br*sin(cr*x + dr*y)))*(gamma/Pr - 4/3))/(ar + br*sin(cr*x + dr*y)) - (bv*cv*mu*sin(cv*x + dv*y)*(br*cr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) - bv*cv*sin(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y)))*(gamma/Pr - 1))/(ar + br*sin(cr*x + dr*y)) - (bu*du*mu*sin(cu*x + du*y)*(br*dr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) - bu*du*sin(cu*x + du*y)*(ar + br*sin(cr*x + dr*y)))*(gamma/Pr - 1))/(ar + br*sin(cr*x + dr*y)) - (bv*dv*mu*sin(cv*x + dv*y)*(br*dr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) - bv*dv*sin(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y)))*(gamma/Pr - 4/3))/(ar + br*sin(cr*x + dr*y)) - (br*dr*mu*cos(cr*x + dr*y)*((2*av)/3 + (2*bv*cos(cv*x + dv*y))/3)*(br*cr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) - bu*cu*sin(cu*x + du*y)*(ar + br*sin(cr*x + dr*y))))/(ar + br*sin(cr*x + dr*y))^^2 - (br*cr*mu*cos(cr*x + dr*y)*((2*au)/3 + (2*bu*cos(cu*x + du*y))/3)*(br*dr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) - bv*dv*sin(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y))))/(ar + br*sin(cr*x + dr*y))^^2 + (br*cr*mu*cos(cr*x + dr*y)*(br*dr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) - bu*du*sin(cu*x + du*y)*(ar + br*sin(cr*x + dr*y)))*(au + bu*cos(cu*x + du*y))*(av + bv*cos(cv*x + dv*y)))/(ar + br*sin(cr*x + dr*y))^^2 + (br*cr*gamma*mu*cos(cr*x + dr*y)*((bp*cp*cos(cp*x + dp*y))/(gamma - 1) - (ar/2 + (br*sin(cr*x + dr*y))/2)*(2*bu*cu*sin(cu*x + du*y)*(au + bu*cos(cu*x + du*y)) + 2*bv*cv*sin(cv*x + dv*y)*(av + bv*cos(cv*x + dv*y))) + (br*cr*cos(cr*x + dr*y)*((au + bu*cos(cu*x + du*y))^^2 + (av + bv*cos(cv*x + dv*y))^^2))/2))/(Pr*(ar + br*sin(cr*x + dr*y))^^2) + (br*dr*gamma*mu*cos(cr*x + dr*y)*((bp*dp*cos(cp*x + dp*y))/(gamma - 1) - (ar/2 + (br*sin(cr*x + dr*y))/2)*(2*bu*du*sin(cu*x + du*y)*(au + bu*cos(cu*x + du*y)) + 2*bv*dv*sin(cv*x + dv*y)*(av + bv*cos(cv*x + dv*y))) + (br*dr*cos(cr*x + dr*y)*((au + bu*cos(cu*x + du*y))^^2 + (av + bv*cos(cv*x + dv*y))^^2))/2))/(Pr*(ar + br*sin(cr*x + dr*y))^^2) - (br*cr*mu*cos(cr*x + dr*y)*(br*cr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) - bu*cu*sin(cu*x + du*y)*(ar + br*sin(cr*x + dr*y)))*(au + bu*cos(cu*x + du*y))*(gamma/Pr - 4/3))/(ar + br*sin(cr*x + dr*y))^^2 - (br*cr*mu*cos(cr*x + dr*y)*(br*cr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) - bv*cv*sin(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y)))*(av + bv*cos(cv*x + dv*y))*(gamma/Pr - 1))/(ar + br*sin(cr*x + dr*y))^^2 - (2*br*cr*dr*mu*sin(cr*x + dr*y)*(au/3 + (bu*cos(cu*x + du*y))/3)*(av + bv*cos(cv*x + dv*y)))/(ar + br*sin(cr*x + dr*y)) - (br*dr*mu*cos(cr*x + dr*y)*(br*dr*cos(cr*x + dr*y)*(au + bu*cos(cu*x + du*y)) - bu*du*sin(cu*x + du*y)*(ar + br*sin(cr*x + dr*y)))*(au + bu*cos(cu*x + du*y))*(gamma/Pr - 1))/(ar + br*sin(cr*x + dr*y))^^2 - (br*dr*mu*cos(cr*x + dr*y)*(br*dr*cos(cr*x + dr*y)*(av + bv*cos(cv*x + dv*y)) - bv*dv*sin(cv*x + dv*y)*(ar + br*sin(cr*x + dr*y)))*(av + bv*cos(cv*x + dv*y))*(gamma/Pr - 4/3))/(ar + br*sin(cr*x + dr*y))^^2 - (2*br^^2*cr*dr*mu*cos(cr*x + dr*y)^^2*(au/3 + (bu*cos(cu*x + du*y))/3)*(av + bv*cos(cv*x + dv*y)))/(ar + br*sin(cr*x + dr*y))^^2 - (br*bv*cr*dv*mu*cos(cr*x + dr*y)*sin(cv*x + dv*y)*(au/3 + (bu*cos(cu*x + du*y))/3))/(ar + br*sin(cr*x + dr*y)) - (br*bv*cv*dr*mu*cos(cr*x + dr*y)*sin(cv*x + dv*y)*(au/3 + (bu*cos(cu*x + du*y))/3))/(ar + br*sin(cr*x + dr*y)) - (br*bu*cr*du*mu*cos(cr*x + dr*y)*sin(cu*x + du*y)*(av + bv*cos(cv*x + dv*y)))/(3*(ar + br*sin(cr*x + dr*y))) - (br*bu*cu*dr*mu*cos(cr*x + dr*y)*sin(cu*x + du*y)*(av + bv*cos(cv*x + dv*y)))/(3*(ar + br*sin(cr*x + dr*y)));
	return S;
}

@nogc Vector!4 solution(double x, double y)
{
	auto q = Vector!4(0);

	double rho = ar + br*sin(cr*x + dr*y);
	double u = au + bu*cos(cu*x + du*y);
	double v = av + bv*cos(cv*x + dv*y);
	double p = ap + bp*cos(cp*x + dp*y);

	q[0] = rho;
	q[1] = rho*u;
	q[2] = rho*v;
	q[3] = p/(gamma - 1) + 0.5*rho*(u^^2 + v^^2);

	return q;
}
