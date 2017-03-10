/+ Copyright (c) 2017 Robert F. Rau II +/
module ebb.ms;

import std.algorithm : canFind, filter, sort;
import std.array : array, split;
import std.conv;
import std.file : dirEntries, DirEntry, mkdir, readText, SpanMode;
import std.getopt;
import std.math;
import std.random;
import std.stdio;
import std.string;

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
import ebb.solve;

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

unittest
{
	UMesh2 mesh;
	mesh.nodes ~= [0, 0];
	mesh.nodes ~= [3.5, 0];
	mesh.nodes ~= [3, 2.5];
	mesh.nodes ~= [1, 1.5];
	mesh.elements ~= [1, 2, 3, 4];
	mesh.q ~= Vector!4(0);
	mesh.cells = new UCell2[1];
	mesh.cells[0].nEdges = 4;
	mesh.buildMesh;

	double x1 = 0;
	double y1 = 1;
	assert(!pointIsInPolygon(x1, y1, mesh, mesh.cells[0]), "Point should be outside polygon");

	double x2 = 3;
	double y2 = 1;
	assert(pointIsInPolygon(x2, y2, mesh, mesh.cells[0]), "Point should be inside polygon");

	double x3 = 2;
	double y3 = 2.5;
	assert(!pointIsInPolygon(x3, y3, mesh, mesh.cells[0]), "Point should be outside polygon");

	double x4 = 4;
	double y4 = 1;
	assert(!pointIsInPolygon(x4, y4, mesh, mesh.cells[0]), "Point should be outside polygon");

	double x5 = 3;
	double y5 = 1;
	assert(pointIsInPolygon(x5, y5, mesh, mesh.cells[0]), "Point should be inside polygon");

	writeln("pointIsInPolygon test done");
}

bool pointIsInPolygon(double x, double y, ref UMesh2 mesh, ref UCell2 cell)
{
	uint edgeIntersections = 0;
	foreach(i; 0..cell.nEdges)
	{
		auto n1 = mesh.edges[cell.edges[i]].nodeIdx[0];
		auto n2 = mesh.edges[cell.edges[i]].nodeIdx[1];

		immutable double x1 = mesh.nodes[n1][0];
		immutable double y1 = mesh.nodes[n1][1];
		immutable double x2 = mesh.nodes[n2][0];
		immutable double y2 = mesh.nodes[n2][1];

		immutable double yMin = min(y1, y2);
		immutable double yMax = max(y1, y2);

		//writeln("yMin = ", yMin, "\t yMax = ", yMax);

		// test point must be within the y-axis bounds
		// of the edge
		if((yMin <= y) && (y <= yMax))
		{
			// special case verticle edge
			if(abs(x2 - x1) <= 1.0e-10)
			{
				if(x <= x1)
				{
					edgeIntersections++;
					continue;
				}
			}

			immutable double xMin = min(x1, x2);
			immutable double xMax = max(x1, x2);

			immutable double m = (y2 - y1)/(x2 - x1);
			immutable double b = y1 - m*x1;
			immutable double x_i = (1.0/m)*(y - b);

			//writeln("x_i = ", x_i, "\tm = ", m);

			// test point must be to the left of the intersection point
			// also if x_i is nan m was likely 0 so there would be no
			// intersection
			if((x <= x_i) && !x_i.isNaN)
			{
				// The intersection point must lie inbetween the edge
				// nodes
				if((xMin <= x_i) && (x_i <= xMax))
				{
					edgeIntersections++;
				}
			}
		}
	}

	// If there are an even number of edge intersections
	// than this point lies outside of the polygon.
	if(edgeIntersections%2 == 0)
	{
		return false;
	}
	else
	{
		return true;
	}
	
}

void addSourceTerm(ref Vector!4[] R, ref UMesh2 mesh, Config config)
{
	foreach(i; mesh.interiorCells)
	{
		double xMin = double.infinity;
		double xMax = -double.infinity;
		double yMin = double.infinity;
		double yMax = -double.infinity;

		// compute the bounding box around the polygon
		foreach(j; 0..mesh.cells[i].nEdges)
		{
			auto n1 = mesh.edges[mesh.cells[i].edges[j]].nodeIdx[0];
			auto n2 = mesh.edges[mesh.cells[i].edges[j]].nodeIdx[1];
			double tmpXmin = min(mesh.nodes[n1][0], mesh.nodes[n2][0]);
			double tmpYmin = min(mesh.nodes[n1][1], mesh.nodes[n2][1]);

			double tmpXmax = max(mesh.nodes[n1][0], mesh.nodes[n2][0]);
			double tmpYmax = max(mesh.nodes[n1][1], mesh.nodes[n2][1]);

			xMin = min(xMin, tmpXmin);
			xMax = max(xMax, tmpXmax);
			yMin = min(yMin, tmpYmin);
			yMax = max(yMax, tmpYmax);
		}

		const uint samplePoints = 100;
		uint currentSampledPoints = 0;
		auto I = Vector!4(0);
		//writeln("Adding source term to cell ", i);
		//writeln("yMin: ", yMin, "\tyMax: ", yMax);
		//writeln("xMin: ", xMin, "\txMax: ", xMax);
		//writeln;

		Mt19937 rando;
		while(currentSampledPoints < samplePoints)
		{
			double x = uniform(xMin, xMax, rando);
			double y = uniform(yMin, yMax, rando);
			if(pointIsInPolygon(x, y, mesh, mesh.cells[i]))
			{
				currentSampledPoints++;

				I += (1.0/samplePoints.to!double)*sourceTerm(x, y, config);
			}
			else
			{
				writeln("point (", x, ", ", y, ") not in poly");
			}
		}

		//I *= mesh.cells[i].area;

		R[i] += I;
		//R[i] -= I;
	}
}

double computeL2(ref Vector!4[] R, ref UMesh2 mesh, uint dim)
{
	double L2 = 0;
	foreach(i; mesh.interiorCells)
	{
		L2 += R[i][dim]^^2;
	}

	L2 = sqrt(L2);

	return L2;
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

void stepMesh(ref UMesh2 mesh, Config config, double t, double dt)
{
	try
	{
		mesh.comm = MPI_COMM_SELF;
		mesh.mpiRank = 0;
		// Setup bc's
		double M = config.ic[0];
		double aoa = config.ic[1] * (PI/180);
		double rho = config.ic[3];
		double U = 1.0;
		double a = U/M;
		double p = (rho*a^^2.0)/gamma;
		double u = U*cos(aoa);
		double v = U*sin(aoa);

		if(config.viscosity)
		{
			config.physicalConfig.mu = (rho*U*config.physicalConfig.L)/config.physicalConfig.Re;
			writeln("mu = ", config.physicalConfig.mu);
		}

		for(uint i = 0; i < mesh.bGroups.length; i++)
		{
			@nogc uint findBcIndex(string tag)
			{
				auto tagIdx = config.boundaries.countUntil!"a.bTag == b"(tag);
				if(tagIdx < 0)
				{
					char[64] str;
					str[] = '\0';
					str[0..tag.length] = tag[];
					printf("Could not find tag %s\n", str.ptr);
					enforce(false, "Could not find matching boundary condition tag");
				}
				return cast(uint)tagIdx;
			}

			uint bcIdx = findBcIndex(mesh.bTags[i]);

			for(uint j = 0; j < mesh.bGroups[i].length; j++)
			{
				enforce(mesh.edges[mesh.bGroups[i][j]].isBoundary, "Edge not boundary edge but should be");
				enforce(mesh.edges[mesh.bGroups[i][j]].boundaryTag == config.boundaries[bcIdx].bTag, "Incorrect boundary tag");

				mesh.edges[mesh.bGroups[i][j]].bIdx = bcIdx;

				M = config.boundaries[bcIdx].boundaryData[0];
				aoa = config.boundaries[bcIdx].boundaryData[1] * (PI/180);
				rho = config.boundaries[bcIdx].boundaryData[3];

				U = 1.0;
				a = U/M;
				p = (rho*a^^2.0)/gamma;
				u = U*cos(aoa);
				v = U*sin(aoa);

				mesh.edges[mesh.bGroups[i][j]].boundaryType = config.boundaries[bcIdx].type;

				if(mesh.edges[mesh.bGroups[i][j]].boundaryType == BoundaryType.FullState)
				{
					mesh.edges[mesh.bGroups[i][j]].q[1][0] = rho;
					mesh.edges[mesh.bGroups[i][j]].q[1][1] = rho*u;
					mesh.edges[mesh.bGroups[i][j]].q[1][2] = rho*v;
					mesh.edges[mesh.bGroups[i][j]].q[1][3] = p/(gamma - 1.0) + 0.5*rho*(u^^2.0 + v^^2.0);
				}
				else if(mesh.edges[mesh.bGroups[i][j]].boundaryType == BoundaryType.ConstPressure)
				{
					enforce(config.boundaries[bcIdx].boundaryData.length == 1, "Constant pressure boundary data should only have one element; the constant pressure");
					mesh.edges[mesh.bGroups[i][j]].bData = config.boundaries[bcIdx].boundaryData;
				}
				else if(mesh.edges[mesh.bGroups[i][j]].boundaryType == BoundaryType.Dirichlet)
				{
					auto mid = mesh.edges[mesh.bGroups[i][j]].mid;
					mesh.edges[mesh.bGroups[i][j]].q[1] = solution(mid[0], mid[1]);
					config.boundaries[bcIdx].dFunc = &solution;
				}
			}
		}

		foreach(i; mesh.interiorCells)
		{
			mesh.q[i] = solution(mesh.cells[i].centroid[0], mesh.cells[i].centroid[1]);
		}

		auto R = new Vector!4[mesh.interiorCells.length];
		final switch(config.limiter)
		{
			foreach(lim; limiterList)
			{
				case lim:
				{
					final switch(config.flux)
					{
						foreach(fl; fluxList)
						{
							case fl:
							{
								final switch(config.integrator)
								{
									foreach(inte; integratorList)
									{
										case inte:
											writeln("Running 2D finite volume solver");
											writeln("-limited: ", config.limited);
											writeln("-order: ", config.order);
											writeln("-limiter: "~lim);
											writeln("-flux: "~fl);
											writeln("-integrator: "~inte);
											double Rmax = -double.infinity;
											double newDt = double.infinity;
											dt = config.dt;
											uint iterations = 0;
											while((abs(Rmax) > 1.0e-8) && (iterations < 10000))
											{
												Rmax = 0;
												//solver(R, mesh.q, mesh, config, newDt, Rmax, true, true);
												ufvmSolver!(mixin(lim), mixin(fl), 4)(R, mesh.q, mesh, config, newDt, Rmax, config.limited, true);
												addSourceTerm(R, mesh, config);

												foreach(i; mesh.interiorCells)
												{
													if(!config.localTimestep)
													{
														mesh.q[i] = mesh.q[i] + dt*R[i];
													}
													else
													{
														mesh.q[i] = mesh.q[i] + mesh.cells[i].dt*R[i];
													}
												}

												dt = newDt;
												if(iterations % config.plotIter == 0)
												{
													writeln("Rmax = ", Rmax, " dt = ", dt);
												}

												iterations++;
											}
											auto rhoL2 = computeL2(R, mesh, 0);
											auto rhouL2 = computeL2(R, mesh, 1);
											auto rhovL2 = computeL2(R, mesh, 2);
											auto rhoEL2 = computeL2(R, mesh, 3);

											writeln("rho   L2: ", rhoL2);
											writeln("rho u L2: ", rhouL2);
											writeln("rho v L2: ", rhovL2);
											writeln("rho E L2: ", rhoEL2);

											break;
									}
								}
								break;
							}
						}
					}
					break;
				}
			}
		}
	}
	catch(CellException ce)
	{
		writeln("Solver encountered an error: ", ce.msg);
		writeln("	In file ", ce.file);
		writeln("	On line ", ce.line);
		MPI_COMM_WORLD.abort(1);
	}
	catch(EdgeException ex)
	{
		writeln("Solver encountered an error: ", ex.msg);
		writeln("	In file ", ex.file);
		writeln("	On line ", ex.line);
		writeln("	pL = ", getPressure(ex.edge.q[0]));
		writeln("	pR = ", getPressure(ex.edge.q[1]));
		writeln("	Flux = ", ex.edge.flux);
		writeln("	qL = ", ex.edge.q[0]);
		writeln("	qR = ", ex.edge.q[1]);
		writeln("	cell L = ", ex.edge.cellIdx[0]);
		writeln("	cell R = ", ex.edge.cellIdx[1]);
		writeln("	normal = ", ex.edge.normal);
		MPI_COMM_WORLD.abort(1);
	}
	catch(Exception ex)
	{
		writeln("Caught unknown exception. Message: ", ex.msg);
		writeln("	In file ", ex.file);
		writeln("	On line ", ex.line);
		MPI_COMM_WORLD.abort(1);
	}
	finally
	{

	}
}

int main(string[] args)
{
	if(args.length < 2)
	{
		writeln("Not enough input arguments");
		writeln("Usage: ebb-reconstruct configFile save_#_*");
		return -1;
	}

	init(args);

	string configFile = args[1];

	import std.path : dirName, asAbsolutePath;
	auto configStr = readText(configFile);
	auto config = loadConfig(configStr);
	auto meshFile = "/";
	meshFile = configFile.dirName.asAbsolutePath.to!string ~ meshFile;
	meshFile ~= config.meshFile;
	UMesh2 mesh = importMesh(meshFile);
	mesh.buildMesh;

	double dt, t;

	stepMesh(mesh, config, t, dt);

	writeln("exiting");

	return shutdown;
}
