/+ Copyright (c) 2017 Robert F. Rau II +/
module ebb.ms;

import std.algorithm : canFind, filter, sort;
import std.array : array, split;
import std.conv;
import std.file : dirEntries, DirEntry, mkdir, readText, SpanMode;
import std.getopt;
import std.math;
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

Vector!4 sourceTerm(double x, double y, Config config)
{
	auto S = Vector!4(0);
	
	double mu = config.physicalConfig.mu;
	double Pr = config.physicalConfig.Pr;
	double ar = 0.9;
	double br = 0.04;
	double cr = -2;
	double dr = 1;

	double au = 0.1;
	double bu = 0.02;
	double cu = 1;
	double du = 1;
	
	double av = 0.05;
	double bv = -0.1;
	double cv = 0.7;
	double dv = 1.3;

	double ap = 1;
	double bp = 0.05;
	double cp = 2;
	double dp = -1;

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

void addSourceTerm(Vector!4[] R, ref UMesh2 mesh, Config config)
{
	foreach(i; mesh.interiorCells)
	{
	}
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
			}
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
											double Rmax;
											ufvmSolver!(mixin(lim), mixin(fl), 4)(R, mesh.q, mesh, config, t, Rmax, config.limited, false);
											addSourceTerm(R, mesh, config);
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
	//if(args.length < 3)
	//{
	//	writeln("Not enough input arguments");
	//	writeln("Usage: ebb-reconstruct configFile save_#_*");
	//	return -1;
	//}

	init(args);

	string configFile = args[1];

	//auto saveNum = slnFiles[0].split("_")[$-2];
	//writeln("saveNum: ", saveNum);
	//string savePrefix = "recon_save_"~saveNum;

	import std.path : dirName, asAbsolutePath;
	auto configStr = readText(configFile);
	auto config = loadConfig(configStr);
	auto meshFile = "/";
	meshFile = configFile.dirName.asAbsolutePath.to!string ~ meshFile;
	meshFile ~= config.meshFile;
	UMesh2 mesh = importMesh(meshFile);
	mesh.buildMesh;

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

	double dt, t;

	immutable size_t dims = 4;
	stepMesh(mesh, config, t, dt);
	Vector!4 S;
	double tmp = 0;
	int num = 1000000;
	foreach(i; 0..num)
	{
		S = sourceTerm(i.to!double, i.to!double, config);
		tmp += (S[0] + S[1] + S[2] + S[3])/num.to!double;
	}

	writeln(tmp);

	writeln("exiting");
	shutdown;
	return tmp.to!int;
}
