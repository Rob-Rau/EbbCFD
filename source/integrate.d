/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.integrate;

import core.stdc.stdio : printf;

import ebb.config;
import ebb.exception;
public import ebb.integrators.forwardeuler;
public import ebb.integrators.rungekutta2;
public import ebb.integrators.rungekutta2tvd;
public import ebb.integrators.rungekutta4;
import ebb.math;
import ebb.mesh;
import ebb.mpid;
import ebb.solve;

import numd.linearalgebra.matrix;
import numd.utility;

import std.math;
import std.meta : aliasSeqOf;

alias integratorList = aliasSeqOf!(["Euler", "RK2", "RK2_TVD", "RK4"]);

@nogc Solution!size runIntegrator(size_t size, alias converged, TimeIntegrator, Communication, T...)(TimeIntegrator integrator, Communication comm, T args)
{
	Solution!size sln;
	bool localTimestep = false;

	if(localTimestep)
	{
		while(!converged(sln))
		{
			sln = integrator.step!true(args);
			integrator.updateTimestep!true();
		}
	}
	else
	{
		while(!converged(sln))
		{
			sln = integrator.step!false(args);
			integrator.updateTimestep!false();
			comm.allreduce(sln.dt, sln.dt, MPI_MAX);
		}
	}

	return sln;
}
