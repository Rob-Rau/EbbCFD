/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.solver;

import core.atomic;
import core.stdc.stdio : fopen, fwrite, fopen, printf, snprintf;
import core.sys.posix.signal;

import ebb.config;
import ebb.euler;
import ebb.exception;
import ebb.flux;
import ebb.gas.config;
import ebb.integrate;
import ebb.math;
import ebb.mesh;
import ebb.io;
import ebb.mpid;
import ebb.solve;
import ebb.utility;

import std.conv;
import std.file;
import std.math;
import std.stdio;
import std.string;

static shared bool interrupted = false;

@nogc @system nothrow extern(C) void handle(int sig)
{
	printf("Signal received\n");
	interrupted.atomicStore(true);
}

@nogc bool nonsteadyConv(size_t size)(double maxTime, Solution!size sln)
{
	static currentTime = 0;
	currentTime += sln.dt;
	if(currentTime >= maxTime)
	{
		return true;
	}
	else
	{
		return interrupted.atomicLoad;
	}
}

@nogc bool steadyConv(size_t size, Communication)(Communication comm, Solution!size sln)
{
	double Rmax = -double.infinity;
	foreach(dQdt; sln.dQdt)
	{
		foreach(dq; dQdt)
		{
			Rmax = fmax(dq, Rmax);
		}
	}

	comm.allreduce(Rmax, Rmax, MPI_MAX);

	if(Rmax < 1.0e-12)
	{
		return true;
	}
	else
	{
		return interrupted.atomicLoad;
	}
}

int main(string[] args)
{
	signal(SIGINT, &handle);
	signal(SIGUSR1, &handle);
	
	uint cells = 100;

	auto mesh = Mesh(cells);

	string solver = "FiniteVolume";
	string integ = "Euler";

	bool steadyState = false;

	mixin(
	q{
		writeln("Running 2D finite volume solver");
		writeln("-solver: "~"{1}");
		writeln("-integrator: "~"{0}");

		alias SolverType = {1}!(solverParams(4UL, PhysicalDimensions.TwoD));
		alias IntegratorType = {0}!SolverType;
		alias PhysicsType = Physics!(solverParams(4UL, PhysicalDimensions.TwoD), GasPhysicalConfig);

		auto comm = Communication!(solverParams(4UL, PhysicalDimensions.TwoD))();
		auto config = GasPhysicalConfig();
		auto physics = new PhysicsType(["roeFlux", "averagedDiffusiveFlux"], config);
		auto slvr = SolverType(cast(immutable Mesh)mesh, physics, comm);
		auto integrator = IntegratorType(slvr, false);

		double startTime = comm.wtime;
		Solution!4U sln;
		if(steadyState)
		{
			sln = runIntegrator!(4U, (sln) => nonsteadyConv!4U(10.0, sln))(integrator, comm);
		}
		else
		{
			sln = runIntegrator!(4U, (sln) => steadyConv!4U(comm, sln))(integrator, comm);
		}
		
		double elapsed = comm.wtime - startTime;
		import std.format : format;
		writeln("elapsed time = ", format("%.16f", elapsed));
	}.
	switchBuilder!(1, "solver", solverList).
	switchBuilder!(0, "integ", integratorList));

	return 0;
}
