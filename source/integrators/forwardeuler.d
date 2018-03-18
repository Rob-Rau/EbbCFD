/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.integrators.forwardeuler;

import core.stdc.stdio : printf;

import ebb.config;
import ebb.exception;
import ebb.math;
import ebb.mesh;
import ebb.mpid;

import numd.linearalgebra.matrix;
import numd.utility;

import std.exception;
import std.math;

/++
	Forward Euler time integration scheme
+/
struct Euler(SpacialSolver)
{
	alias Sln = Solution!(solver.size);

	private SpacialSolver solver;

	private Vector!(solver.size)[] dQdt;
	private Vector!(solver.size)[] Q;

	private Sln sln;
	this(SpacialSolver solver, bool localTimestepping)
	{
		this.solver = solver;
		sln = new Sln();
		dQdt = new Vector!(solver.size)[solver.integrationSize];
		Q = new Vector!(solver.size)[solver.integrationSize];
		if(localTimestepping)
		{
			sln.dts = new double[solver.integrationSize];
		}
	}

	@nogc updateTimestep(bool localTimestepping)()
	{
		static if(localTimestepping)
		{
			solver.updateTimestep(sln.dts);
		}
		else
		{
			sln.dt = solver.updateTimestep;
		}
	}

	@nogc Sln step(bool localTimestepping, T...)(T args)
	{
		solver.solve!true(dQdt, Q, args);

		for(size_t i = 0; i < solver.integrationSize; i++)
		{
			static if(localTimestepping)
			{
				Q[i] = Q[i] + sln.dts[i]*dQdt[i];
			}
			else
			{
				Q[i] = Q[i] + sln.dt*dQdt[i];
			}
		}

		sln.dQdt = dQdt;
		sln.Q = Q;
		sln.iteration++;
		return sln;
	}
}
