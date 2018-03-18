/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.integrators.rungekutta2;

import core.stdc.stdio : printf;

import ebb.config;
import ebb.exception;
import ebb.integrators.forwardeuler;
import ebb.math;
import ebb.mesh;
import ebb.mpid;

import numd.linearalgebra.matrix;
import numd.utility;

import std.math;
import std.meta : aliasSeqOf;

/++
	Second order Runge-Kutta time integration scheme
+/
struct RK2(SpacialSolver)
{
	alias Sln = Solution!(solver.size);
	private SpacialSolver solver;
	private Sln sln;
	private Vector!(solver.size)[] Q;
	private Vector!(solver.size)[] tmp;
	private Vector!(solver.size)[] k1;
	private Vector!(solver.size)[] k2;

	this(SpacialSolver solver, bool localTimestepping)
	{
		this.solver = solver;
		sln = new Sln();
		Q = new Vector!(solver.size)[solver.integrationSize];
		tmp = new Vector!(solver.size)[solver.integrationSize];
		k1 = new Vector!(solver.size)[solver.integrationSize];
		k2 = new Vector!(solver.size)[solver.integrationSize];
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
		solver.solve!true(k1, Q, args);

		for(size_t i = 0; i < solver.integrationSize; i++)
		{
			static if(localTimestepping)
			{
				tmp[i] = Q[i] + (sln.dts[i]*k1[i]);
			}
			else
			{
				tmp[i] = Q[i] + (sln.dt*k1[i]);
			}
		}

		solver.solve!false(k2, tmp, args);

		for(size_t i = 0; i < solver.integrationSize; i++)
		{
			static if(localTimestepping)
			{
				Q[i] = Q[i] + (sln.dts[i]/2.0)*(k1[i] + k2[i]);
			}
			else
			{
				Q[i] = Q[i] + (sln.dt/2.0)*(k1[i] + k2[i]);
			}
		}

		sln.dQdt = k2;
		sln.Q = Q;
		sln.iteration++;
		return sln;
	}
}
