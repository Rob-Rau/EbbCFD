/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.integrators.rungekutta2tvd;

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
	This is a version of the Runge-Kutta scheme with total variation dimishing properties.
+/
struct RK2_TVD(SpacialSolver)
{
	alias Sln = Solution!(solver.size);
	private SpacialSolver solver;
	private Sln sln;
	private Vector!(solver.size)[] Q;
	private Vector!(solver.size)[] dQdt;
	private Vector!(solver.size)[] qFE;
	private Vector!(solver.size)[] k2;

	this(SpacialSolver solver, bool localTimestepping)
	{
		this.solver = solver;
		sln = new Sln();
		Q = new Vector!(solver.size)[solver.integrationSize];
		dQdt = new Vector!(solver.size)[solver.integrationSize];
		qFE = new Vector!(solver.size)[solver.integrationSize];
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
		solver.solve!true(dQdt, Q, args);

		for(size_t i = 0; i < solver.integrationSize; i++)
		{
			static if(localTimestepping)
			{
				qFE[i] = Q[i] + sln.dts[i]*dQdt[i];
			}
			else
			{
				qFE[i] = Q[i] + sln.dt*dQdt[i];
			}
		}

		solver.solve!false(k2, qFE, args);

		for(size_t i = 0; i < solver.integrationSize; i++)
		{
			static if(localTimestepping)
			{
				Q[i] = 0.5*(Q[i] + qFE[i] + sln.dts[i]*k2[i]);
			}
			else
			{
				Q[i] = 0.5*(Q[i] + qFE[i] + sln.dt*k2[i]);
			}
		}

		sln.dQdt = k2;
		sln.Q = Q;
		sln.iteration++;
		return sln;
	}
}
