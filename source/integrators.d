module ebb.integrators;

import std.meta : aliasSeqOf;

import ebb.config;
import ebb.exception;
import ebb.mesh;

import numd.utility;
import numd.linearalgebra.matrix;

alias integratorList = aliasSeqOf!(["Euler", "RK2", "RK4"]);

struct Euler
{
	@nogc static void init(ref UMesh2 mesh)
	{

	}

	@nogc static void step(alias solver)(Vector!4[] R, Vector!4[] q, ref UMesh2 mesh, Config config, ref double dt, ref double Rmax, SolverException ex)
	{
		double newDt = double.infinity;

		solver(R, mesh.q, mesh, config, newDt, Rmax, ex);

		for(uint i = 0; i < mesh.cells.length; i++)
		{
			if(!config.localTimestep)
			{
				q[i] = mesh.q[i] + dt*R[i];
			}
			else
			{
				q[i] = mesh.q[i] + mesh.cells[i].dt*R[i];
			}
		}

		dt = newDt;
	}
}

struct RK4
{
	private static bool initialized = false;
	private static Vector!4[] tmp;
	private static Vector!4[] k1;
	private static Vector!4[] k2;
	private static Vector!4[] k3;
	private static Vector!4[] k4;

	@nogc static void init(ref UMesh2 mesh)
	{
		import std.experimental.allocator.mallocator : Mallocator;
		if(!initialized)
		{
			tmp = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
			k1 = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
			k2 = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
			k3 = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
			k4 = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
			initialized = true;
		}
	}

	@nogc static ~this()
	{
		if(initialized)
		{
			import std.experimental.allocator.mallocator : Mallocator;
			Mallocator.instance.deallocate(tmp);
			Mallocator.instance.deallocate(k1);
			Mallocator.instance.deallocate(k2);
			Mallocator.instance.deallocate(k3);
			Mallocator.instance.deallocate(k4);
			initialized = false;
		}
	}

	@nogc static void step(alias solver)(Vector!4[] R, Vector!4[] q, ref UMesh2 mesh, Config config, ref double dt, ref double Rmax, SolverException ex)
	{
		import core.stdc.stdio : printf;

		double newDt = double.infinity;

		solver(k1, mesh.q, mesh, config, newDt, Rmax, ex);

		for(uint i = 0; i < tmp.length; i++)
		{
			tmp[i] = mesh.q[i] + ((dt/2.0)*k1[i]);
		}
		solver(k2, tmp, mesh, config, newDt, Rmax, ex);

		for(uint i = 0; i < tmp.length; i++)
		{
			tmp[i] = mesh.q[i] + ((dt/2.0)*k2[i]);
		}
		solver(k3, tmp, mesh, config, newDt, Rmax, ex);

		for(uint i = 0; i < tmp.length; i++)
		{
			tmp[i] = mesh.q[i] + (dt*k3[i]);
		}
		solver(k4, tmp, mesh, config, newDt, Rmax, ex);

		for(uint i = 0; i < mesh.q.length; i++)
		{
			if(!config.localTimestep)
			{
				q[i] = mesh.q[i] + (dt/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
			}
			else
			{
				q[i] = mesh.q[i] + (mesh.cells[i].dt/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
			}
		}

		dt = newDt;
	}
}

struct RK2
{
	private static bool initialized = false;
	private static Vector!4[] tmp;
	private static Vector!4[] k1;
	private static Vector!4[] k2;

	@nogc static void init(ref UMesh2 mesh)
	{
		import std.experimental.allocator.mallocator : Mallocator;
		if(!initialized)
		{
			tmp = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
			k1 = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
			k2 = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
			initialized = true;
		}
	}

	@nogc static ~this()
	{
		if(initialized)
		{
			import std.experimental.allocator.mallocator : Mallocator;
			Mallocator.instance.deallocate(tmp);
			Mallocator.instance.deallocate(k1);
			Mallocator.instance.deallocate(k2);
			initialized = false;
		}
	}

	@nogc static void step(alias solver)(Vector!4[] R, Vector!4[] q, ref UMesh2 mesh, Config config, ref double dt, ref double Rmax, SolverException ex)
	{
		import core.stdc.stdio : printf;

		double newDt = double.infinity;

		solver(k1, mesh.q, mesh, config, newDt, Rmax, ex);

		for(uint i = 0; i < tmp.length; i++)
		{
			tmp[i] = mesh.q[i] + (dt*k1[i]);
		}
		solver(k2, tmp, mesh, config, newDt, Rmax, ex);

		for(uint i = 0; i < mesh.q.length; i++)
		{
			if(!config.localTimestep)
			{
				q[i] = mesh.q[i] + (dt/2.0)*(k1[i] + k2[i]);
			}
			else
			{
				q[i] = mesh.q[i] + (mesh.cells[i].dt/2.0)*(k1[i] + k2[i]);
			}
		}

		dt = newDt;
	}
}