module ebb.integrators;

import core.stdc.stdio : printf;
 
import std.meta : aliasSeqOf;
import std.math;

import ebb.config;
import ebb.exception;
import ebb.mesh;

import numd.utility;
import numd.linearalgebra.matrix;

//alias integratorList = aliasSeqOf!(["BEuler", "Euler", "RK2", "RK2_TVD", "RK4"]);
alias integratorList = aliasSeqOf!(["Euler", "RK2", "RK2_TVD", "RK4"]);

@nogc Vector!4 abs(Vector!4 data)
{
	return Vector!4(std.math.abs(data[0]), std.math.abs(data[1]), std.math.abs(data[2]), std.math.abs(data[3]));
}

@nogc S sum(string member = "", string op = "", S = double, T)(T[] data)
{
	if(data.length <= 6)
	{
		mixin("S s = "~op~"(data[0]"~member~");");
		for(uint i = 1; i < data.length; i++)
		{
			mixin("s += "~op~"(data[i]"~member~");");
		}
		return s;
	}
	else
	{
		uint m = cast(uint)floor(cast(double)data.length/2.0);
		return sum!(member, op, S, T)(data[0..m-1]) + sum!(member, op, S, T)(data[m..$-1]);
	}
}

/++
	Backward Euler time integration scheme
+/
struct BEuler
{
	/+
	@nogc double secantRoot(RootEquation eqn, double xk)
	{
		int iterations = 0;
		double xkNext;
		double xkPrev = xk - 1;
		double f = 1.0;
		
		while(abs(f.re) > Tolerance)
		{
			f = eqn(xk);
			auto fp = eqn(xkPrev);
			xkNext = xk - f*((xk - xkPrev)/(f - fp));
			xkPrev = xk;
			xk = xkNext;
			//writeln(f.re);
			iterations++;
		}
		//writeln("Secant root finder converged in ", iterations, " iterations.");
		return xk;
	}
	+/
	private static bool initialized = false;
	private static Vector!4[] qNext;
	private static Vector!4[] qLast;
	private static Vector!4[] tmp;
	private static Vector!4[] f1;
	private static Vector!4[] f2;

	@nogc static void init(ref UMesh2 mesh)
	{
		import std.experimental.allocator.mallocator : Mallocator;
		if(!initialized)
		{
			tmp = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
			f1 = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
			f2 = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
			qNext = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
			qLast = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
			initialized = true;
		}
	}

	@nogc static ~this()
	{
		if(initialized)
		{
			import std.experimental.allocator.mallocator : Mallocator;
			Mallocator.instance.deallocate(tmp);
			Mallocator.instance.deallocate(f1);
			Mallocator.instance.deallocate(f2);
			Mallocator.instance.deallocate(qLast);
			Mallocator.instance.deallocate(qNext);
			initialized = false;
		}
	}

	@nogc static void step(alias solver)(Vector!4[] R, ref Vector!4[] q, ref UMesh2 mesh, Config config, ref double dt, ref double Rmax, SolverException ex)
	{
		double newDt = double.infinity;
		double maxDiff = -double.infinity;
		immutable double tol = 1e-10;
		uint iterations = 0;

		//for(uint i = 0; i < q.length; i++)
		foreach(i; mesh.interiorCells)
		{
			for(uint j = 0; j < 4; j++)
			{
				qLast[i][j] = mesh.q[i][j];
			}
		}		
		//qLast[] = mesh.q[];

		printf("dt = %.10e\n", dt);
		dt = dt/config.CFL*0.5;
		Euler.step!solver(R, q, mesh, config, dt, Rmax, ex);
		printf("new dt = %.10e\n", dt);

		@nogc void func(Vector!4[] f, Vector!4[] qn)
		{
			printf("in func\n");
			solver(R, qn, mesh, config, newDt, Rmax, true, false, ex);
			//for(uint i = 0; i < f.length; i++)
			foreach(i; mesh.interiorCells)
			{
				f[i] = qn[i] - mesh.q[i] - dt*R[i];
			}
		}

		while(std.math.abs(maxDiff) > tol)
		{
			func(f1, q);
			func(f2, qLast);
			foreach(i; mesh.interiorCells)
			{
				for(uint j = 0; j < 4; j++)
				{
					if(q[i][j] == qLast[i][j])
					{
						//printf("qs are the same\n");
						qNext[i][j] = q[i][j];
					}
					else
					{
						qNext[i][j] = q[i][j] - f1[i][j]*((q[i][j] - qLast[i][j])/(f1[i][j] - f2[i][j]));
					}
					/+
					if((f1[i][j] - f2[i][j]) == 0.0)
					{
						printf("zero denominator\n");
					}
					+/
					//qNext[i][j] = q[i][j] - f1[i][j]*((q[i][j] - qLast[i][j])/(f1[i][j] - f2[i][j]));
				}
			}

			qLast[] = q[];
			q[] = qNext[];
			iterations++;
			//maxDiff = f1.sum!("", "abs", Vector!4)[].sum;
			maxDiff = -double.infinity;
			foreach(i; mesh.interiorCells)
			{
				for(uint j = 0; j < 4; j++)
				{
					maxDiff = fmax(maxDiff, std.math.abs(f1[i][j]));
				}
			}

			printf("maxDiff = %.10e\n", maxDiff);
			printf("iterations = %d\n", iterations);
			if(iterations > 1000)
			{
				printf("1000 iterations!\n");
				break;
			}
/+
			f = eqn(xk);
			auto fp = eqn(xkPrev);
			xkNext = xk - f*((xk - xkPrev)/(f - fp));
			xkPrev = xk;
			xk = xkNext;
			//writeln(f.re);
			iterations++;
+/
		}

		printf("beuler step took %d iterations\n", iterations);
		/+
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
		+/
		dt = newDt;
	}
}

/++
	Forward Euler time integration scheme
+/
struct Euler
{
	@nogc static void init(ref UMesh2 mesh)
	{

	}

	@nogc static void step(alias solver)(Vector!4[] R, ref Vector!4[] q, ref UMesh2 mesh, Config config, ref double dt, ref double Rmax, SolverException ex)
	{
		double newDt = double.infinity;

		solver(R, mesh.q, mesh, config, newDt, Rmax, true, true, ex);

		foreach(i; mesh.interiorCells)
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

/++
	Fourth order Runge-Kutta time integration scheme
+/
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

	@nogc static void step(alias solver)(Vector!4[] R, ref Vector!4[] q, ref UMesh2 mesh, Config config, ref double dt, ref double Rmax, SolverException ex)
	{
		import core.stdc.stdio : printf;

		double newDt = double.infinity;

		solver(k1, mesh.q, mesh, config, newDt, Rmax, true, true, ex);

		foreach(i; mesh.interiorCells)
		{
			tmp[i] = mesh.q[i] + ((dt/2.0)*k1[i]);
		}
		solver(k2, tmp, mesh, config, newDt, Rmax, config.multistageLimiting, false, ex);

		foreach(i; mesh.interiorCells)
		{
			tmp[i] = mesh.q[i] + ((dt/2.0)*k2[i]);
		}
		solver(k3, tmp, mesh, config, newDt, Rmax, config.multistageLimiting, false, ex);

		foreach(i; mesh.interiorCells)
		{
			tmp[i] = mesh.q[i] + (dt*k3[i]);
		}
		solver(k4, tmp, mesh, config, newDt, Rmax, config.multistageLimiting, false, ex);

		foreach(i; mesh.interiorCells)
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

struct RK2_TVD
{
	private static bool initialized = false;
	private static Vector!4[] tmp;
	private static Vector!4[] qFE;
	private static Vector!4[] k2;

	@nogc static void init(ref UMesh2 mesh)
	{
		import std.experimental.allocator.mallocator : Mallocator;
		if(!initialized)
		{
			tmp = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
			qFE = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
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
			Mallocator.instance.deallocate(qFE);
			Mallocator.instance.deallocate(k2);
			initialized = false;
		}
	}

	@nogc static void step(alias solver)(Vector!4[] R, ref Vector!4[] q, ref UMesh2 mesh, Config config, ref double dt, ref double Rmax, SolverException ex)
	{
		import core.stdc.stdio : printf;

		double newDt = double.infinity;

		//Euler.step!solver(R, qFE, mesh, config, newDt, Rmax, ex);
		solver(R, mesh.q, mesh, config, newDt, Rmax, true, true, ex);

		foreach(i; mesh.interiorCells)
		{
			if(!config.localTimestep)
			{
				qFE[i] = mesh.q[i] + dt*R[i];
			}
			else
			{
				qFE[i] = mesh.q[i] + mesh.cells[i].dt*R[i];
			}
		}

		//solver(k1, mesh.q, mesh, config, newDt, Rmax, ex);
		/*
		for(uint i = 0; i < tmp.length; i++)
		{
			tmp[i] = mesh.q[i] + (dt*k1[i]);
		}
		*/
		solver(k2, qFE, mesh, config, newDt, Rmax, false, false, ex);

		foreach(i; mesh.interiorCells)
		{
			if(!config.localTimestep)
			{
				q[i] = 0.5*(mesh.q[i] + qFE[i] + dt*k2[i]);
				//q[i] = mesh.q[i] + (dt/2.0)*(k1[i] + k2[i]);
			}
			else
			{
				q[i] = 0.5*(mesh.q[i] + qFE[i] + mesh.cells[i].dt*k2[i]);
				//q[i] = mesh.q[i] + (mesh.cells[i].dt/2.0)*(k1[i] + k2[i]);
			}
		}

		dt = newDt;
	}
}

/++
	Second order Runge-Kutta time integration scheme
+/
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

	@nogc static void step(alias solver)(Vector!4[] R, ref Vector!4[] q, ref UMesh2 mesh, Config config, ref double dt, ref double Rmax, SolverException ex)
	{
		import core.stdc.stdio : printf;

		double newDt = double.infinity;

		solver(k1, mesh.q, mesh, config, newDt, Rmax, true, true, ex);

		foreach(i; mesh.interiorCells)
		{
			tmp[i] = mesh.q[i] + (dt*k1[i]);
		}
		solver(k2, tmp, mesh, config, newDt, Rmax, config.multistageLimiting, false, ex);

		foreach(i; mesh.interiorCells)
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