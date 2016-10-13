/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.solver;

import core.atomic;
import core.stdc.stdio : fopen, fwrite, fopen, printf, snprintf;
import core.sys.posix.signal;

import std.algorithm : canFind, min, max, reduce;
import std.conv;
import std.json;
import std.file;
import std.math;
import std.stdio;

import rpp.client.rpc;

import numd.utility;
import numd.linearalgebra.matrix;

import ebb.limiters;
import ebb.euler;
import ebb.flux;
import ebb.mesh;

static shared bool interupted = false;

@nogc nothrow @system extern(C) void handle(int sig)
{
	printf("Signal received\n");
	atomicStore(interupted, true);
}

class SolverException : Exception
{
	this(string msg)
	{
		super(msg);
	}

	enum SExceptionType
	{
		EdgeException,
		CellException
	}

	@nogc struct EdgeException
	{
		double pL, pR;
		Vector!4 flux, qL, qR;
		Vector!2 normal;
		uint cellL, cellR;

		@nogc this(double pL, double pR, Vector!4 flux, Vector!4 qL, Vector!4 qR, Vector!2 n, uint cellL, uint cellR)
		{
			this.pL = pL;
			this.pR = pR;
			this.flux = flux;
			this.qL = qL;
			this.qR = qR;
			this.normal = n;
			this.cellL = cellL;
			this.cellR = cellR;
		}
	}

	@nogc struct CellException
	{
		double p;
		Vector!4 q;
		uint cellIdx;

		@nogc this(double p, Vector!4 q, uint cellIdx)
		{
			this.p = p;
			this.q = q;
			this.cellIdx = cellIdx;
		}
	}

	SExceptionType exceptionType;
	CellException cExcept;
	EdgeException eExcept;

	@nogc void SetException(SExceptionType type, string msg, EdgeException ex, string file = __FILE__, size_t line = __LINE__)
	{
		exceptionType = type;
		eExcept = ex;
		this.msg = msg;
		this.file = file;
		this.line = line;
	}

	@nogc void SetException(SExceptionType type, string msg, CellException ex, string file = __FILE__, size_t line = __LINE__)
	{
		exceptionType = type;
		cExcept = ex;
		this.msg = msg;
		this.file = file;
		this.line = line;
	}

	/+
	printf("pL = %f\n", getPressure(mesh.q[i]));
	printf("qL = [%f, %f, %f, %f]\n", mesh.q[i][0], mesh.q[i][1], mesh.q[i][2], mesh.q[i][3]);
	printf("cell = %d\n", i);
	+/
	/+
	printf("pL = %f\n", getPressure(mesh.edges[i].q[0]));
	printf("pR = %f\n", getPressure(mesh.edges[i].q[1]));
	printf("Flux = [%f, %f, %f, %f]\n", mesh.edges[i].flux[0], mesh.edges[i].flux[1], mesh.edges[i].flux[2], mesh.edges[i].flux[3]);
	printf("qL = [%f, %f, %f, %f]\n", mesh.edges[i].q[0][0], mesh.edges[i].q[0][1], mesh.edges[i].q[0][2], mesh.edges[i].q[0][3]);
	printf("qR = [%f, %f, %f, %f]\n", mesh.edges[i].q[1][0], mesh.edges[i].q[1][1], mesh.edges[i].q[1][2], mesh.edges[i].q[1][3]);
	printf("cell L = %d\n", mesh.edges[i].cellIdx[0]);
	printf("cell R = %d\n", mesh.edges[i].cellIdx[1]);
	printf("normal = [%f, %f]\n", mesh.edges[i].normal[0], mesh.edges[i].normal[1]);
	+/
}

alias integratorList = aliasSeqOf!(["Euler", "RK2", "RK4"]);

struct Euler
{
	@nogc static void init(ref UMesh2 mesh)
	{

	}

	@nogc static void step(alias solver)(Vector!4[] R, ref UMesh2 mesh, Config config, ref double dt, ref double Rmax, SolverException ex)
	{
		double newDt = double.infinity;

		solver(R, mesh.q, mesh, config, newDt, Rmax, ex);

		for(uint i = 0; i < mesh.cells.length; i++)
		{
			mesh.q[i] = mesh.q[i] + dt*R[i];
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

	@nogc static void step(alias solver)(Vector!4[] R, ref UMesh2 mesh, Config config, ref double dt, ref double Rmax, SolverException ex)
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
			mesh.q[i] = mesh.q[i] + (dt/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
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

	@nogc static void step(alias solver)(Vector!4[] R, ref UMesh2 mesh, Config config, ref double dt, ref double Rmax, SolverException ex)
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
			mesh.q[i] = mesh.q[i] + (dt/2.0)*(k1[i] + k2[i]);
		}

		dt = newDt;
	}
}

@nogc void runIntegrator(alias setup, alias solver, alias integrator)(ref UMesh2 mesh, Config config, string saveFile, SolverException ex)
{
	import std.experimental.allocator.mallocator : Mallocator;
	import std.bitmanip : write;

	double residRhoLast = double.infinity;

	uint residRhoIncIters = 0;

	uint iterations = 0;
	uint saveItr = 0;
	double t = 0;
	double dt = config.dt;

	auto rotMat = Matrix!(2, 2)(cos(config.ic[1] * (PI/180)), -sin(config.ic[1] * (PI/180)), sin(config.ic[1] * (PI/180)), cos(config.ic[1] * (PI/180)));
	auto ld = Vector!2(0);

	immutable uint buffSize = 3*1024*1024*double.sizeof;
	size_t buffPos = 0;
	auto forceFile = fopen("boundaryForces.frc", "wb");

	double residRho = 0.0;
	double residU = 0.0;
	double residV = 0.0;
	double residE = 0.0;
	double residMax = 0.0;

	double[] lastRho = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] thisRho = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] lastU = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] thisU = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] lastV = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] thisV = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] lastE = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] thisE = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] tmp = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	ubyte[] forceBuffer = cast(ubyte[])Mallocator.instance.allocate(buffSize);
	Vector!4[] R = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);

	scope(exit) Mallocator.instance.deallocate(forceBuffer);
	scope(exit) Mallocator.instance.deallocate(R);
	scope(exit) Mallocator.instance.deallocate(lastRho);
	scope(exit) Mallocator.instance.deallocate(thisRho);
	scope(exit) Mallocator.instance.deallocate(lastU);
	scope(exit) Mallocator.instance.deallocate(thisU);
	scope(exit) Mallocator.instance.deallocate(lastV);
	scope(exit) Mallocator.instance.deallocate(thisV);
	scope(exit) Mallocator.instance.deallocate(lastE);
	scope(exit) Mallocator.instance.deallocate(thisE);
	scope(exit) Mallocator.instance.deallocate(tmp);

	// let the integrator do any neccessary initialization
	integrator.init(mesh);

	double lastRmax = 0.0;

	// Setup IC's and BC's
	setup(mesh, config, lastRho, lastU, lastV, lastE, t, dt, saveFile, ex);

	while(!approxEqual(t, config.tEnd) && !atomicLoad(interupted))
	{
		double Rmax = 0;

		integrator.step!solver(R, mesh, config, dt, Rmax, ex);

		for(uint i = 0; i < mesh.cells.length; i++)
		{
			thisRho[i] = mesh.q[i][0];
			thisU[i] = mesh.q[i][1];
			thisV[i] = mesh.q[i][2];
			thisE[i] = mesh.q[i][3];
			
			if(mesh.q[i][0].isNaN || mesh.q[i][1].isNaN || mesh.q[i][2].isNaN || mesh.q[i][3].isNaN)
			{
				printf("pL = %f\n", getPressure(mesh.q[i]));
				printf("qL = [%f, %f, %f, %f]\n", mesh.q[i][0], mesh.q[i][1], mesh.q[i][2], mesh.q[i][3]);
				printf("cell = %d\n", i);
				printf("iteration = %d\n", iterations);
				printf("dt = %f\n", dt);
				printf("t = %f\n", t);
				ex.msg = "Got nan on cell average value";
				ex.file = __FILE__;
				ex.line = __LINE__;
				throw ex;
			}
		}

		import std.algorithm : sum;

		tmp[] = (thisRho[] - lastRho[])^^2;
		residRho = sqrt(mesh.cells.length*tmp.sum)/thisRho.sum;
		tmp[] = (thisU[] - lastU[])^^2;
		residU = sqrt(mesh.cells.length*tmp.sum)/thisU.sum;
		tmp[] = (thisV[] - lastV[])^^2;
		residV = sqrt(mesh.cells.length*tmp.sum)/thisV.sum;
		tmp[] = (thisE[] - lastE[])^^2;
		residE = sqrt(mesh.cells.length*tmp.sum)/thisE.sum;

		lastRho[] = thisRho[];
		lastU[] = thisU[];
		lastV[] = thisV[];
		lastE[] = thisE[];

		auto f = mesh.computeBoundaryForces(config.forceBoundary);
		ld = rotMat*f;

		forceBuffer.write!double(t, &buffPos);
		forceBuffer.write!double(ld[0], &buffPos);
		forceBuffer.write!double(ld[1], &buffPos);

		if(buffPos == buffSize)
		{
			fwrite(forceBuffer.ptr, ubyte.sizeof, buffSize, forceFile);
			buffPos = 0;
		}

		residMax = max(residRho, residU, residV, residE);
		if(residRhoLast < residMax)
		{
			residRhoIncIters++;
		}
		else
		{
			residRhoIncIters = 0;
		}

		if(residRhoIncIters >= 2000)
		{
			config.CFL *= 0.99;
			printf("Max RMS residual hasn't decreased in 2000 iterations, decreasing CFL to %f\n", config.CFL);
			residRhoIncIters = 0;
		}

		residRhoLast = residMax;

		if(iterations % config.plotIter == 0)
		{
			printf("lift force = %f\t drag force = %f\t t = %f\n", ld[1], ld[0], t);
			printf("rho_RMS = %.10e\tu_RMS = %.10e\tv_RMS = %.10e\tE_RMS = %.10e\tFlux_R = %.10e\t dt = %10.10f\n", residRho, residU, residV, residE, Rmax, dt);
		}
		
		if(config.saveIter != -1)
		{
			if(iterations % config.saveIter == 0)
			{
				char[512] filename;
				filename[] = 0;
				snprintf(filename.ptr, 512, "save_%d.esln", saveItr);
				//saveMatlabSolution(mesh, filename.ptr);
				saveSolution(mesh, filename.ptr, t, dt);
				saveItr++;
			}
		}
		t += dt;
		iterations++;
		lastRmax = Rmax;
	}

	if(buffPos != 0)
	{
		fwrite(forceBuffer.ptr, ubyte.sizeof, buffPos, forceFile);
		buffPos = 0;
	}
	fclose(forceFile);

	saveSolution(mesh, cast(char*)"final.esln", t, dt);
	printf("lift force = %f\t drag force = %f\t t = %f\n", ld[1], ld[0], t);
	printf("rho_RMS = %.10e\tu_RMS = %.10e\tv_RMS = %.10e\tE_RMS = %.10e\tFlux_R = %.10e\t dt = %10.10f\n", residRho, residU, residV, residE, lastRmax, dt);

}

@nogc void ufvmSetup(ref UMesh2 mesh, Config config, double[] lastRho, double[] lastU, double[] lastV, double[] lastE, ref double t, ref double dt, string saveFile, SolverException ex)
{
	double M = config.ic[0];
	double aoa = config.ic[1] * (PI/180);
	double p = config.ic[2];
	double rho = config.ic[3];
	double a = sqrt(gamma*(p/rho));
	double U = M*a;
	double u = U*cos(aoa);
	double v = U*sin(aoa);

	printf("M = %f\tU = %f\n", M, U);

	if(saveFile == "")
	{
		// Setup initial conditions
		for(uint i = 0; i < mesh.cells.length; i++)
		{
			mesh.q[i] = buildQ(rho, u, v, p);

			lastRho[i] = mesh.q[i][0];
			lastU[i] = mesh.q[i][1];
			lastV[i] = mesh.q[i][2];
			lastE[i] = mesh.q[i][3];
		}
	}
	else
	{
		if(loadSolution(mesh, t, dt, saveFile))
		{
			for(uint i = 0; i < mesh.cells.length; i++)
			{
				lastRho[i] = mesh.q[i][0];
				lastU[i] = mesh.q[i][1];
				lastV[i] = mesh.q[i][2];
				lastE[i] = mesh.q[i][3];
			}
		}
		else
		{
			ex.msg = "Failed to load solution";
			ex.line = __LINE__;
			ex.file = __FILE__;
			throw ex;
		}
	}

	// Setup bc's
	for(uint i = 0; i < mesh.bGroups.length; i++)
	{
		@nogc uint findBcIndex(string tag)
		{
			for(uint j = 0; j < config.bTags.length; j++)
			{
				if(config.bTags[j] == tag)
				{
					return j;
				}
			}
			ex.msg = "Could not find match boundary condition tag";
			ex.file = __FILE__;
			ex.line = __LINE__;
			throw ex;
		}

		uint bcIdx = findBcIndex(mesh.bTags[i]);

		for(uint j = 0; j < mesh.bGroups[i].length; j++)
		{
			if(!mesh.edges[mesh.bGroups[i][j]].isBoundary)
			{
				ex.msg = "Edge not boundary edge but should be";
				ex.file = __FILE__;
				ex.line = __LINE__;
				throw ex;
			}

			if(mesh.edges[mesh.bGroups[i][j]].boundaryTag != config.bTags[bcIdx])
			{
				ex.msg = "Incorrect boundary tag";
				ex.file = __FILE__;
				ex.line = __LINE__;
				throw ex;
			}

			M = config.bc[bcIdx][0];
			aoa = config.bc[bcIdx][1] * (PI/180);
			p = config.bc[bcIdx][2];
			rho = config.bc[bcIdx][3];
			a = sqrt(gamma*(p/rho));
			U = M*a;
			u = U*cos(aoa);
			v = U*sin(aoa);
			mesh.edges[mesh.bGroups[i][j]].boundaryType = config.bTypes[bcIdx];

			if(mesh.edges[mesh.bGroups[i][j]].boundaryType == BoundaryType.FullState)
			{
				mesh.edges[mesh.bGroups[i][j]].q[1] = buildQ(rho, u, v, p);
			}
		}
	}
}

// Unstructured finite volume solver
@nogc void ufvmSolver(alias S, alias F, size_t dims)(ref Vector!4[] R, Vector!4[] q, ref UMesh2 mesh, Config config, ref double newDt, ref double Rmax, SolverException ex)
{
	// Build gradients
	if(config.order > 1)
	{
		for(uint i = 0; i < mesh.cells.length; i++)
		{
			Vector!6[4] du;
			for(uint j = 0; j < 4; j++)
			{
				du[j] = Vector!6(0); 
			}
			mesh.cells[i].minQ = q[i];
			mesh.cells[i].maxQ = q[i];
			mesh.cells[i].lim = Vector!4(1.0);
			for(uint j = 0; j < mesh.cells[i].nNeighborCells; j++)
			{
				uint idx = mesh.cells[i].neighborCells[j];
				
				for(uint k = 0; k < 4; k++)
				{
					mesh.cells[i].minQ[k] = fmin(mesh.cells[i].minQ[k], q[idx][k]);
					mesh.cells[i].maxQ[k] = fmax(mesh.cells[i].maxQ[k], q[idx][k]);

					du[k][j] = q[idx][k] - q[i][k];
				}
			}

			for(uint j = 0; j < 4; j++)
			{
				mesh.cells[i].gradient[j] = mesh.cells[i].gradMat*du[j];
			}

			if(config.limited)
			{
				for(uint j = 0; j < mesh.cells[i].nEdges; j++)
				{
					uint eIdx = mesh.cells[i].edges[j];
					auto qM = q[i];
					auto grad = mesh.cells[i].gradient;
					auto centroid = mesh.cells[i].centroid;
					auto mid = mesh.edges[eIdx].mid;
					auto dx = mid[0] - centroid[0];
					auto dy = mid[1] - centroid[1];
					auto minQ = mesh.cells[i].minQ;
					auto maxQ = mesh.cells[i].maxQ;
					
					auto qE = Vector!4(0);
					for(uint k = 0; k < dims; k++)
					{
						double s = 1.0;
						qE[k] = qM[k] + (grad[k][0]*dx + grad[k][1]*dy);

						if(qE[k] < minQ[k])
						{
							s = (minQ[k] - qM[k])/(grad[k][0]*dx + grad[k][1]*dy);
						}
						else if(qE[k] > maxQ[k])
						{
							s = (maxQ[k] - qM[k])/(grad[k][0]*dx + grad[k][1]*dy);
						}

						if(s < 0)
						{
							printf("computed negative limiter\n");
						}

						if(s > 1.0)
						{
							printf("computed limiter greater than 1.0\n");
						}
						
						mesh.cells[i].lim[k] = fmin(mesh.cells[i].lim[k], s);
					}
				}
			}
		}
	}

	for(uint i = 0; i < mesh.edges.length; i++)
	{
		if(mesh.edges[i].isBoundary)
		{
			switch(mesh.edges[i].boundaryType)
				with(BoundaryType)
			{
				case FullState:
					if(config.order == 1)
					{
						mesh.edges[i].q[0] = q[mesh.edges[i].cellIdx[0]];
					}
					else
					{
						auto qM = q[mesh.edges[i].cellIdx[0]];
						auto grad = mesh.cells[mesh.edges[i].cellIdx[0]].gradient;
						auto centroid = mesh.cells[mesh.edges[i].cellIdx[0]].centroid;
						auto mid = mesh.edges[i].mid;
						auto dx = mid[0] - centroid[0];
						auto dy = mid[1] - centroid[1];
						//auto lim = mesh.cells[mesh.edges[i].cellIdx[0]].lim;
						
						for(uint j = 0; j < dims; j++)
						{
							mesh.edges[i].q[0][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[0]].lim[j]*(grad[j][0]*dx + grad[j][1]*dy);
							//mesh.edges[i].q[0][j] = qM[j] + lim*(grad[j][0]*dx + grad[j][1]*dy);
						}

						if(getPressure(mesh.edges[i].q[0]) < 0)
						{
							double lim = 1.0;
							for(uint j = 0; j < dims; j++)
							{
								lim = fmin(lim, mesh.cells[mesh.edges[i].cellIdx[0]].lim[j]);
							}
							for(uint j = 0; j < dims; j++)
							{
								//mesh.edges[i].q[0][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[0]].lim[j]*(grad[j][0]*dx + grad[j][1]*dy);
								mesh.edges[i].q[0][j] = qM[j] + lim*(grad[j][0]*dx + grad[j][1]*dy);
							}
						}
					}

					auto qL = mesh.edges[i].q[0];
					auto qR = mesh.edges[i].q[1];

					mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax);
					if(mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN)
					{
						ex.msg = "Got nan on FullState boundary";
						ex.file = __FILE__;
						ex.line = __LINE__;
						throw ex;
					}
					break;
				case InviscidWall:
					if(config.order == 1)
					{
						mesh.edges[i].q[0] = q[mesh.edges[i].cellIdx[0]];
					}
					else
					{
						auto qM = q[mesh.edges[i].cellIdx[0]];
						auto grad = mesh.cells[mesh.edges[i].cellIdx[0]].gradient;
						auto centroid = mesh.cells[mesh.edges[i].cellIdx[0]].centroid;
						auto mid = mesh.edges[i].mid;
						auto dx = mid[0] - centroid[0];
						auto dy = mid[1] - centroid[1];
						//auto lim = mesh.cells[mesh.edges[i].cellIdx[0]].lim;
						
						for(uint j = 0; j < dims; j++)
						{
							mesh.edges[i].q[0][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[0]].lim[j]*(grad[j][0]*dx + grad[j][1]*dy);
							//mesh.edges[i].q[0][j] = qM[j] + lim*(grad[j][0]*dx + grad[j][1]*dy);
						}

						if(getPressure(mesh.edges[i].q[0]) < 0)
						{
							double lim = 1.0;
							for(uint j = 0; j < dims; j++)
							{
								lim = fmin(lim, mesh.cells[mesh.edges[i].cellIdx[0]].lim[j]);
							}

							for(uint j = 0; j < dims; j++)
							{
								//mesh.edges[i].q[0][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[0]].lim[j]*(grad[j][0]*dx + grad[j][1]*dy);
								mesh.edges[i].q[0][j] = qM[j] + lim*(grad[j][0]*dx + grad[j][1]*dy);
							}
						}
					}

					Vector!2 velP = (1/mesh.edges[i].q[0][0])*Vector!2(mesh.edges[i].q[0][1], mesh.edges[i].q[0][2]);
					auto vel = (velP - (velP.dot(mesh.edges[i].normal))*mesh.edges[i].normal).magnitude;
					double p = (gamma - 1)*(mesh.edges[i].q[0][3] - 0.5*mesh.edges[i].q[0][0]*vel^^2);
					double a = sqrt(gamma*(p/mesh.edges[i].q[0][0]));
					if(p < 0)
					{
						p = 1.0e-12;
						//printf("pressure less than 0 at wall\n");
					}
					mesh.edges[i].flux = Vector!4(0, p*mesh.edges[i].normal[0], p*mesh.edges[i].normal[1], 0);
					mesh.edges[i].sMax = std.math.abs(a);

					if(mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN)
					{
						ex.msg = "Got nan on wall boundary";
						ex.file = __FILE__;
						ex.line = __LINE__;
						throw ex;
					}
					break;
				default:
					ex.msg = "Unsupported boundary type";
					ex.file = __FILE__;
					ex.line = __LINE__;
					throw ex;
			}
		}
		else
		{
			if(config.order == 1)
			{
				mesh.edges[i].q[0] = q[mesh.edges[i].cellIdx[0]];
				mesh.edges[i].q[1] = q[mesh.edges[i].cellIdx[1]];
			}
			else
			{
				for(uint k = 0; k < 2; k++)
				{
					auto qM = q[mesh.edges[i].cellIdx[k]];
					auto grad = mesh.cells[mesh.edges[i].cellIdx[k]].gradient;
					auto centroid = mesh.cells[mesh.edges[i].cellIdx[k]].centroid;
					auto mid = mesh.edges[i].mid;
					auto dx = mid[0] - centroid[0];
					auto dy = mid[1] - centroid[1];
					//auto lim = mesh.cells[mesh.edges[i].cellIdx[k]].lim;
					
					for(uint j = 0; j < dims; j++)
					{
						mesh.edges[i].q[k][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[k]].lim[j]*(grad[j][0]*dx + grad[j][1]*dy);
						//mesh.edges[i].q[k][j] = qM[j] + lim*(grad[j][0]*dx + grad[j][1]*dy);
					}

					if(getPressure(mesh.edges[i].q[k]) < 0)
					{
						double lim = 1.0;
						for(uint j = 0; j < dims; j++)
						{
							lim = fmin(lim, mesh.cells[mesh.edges[i].cellIdx[k]].lim[j]);
						}

						for(uint j = 0; j < dims; j++)
						{
							//mesh.edges[i].q[k][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[k]].lim[j]*(grad[j][0]*dx + grad[j][1]*dy);
							mesh.edges[i].q[k][j] = qM[j] + lim*(grad[j][0]*dx + grad[j][1]*dy);
						}
					}
				}
			}

			auto qL = mesh.edges[i].q[0];
			auto qR = mesh.edges[i].q[1];

			mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax);

			if(mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN)
			{
				ex.SetException(SolverException.SExceptionType.EdgeException,
								"Got NaN on interior edge",
								SolverException.EdgeException(getPressure(mesh.edges[i].q[0]), 
															  getPressure(mesh.edges[i].q[1]),
															  mesh.edges[i].flux,
															  mesh.edges[i].q[0],
															  mesh.edges[i].q[1],
															  mesh.edges[i].normal,
															  mesh.edges[i].cellIdx[0],
															  mesh.edges[i].cellIdx[1]));
				throw ex;
			}
		}
	}

	newDt = double.infinity;

	for(uint i = 0; i < mesh.cells.length; i++)
	{
		R[i] = Vector!4(0);
		double sAve = 0;
		// integrate fluxes over cell edges
		for(uint j = 0; j < mesh.cells[i].nEdges; j++)
		{
			R[i] += mesh.cells[i].fluxMultiplier[j]*mesh.edges[mesh.cells[i].edges[j]].len*mesh.edges[mesh.cells[i].edges[j]].flux;
			sAve += mesh.edges[mesh.cells[i].edges[j]].len*mesh.edges[mesh.cells[i].edges[j]].sMax;
		}
		sAve /= mesh.cells[i].perim;

		newDt = fmin(newDt, (config.CFL*mesh.cells[i].d)/sAve);

		for(uint j = 0; j < dims; j++)
		{
			if(std.math.abs(R[i][j]) > Rmax)
			{
				Rmax = std.math.abs(R[i][j]);
			}
		}

		R[i] *= -(1.0/mesh.cells[i].area);
	}

}

/+
{
	"mesh": "../box.mesh",
	"dt": 0.01,
	"tEnd": 750.0,
	"limiter": "minmodS",
	"flux": "rusanovFlux",
	"plotIter": 50,
	"saveIter": 50,
	"initialConditions": [0.85, 0, 1, 1.4], (M, aoa, p, rho)
	"boudaryConditions": [
		{
			"tag": "Bottom",
			"type": "const",
			"q": [0.85, 0, 1, 1.4]
		},
		{
			"tag": "Top",
			"type": "const",
			"q": [0.85, 0, 1, 1.4]
		},
		{
			"tag": "Left",
			"type": "const",
			"q": [0.85, 0, 1, 1.4]
		},
		{
			"tag": "Right",
			"type": "const",
			"q": [0.85, 0, 1, 1.4]
		},
		{
			"tag": "Airfoil",
			"type": "wall",
			"q": [0.85, 0, 1, 1.4]
		}
	]
}
+/
struct Config
{
	string meshFile;
	double dt;
	double tEnd;
	string limiter;
	string flux;
	long saveIter;
	long plotIter;
	string forceBoundary;
	string integrator;
	bool limited;
	long order;
	double CFL;
	double[4] ic;
	string[] bTags;
	BoundaryType[] bTypes;
	double[][] bc;
}

Config loadConfig(string conf)
{
	JSONValue jConfig = parseJSON(conf);
	Config config;

	double getDouble(T)(T val)
	{
		if(val.type == JSON_TYPE.INTEGER)
		{
			return val.integer.to!double;
		}
		else if(val.type == JSON_TYPE.FLOAT)
		{
			return val.floating;
		}
		else
		{
			assert(false, "invalid type");
		}
	}

	config.meshFile = jConfig["mesh"].str;
	config.limiter = jConfig["limiter"].str;
	config.flux = jConfig["flux"].str;
	config.dt = jConfig["dt"].floating;
	config.tEnd = getDouble(jConfig["tEnd"]);
	config.saveIter = jConfig["saveIter"].integer;
	config.plotIter = jConfig["plotIter"].integer;
	config.forceBoundary = jConfig["forceBoundary"].str;

	try
	{
		config.CFL = getDouble(jConfig["CFL"]);
	}
	catch(Exception ex)
	{
		writeln("CFL not provided, setting to 0.5");
		config.CFL = 0.5;
	}

	try
	{
		config.integrator = jConfig["integrator"].str;
	}
	catch(Exception ex)
	{
		writeln("Integrator not provided, using forward euler");
		config.integrator = "Euler";
	}

	try
	{
		config.limited = (jConfig["limited"].type == JSON_TYPE.TRUE);
	}
	catch(Exception ex)
	{
		writeln("Limited option not provided, setting to true");
		config.limited = true;
	}

	try
	{
		config.order = jConfig["order"].integer;

		if((config.order != 1) && (config.order != 2))
		{
			writeln("Invalid order supplied, setting order 2");
			config.order = 2;
		}
	}
	catch(Exception ex)
	{
		writeln("Order option not provided, setting order 2");
		config.order = 2;
	}

	auto ics = jConfig["initialConditions"].array;
	config.ic[0] = getDouble(ics[0]);
	config.ic[1] = getDouble(ics[1]);
	config.ic[2] = getDouble(ics[2]);
	config.ic[3] = getDouble(ics[3]);

	auto bcs = jConfig["boudaryConditions"].array;
	for(uint i = 0; i < bcs.length; i++)
	{
		config.bTags ~= bcs[i]["tag"].str;
		immutable string bType = bcs[i]["type"].str;
		if(bType == "fullState")
		{
			config.bTypes ~= BoundaryType.FullState;
		}
		else if(bType == "inviscidWall")
		{
			config.bTypes ~= BoundaryType.InviscidWall;
		}
		else if(bType == "constP")
		{
			config.bTypes ~= BoundaryType.ConstPressure;
		}

		auto state = bcs[i]["q"].array;
		config.bc ~= [getDouble(state[0]), getDouble(state[1]), getDouble(state[2]), getDouble(state[3])];
	}

	return config;
}

void startComputation(Config config, string saveFile)
{
	try
	{
		UMesh2 umesh;

		double dt = config.dt;
		double t = 0;

		auto ex = new SolverException("No error");

		if(config.meshFile.canFind(".gri"))
		{
			umesh = parseXflowMesh(config.meshFile);
		}
		else
		{
			writeln("Unsupported mesh format, exiting");
			return;
		}

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
											runIntegrator!(ufvmSetup, ufvmSolver!(mixin(lim), mixin(fl), 4), mixin(inte))(umesh, config, saveFile, ex);
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
	catch(SolverException ex)
	{
		writeln("Solver encountered an error: ", ex.msg);
		if(ex.exceptionType == SolverException.SExceptionType.EdgeException)
		{
			writeln("pL = ", ex.eExcept.pL);
			writeln("pR = ", ex.eExcept.pR);
			writeln("Flux = ", ex.eExcept.flux);
			writeln("qL = ", ex.eExcept.qL);
			writeln("qR = ", ex.eExcept.qR);
			writeln("cell L = ", ex.eExcept.cellL);
			writeln("cell R = ", ex.eExcept.cellR);
			writeln("normal = ", ex.eExcept.normal);
		}
		writeln("exiting");
	}
}

void main(string[] args)
{
	import std.getopt;
	string configFile;
	string plotAddr = "127.0.0.1";
	string saveFile = "";
	ushort plotPort = 54000;

	signal(SIGINT, &handle);
	auto res = getopt(args, "c|config", "config file to read", &configFile, "pa|plotAddr", "IP address to plot to", &plotAddr, 
							"pp|plotPort", "Port to plot to", &plotPort, "s|save", "Save file to start from", &saveFile);

	//initRPP(plotAddr, plotPort);

	auto configStr = readText(configFile);
	auto config = loadConfig(configStr);

	startComputation(config, saveFile);
}
