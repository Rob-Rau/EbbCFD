module ebb.meshopt;

import core.atomic;
import core.stdc.stdio : fopen, fwrite, fopen, printf, snprintf;
import core.sys.posix.signal;

import std.algorithm : canFind, countUntil, min, max, reduce, sum;
import std.complex;
import std.conv;
import std.file;
import std.math;
import std.stdio;

import numd.utility;
import numd.linearalgebra.matrix;
import numd.optimization.ObjectiveFunction;
import numd.optimization.Gradient;
import numd.optimization.Derivative;
import numd.optimization.ComplexStep;
import numd.optimization.FiniteDifference;

import ebb.config;
import ebb.euler;
import ebb.exception;
import ebb.flux;
import ebb.integrators;
import ebb.limiters;
import ebb.mesh;
import ebb.io;
import ebb.mpid;
import ebb.solver;

void startOptimization(Config config, string saveFile, uint p, int runIterations, uint id)
{
	try
	{
		auto umesh = UMesh2(MPI_COMM_SELF, id);

		double dt = config.dt;
		double t = 0;

		auto ex = new SolverException("No error");

		if(id == 0)
		{
			if(config.meshFile.canFind(".gri"))
			{
				umesh = parseXflowMesh(config.meshFile);
				umesh.comm = MPI_COMM_SELF;
				umesh.mpiRank = id;
			}
			else
			{
				writeln("Unsupported mesh format, exiting");
				return;
			}
		}

		if(p > 1)
		{
			umesh = partitionMesh(umesh, p, id, MPI_COMM_WORLD);
			umesh.comm = MPI_COMM_WORLD;
			umesh.mpiRank = id;
		}

		umesh.buildMesh;

		import std.array : split;
		saveMatlabMesh(umesh, config.meshFile.split('.')[0]~"_"~id.to!string~".mmsh");

		ObjectiveFunction meshOpt;
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
											meshOpt = new MeshOpt!(ufvmSetup, ufvmSolver!(mixin(lim), mixin(fl), 4), mixin(inte))(config);
											//runIntegrator!(ufvmSetup, ufvmSolver!(mixin(lim), mixin(fl), 4), mixin(inte))(umesh, config, saveFile, runIterations, ex);
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
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

class MeshOpt(alias setup, alias solver, alias integrator) : ObjectiveFunction
{
	import std.experimental.allocator.mallocator : Mallocator;
	import std.bitmanip : write;

	UMesh2 bigMesh;
	int mpiRank;

	double residRho = 0.0;
	double residU = 0.0;
	double residV = 0.0;
	double residE = 0.0;
	double residMax = 0.0;

	double[] lastRho;
	double[] thisRho;
	double[] lastU;
	double[] thisU;
	double[] lastV;
	double[] thisV;
	double[] lastE;
	double[] thisE;
	double[] tmp;
	ubyte[] forceBuffer;
	Vector!4[] R;

	double lastRmax = 0.0;

	this(Config config)
	{
		Constraints = 1;
		if(mpiRank == 0)
		{
			if(config.meshFile.canFind(".gri"))
			{
				bigMesh = parseXflowMesh(config.meshFile);
				bigMesh.comm = MPI_COMM_SELF;
				bigMesh.mpiRank = mpiRank;
			}
			else
			{
				writeln("Unsupported mesh format, exiting");
				return;
			}
		}

		double residRhoLast = double.infinity;

		uint residRhoIncIters = 0;

		uint iterations = 0;
		uint saveItr = 0;
		double t = 0;
		double dt = config.dt;

		auto rotMat = Matrix!(2, 2)(cos(config.ic[1] * (PI/180)), -sin(config.ic[1] * (PI/180)), sin(config.ic[1] * (PI/180)), cos(config.ic[1] * (PI/180)));
		auto ld = Vector!2(0);

/+
		FILE* forceFile;
		if(mesh.mpiRank == 0)
		{
			forceFile = fopen("boundaryForces.frc", "wb");
		}
+/

		// Setup IC's and BC's
		//setup(mesh, config, lastRho, lastU, lastV, lastE, t, dt, saveFile, ex);
		//MPI_Barrier(mesh.comm);
		bool done = false;

	}

	final override Complex!double Compute(Complex!double[] designVar)
	{
		UMesh2 mesh;
		// let the integrator do any neccessary initialization
		integrator.init(mesh);

		immutable uint buffSize = 3*1024*1024*double.sizeof;
		size_t buffPos = 0;
		
		lastRho = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		thisRho = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		lastU = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		thisU = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		lastV = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		thisV = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		lastE = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		thisE = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		tmp = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		forceBuffer = cast(ubyte[])Mallocator.instance.allocate(buffSize);
		R = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
		
		scope(exit) Mallocator.instance.deallocate(cast(void[])lastRho);
		scope(exit) Mallocator.instance.deallocate(cast(void[])thisRho);
		scope(exit) Mallocator.instance.deallocate(cast(void[])lastU);
		scope(exit) Mallocator.instance.deallocate(cast(void[])thisU);
		scope(exit) Mallocator.instance.deallocate(cast(void[])lastV);
		scope(exit) Mallocator.instance.deallocate(cast(void[])thisV);
		scope(exit) Mallocator.instance.deallocate(cast(void[])lastE);
		scope(exit) Mallocator.instance.deallocate(cast(void[])thisE);
		scope(exit) Mallocator.instance.deallocate(cast(void[])tmp);
		scope(exit) Mallocator.instance.deallocate(cast(void[])forceBuffer);
		scope(exit) Mallocator.instance.deallocate(cast(void[])R);
		
		
		lastRho[] = 0.0;
		thisRho[] = 0.0;
		lastU[] = 0.0;
		thisU[] = 0.0;
		lastV[] = 0.0;
		thisV[] = 0.0;
		lastE[] = 0.0;
		thisE[] = 0.0;
		tmp[] = 0.0;

		return complex(0.0, 0.0);
	}

	final override Complex!double[] Constraint(Complex!double[] designVar)
	{
		Complex!double[] c = new Complex!double[1];

		return c;
	}
}

import mpi;
import mpi.util;

int main(string[] args)
{
	import std.getopt;
	int argc = cast(int)args.length;
	auto argv = args.toArgv;

	int p;
	int id;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	
/+
	string configFile;
	string saveFile = "";

	

	vec4dataType = toMPIType!(Vector!4);

	double startTime = MPI_Wtime();

	signal(SIGINT, &handle);
	signal(SIGUSR1, &handle);
	auto res = getopt(args, "c|config", "config file to read", &configFile, 
							"s|save", "Save file to start from", &saveFile);

	auto configStr = readText(configFile);
	auto config = loadConfig(configStr);

	if(saveFile != "")
	{
		saveFile ~= id.to!string ~ ".esln";
	}

	startComputation(config, saveFile, p, id);

	MPI_Barrier(MPI_COMM_WORLD);
	double elapsed = MPI_Wtime() - startTime;
	if(id == 0)
	{
		writeln("total time: ", elapsed);
	}
	
	return MPI_Shutdown;
	+/
	return MPI_Shutdown;
}