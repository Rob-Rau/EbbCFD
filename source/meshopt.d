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
	this(Config config)
	{
		Constraints = 1;
	}

	final override Complex!double Compute(Complex!double[] designVar)
	{


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
	return 0;
}