module ebb.meshopt;

import core.atomic;
import core.stdc.stdio : fopen, fwrite, fopen, printf, snprintf;
import core.sys.posix.signal;

import std.algorithm : canFind, countUntil, min, max, reduce;
import std.complex;
import std.conv;
import std.file;
import std.math;
import std.stdio;

import numd.utility;
import numd.linearalgebra.matrix;
import numd.optimization.ObjectiveFunction;
import numd.optimization.Optimizer;
import numd.optimization.Gradient;
import numd.optimization.Derivative;
import numd.optimization.ComplexStep;
import numd.optimization.FiniteDifference;
import numd.optimization.SQP;

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

void startOptimization(Config config, string saveFile, uint p, uint id)
{
	try
	{
		auto umesh = UMesh2(MPI_COMM_SELF, id);

		double dt = config.dt;
		double t = 0;

		auto ex = new SolverException("No error");
/+
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
+/
		//umesh.buildMesh;

		//import std.array : split;
		//saveMatlabMesh(umesh, config.meshFile.split('.')[0]~"_"~id.to!string~".mmsh");

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
											meshOpt = new MeshOpt!(ufvmSetup, ufvmSolver!(mixin(lim), mixin(fl), 4), mixin(inte))(config, p, id);
											MPI_Barrier(MPI_COMM_WORLD);
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

		auto sqp = new SQP;

		meshOpt.DerivativeType = "numd.optimization.FiniteDifference.FiniteDifference";
		meshOpt.StepSize = 1.0e-2;

		sqp.InitialGuess = new double[p];
		sqp.InitialGuess[] = 1.0/(cast(double)p);
		sqp.PointFilename = "SQPpoints.csv";
		sqp.ErrorFilename = "SQPerror.csv";
		sqp.FileOutput = false;
		logln("Starting optimization");
		MPI_Barrier(MPI_COMM_WORLD);
		Result result = sqp.Optimize(meshOpt);
		writeln();
		writeln("SQP:");
		writefln("\tOptimal Point =  [%(%20.20f, %)]", result.DesignVariables);
		writefln("\tDrag = %20.20f", result.ObjectiveFunctionValue);
		writeln("\tConverged in ", result.Iterations, " iterations.");
		writeln("\tComputation time: ", result.ComputationTime, " usecs.");
		writeln("\tMinor iterations: ", result.MinorIterations);
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

	double t;
	double dt;

	double residRhoLast;

	Matrix!(2, 2) rotMat;
	Vector!2 ld;

	Config config;

	SolverException ex;

	FILE* forceFile;
	uint totalIterations;
	uint saveItr;
	int p;

	this(Config config, int p, int rank)
	{
		this.p = p;
		mpiRank = rank;

		this.config = config;
		dt = config.dt;
		
		ex = new SolverException("No error");

		Constraints = 1;
		if(mpiRank == 0)
		{
			if(config.meshFile.canFind(".gri"))
			{
				bigMesh = parseXflowMesh(config.meshFile);
				bigMesh.comm = MPI_COMM_SELF;
				bigMesh.mpiRank = mpiRank;
				logln("big mesh build");
				//bigMesh.buildMesh;
			}
			else
			{
				logln("Unsupported mesh format, exiting");
				return;
			}
		}

		residRhoLast = double.infinity;

		uint residRhoIncIters = 0;

		uint iterations = 0;
		uint saveItr = 0;
		double t = 0;
		double dt = config.dt;

		rotMat = Matrix!(2, 2)(cos(config.ic[1] * (PI/180)), -sin(config.ic[1] * (PI/180)), sin(config.ic[1] * (PI/180)), cos(config.ic[1] * (PI/180)));
		ld = Vector!2(0);

		if(mpiRank == 0)
		{
			forceFile = fopen("boundaryForces.frc", "wb");
		}

		totalIterations = 0;
		saveItr = 0;
		// Setup IC's and BC's
		//setup(mesh, config, lastRho, lastU, lastV, lastE, t, dt, saveFile, ex);
		//MPI_Barrier(mesh.comm);
		bool done = false;
	}

	final override Complex!double Compute(Complex!double[] designVar)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		logln("in compute");
		UMesh2 mesh;

		logln("paritioning mesh");
		mesh = partitionMesh(bigMesh, p, mpiRank, MPI_COMM_WORLD);
		mesh.comm = MPI_COMM_WORLD;
		mesh.mpiRank = mpiRank;

		logln("building mesh");
		mesh.buildMesh;
		// let the integrator do any neccessary initialization
		logln("Re initing integrator");
		integrator.reinit(mesh);

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

		MPI_Barrier(mesh.comm);
		double startTime = MPI_Wtime;
		uint iterations;

		while((iterations < 10) && !atomicLoad(interrupted))
		{
			double Rmax = 0;
			double newDt = dt;

			integrator.step!solver(R, mesh.q, mesh, config, newDt, Rmax, ex);

			if(config.dynamicDt)
			{
				dt = newDt;
			}

			foreach(i; mesh.interiorCells)
			{
				thisRho[i] = mesh.q[i][0];
				thisU[i] = mesh.q[i][1];
				thisV[i] = mesh.q[i][2];
				thisE[i] = mesh.q[i][3];
				
				if(mesh.q[i][0].isNaN || mesh.q[i][1].isNaN || mesh.q[i][2].isNaN || mesh.q[i][3].isNaN)
				{
					printf("p = %f\n", getPressure(mesh.q[i]));
					printf("q = [%f, %f, %f, %f]\n", mesh.q[i][0], mesh.q[i][1], mesh.q[i][2], mesh.q[i][3]);
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

			@nogc double computeRMSResidual(double[] now, double[] last)
			{
				tmp[] = (now[] - last[])^^2;
				double tmpSum = tmp.sum;
				MPI_Allreduce(&tmpSum, &tmpSum, 1, MPI_DOUBLE, MPI_SUM, mesh.comm);
				uint totLen = cast(uint)mesh.interiorCells.length;
				MPI_Allreduce(&totLen, &totLen, 1, MPI_UINT32_T, MPI_SUM, mesh.comm);
				double valSum = now.sum;
				MPI_Allreduce(&valSum, &valSum, 1, MPI_DOUBLE, MPI_SUM, mesh.comm);
				return sqrt(cast(double)totLen*tmpSum)/valSum;
			}

			residRho = computeRMSResidual(thisRho, lastRho);
			residU = computeRMSResidual(thisU, lastU);
			residV = computeRMSResidual(thisV, lastV);
			residE = computeRMSResidual(thisE, lastE);

			lastRho[] = thisRho[];
			lastU[] = thisU[];
			lastV[] = thisV[];
			lastE[] = thisE[];

			auto f = mesh.computeBoundaryForces(config.forceBoundary);

			if(mesh.mpiRank == 0)
			{
				ld = rotMat*f;

				forceBuffer.write!double(t, &buffPos);
				forceBuffer.write!double(ld[0], &buffPos);
				forceBuffer.write!double(ld[1], &buffPos);
			}

			if(buffPos == buffSize)
			{
				fwrite(forceBuffer.ptr, ubyte.sizeof, buffSize, forceFile);
				buffPos = 0;
			}

			residMax = max(residRho, residU, residV, residE);

			if(totalIterations % config.plotIter == 0)
			{
				if(mesh.mpiRank == 0)
				{
					printf("lift force = %f\t drag force = %f\t t = %f\n", ld[1], ld[0], t);
					printf("rho_RMS = %.10e\tu_RMS = %.10e\tv_RMS = %.10e\tE_RMS = %.10e\tFlux_R = %.10e\t dt = %10.10f\n", residRho, residU, residV, residE, Rmax, dt);
				}
			}
			
			if(config.saveIter != -1)
			{
				if(totalIterations % config.saveIter == 0)
				{
					char[512] filename;
					filename[] = 0;
					snprintf(filename.ptr, 512, "save_%d_%d.esln", saveItr, mesh.mpiRank);
					saveSolution(mesh, filename.ptr, t, dt);
					filename[] = 0;
					snprintf(filename.ptr, 512, "save_%d_%d.lsln", saveItr, mesh.mpiRank);
					saveLimits(mesh, filename.ptr, t, dt);
					saveItr++;
				}
			}
			t += dt;
			iterations++;
			totalIterations++;
			lastRmax = Rmax;
		}
		double elapsed = MPI_Wtime - startTime;
		writeln("10 iterations took ", elapsed, " seconds");
		return complex(elapsed, 0.0);
	}

	final override Complex!double[] Constraint(Complex!double[] designVar)
	{
		Complex!double[] c = new Complex!double[1];

		double partSum = designVar.sum!".re";
		c[0].re = 1 - partSum;
		c[0].im = 0;
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
	
	string configFile;
	string saveFile = "";

	vec4dataType = toMPIType!(Vector!4);

	double startTime = MPI_Wtime();

	signal(SIGINT, &handle);
	signal(SIGUSR1, &handle);
	auto res = getopt(args, "c|config", "config file to read", &configFile);

	auto configStr = readText(configFile);
	auto config = loadConfig(configStr);

	if(saveFile != "")
	{
		saveFile ~= id.to!string ~ ".esln";
	}

	writeln("Starting opt");
	startOptimization(config, saveFile, p, id);

	MPI_Barrier(MPI_COMM_WORLD);
	double elapsed = MPI_Wtime() - startTime;
	if(id == 0)
	{
		writeln("total time: ", elapsed);
	}
	
	return MPI_Shutdown;
}