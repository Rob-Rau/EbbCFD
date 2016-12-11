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
import numd.optimization.ArrayOps;
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

		AbstractMeshOpt meshOpt;
		auto sqp = new SQP;
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
											meshOpt = new MeshOpt!(ufvmSetup, ufvmSolver!(mixin(lim), mixin(fl), 4), mixin(inte))(config, p, id, sqp);
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

		meshOpt.DerivativeType = "numd.optimization.FiniteDifference.FiniteDifferenceEqualized";
		//meshOpt.StepSize = 5.0e-3;
		double stepSize = 300.0/meshOpt.bigMesh.cells.length.to!double;
		MPI_Bcast(&stepSize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		meshOpt.StepSize = stepSize;
		meshOpt.runIterations = 270;

		sqp.DebugMode = true;
		sqp.InitialGuess = new double[p];
		sqp.InitialGuess[] = 1.0/(cast(double)p);
		sqp.PointFilename = "SQPpoints.csv";
		sqp.ErrorFilename = "SQPerror.csv";
		sqp.FileOutput = false;
		sqp.Tolerance = 2.0e-3;
		sqp.id = id;
		sqp.Eta = 0.05;
		logln("Starting optimization");
		MPI_Barrier(MPI_COMM_WORLD);

		while((meshOpt.t < config.tEnd) && !atomicLoad(interrupted))
		{
			Result result = sqp.Optimize(meshOpt);

			if(id == 0)
			{
				writeln();
				writefln("\tOptimal weights =  [%(%20.20f, %)]", result.DesignVariables);
				writefln("\tAverage iteration time = %20.20f", result.ObjectiveFunctionValue);
				writeln("\tConverged in ", result.Iterations, " iterations.");
				writeln("\tComputation time: ", result.ComputationTime/(1000.0*1000.0), " secs.\n");
				writeln("Mesh partitions optimized, running and monitoring simulation");
			}

			//UMesh2 partitionMesh(ref UMesh2 bigMesh, uint p, uint id, MPI_Comm comm, double[] partWeights)
			auto mesh = partitionMesh(meshOpt.bigMesh, p, id, MPI_COMM_WORLD, meshOpt.bestWeights);
			mesh.comm = MPI_COMM_WORLD;
			mesh.mpiRank = id;
			mesh.buildMesh;
			
			import std.experimental.allocator.mallocator : Mallocator;
			meshOpt.lastRho = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
			meshOpt.thisRho = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
			meshOpt.lastU = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
			meshOpt.thisU = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
			meshOpt.lastV = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
			meshOpt.thisV = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
			meshOpt.lastE = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
			meshOpt.thisE = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
			meshOpt.tmp = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
			meshOpt.R = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);

			// Scatter bigMesh to smaller meshes
			ufvmSetup(mesh, config, meshOpt.lastRho, meshOpt.lastU, meshOpt.lastV, meshOpt.lastE, meshOpt.t, meshOpt.dt, "", meshOpt.ex);
			if(id == 0)
			{
				foreach(int proc, localToGlobalMap; meshOpt.bigMesh.localToGlobalMaps)
				{
					if(proc != 0)
					{
						//logln("Sending to proc ", proc);
						auto localState = new Vector!4[localToGlobalMap.length];

						foreach(i, globalIdx; localToGlobalMap)
						{
							localState[i] = meshOpt.bigMesh.q[globalIdx];
						} 
						MPI_COMM_WORLD.sendArray(localState, proc, mesh.meshTag);
					}
					else
					{
						foreach(j, i; mesh.interiorCells)
						{
							mesh.q[i] = meshOpt.bigMesh.q[localToGlobalMap[j]];
						}
					}
				}
			}
			else
			{
				auto localState = MPI_COMM_WORLD.recvArray!(Vector!4)(0, mesh.meshTag);

				foreach(j, i; mesh.interiorCells)
				{
					mesh.q[i] = localState[j];
				}
			}

			meshOpt.reinitIntegrator(mesh);

			uint iterations = 0;
			double startTime = MPI_Wtime;
			double elapsed;
			while((meshOpt.t < config.tEnd) && !atomicLoad(interrupted))
			{
				meshOpt.solverIteration(mesh);
				iterations++;
				
				if((iterations % meshOpt.runIterations) == 0)
				{
					elapsed = MPI_Wtime - startTime;
					MPI_Bcast(&elapsed, 1, MPI_DOUBLE, 0, mesh.comm);

					if(id == 0) logln("weights: ", meshOpt.bestWeights, ";  Average solver iteration time: ", elapsed/meshOpt.runIterations.to!double);
					if(elapsed > 1.10*meshOpt.minTime)
					{
						if(id == 0) logln("weights: ", meshOpt.bestWeights, ";  elapsed time jumped up by 10%, restarting optimization");
						break;
					}
					startTime = MPI_Wtime;
				}
			}

			// Gather results and put it in bigMesh
			if(id != 0)
			{
				auto localState = new Vector!4[mesh.interiorCells.length];
				foreach(i; mesh.interiorCells)
				{
					localState[i] = mesh.q[i];
				}

				MPI_COMM_WORLD.sendArray(localState, 0, mesh.meshTag);
			}
			else
			{
				foreach(int proc, localToGlobalMap; meshOpt.bigMesh.localToGlobalMaps)
				{
					if(proc != 0)
					{
						auto localState = MPI_COMM_WORLD.recvArray!(Vector!4)(proc, mesh.meshTag);
						foreach(localIdx, globalIdx; localToGlobalMap)
						{
							meshOpt.bigMesh.q[globalIdx] = localState[localIdx];
						}
					}
					else
					{
						foreach(localIdx, globalIdx; localToGlobalMap)
						{
							meshOpt.bigMesh.q[globalIdx] = mesh.q[localIdx];
						}
					}
				}
			}


			Mallocator.instance.deallocate(cast(void[])meshOpt.lastRho);
			Mallocator.instance.deallocate(cast(void[])meshOpt.thisRho);
			Mallocator.instance.deallocate(cast(void[])meshOpt.lastU);
			Mallocator.instance.deallocate(cast(void[])meshOpt.thisU);
			Mallocator.instance.deallocate(cast(void[])meshOpt.lastV);
			Mallocator.instance.deallocate(cast(void[])meshOpt.thisV);
			Mallocator.instance.deallocate(cast(void[])meshOpt.lastE);
			Mallocator.instance.deallocate(cast(void[])meshOpt.thisE);
			Mallocator.instance.deallocate(cast(void[])meshOpt.tmp);
			Mallocator.instance.deallocate(cast(void[])meshOpt.R);

		}
		if(id == 0)
		{
			saveSolution(meshOpt.bigMesh, cast(char*)"final.esln".ptr, meshOpt.t, meshOpt.dt);
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

abstract class AbstractMeshOpt : ObjectiveFunction
{
	uint runIterations;
	UMesh2 bigMesh;
	double t = 0;
	double dt;
	void solverIteration(ref UMesh2 mesh);
	void reinitIntegrator(ref UMesh2 mesh);
	double minTime = double.infinity;
	double[] bestWeights;

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
	SolverException ex;
}

class MeshOpt(alias setup, alias solver, alias integrator) : AbstractMeshOpt 
{
	import std.experimental.allocator.mallocator : Mallocator;
	import std.bitmanip : write;

	int mpiRank;

	double residRho = 0.0;
	double residU = 0.0;
	double residV = 0.0;
	double residE = 0.0;
	double residMax = 0.0;
	double lastRmax = 0.0;

	double residRhoLast;

	Matrix!(2, 2) rotMat;
	Vector!2 ld;

	Config config;

	FILE* forceFile;
	uint totalIterations;
	uint saveItr;
	int p;

	immutable uint buffSize = 3*1024*1024*double.sizeof;
	size_t buffPos = 0;

	SQP sqp;
	this(Config config, int p, int rank, SQP sqp)
	{
		this.sqp = sqp;
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
				bigMesh.buildMesh;
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
		bestWeights = new double[p];
		

		lastRho = cast(double[])Mallocator.instance.allocate(bigMesh.cells.length*double.sizeof);
		thisRho = cast(double[])Mallocator.instance.allocate(bigMesh.cells.length*double.sizeof);
		lastU = cast(double[])Mallocator.instance.allocate(bigMesh.cells.length*double.sizeof);
		thisU = cast(double[])Mallocator.instance.allocate(bigMesh.cells.length*double.sizeof);
		lastV = cast(double[])Mallocator.instance.allocate(bigMesh.cells.length*double.sizeof);
		thisV = cast(double[])Mallocator.instance.allocate(bigMesh.cells.length*double.sizeof);
		lastE = cast(double[])Mallocator.instance.allocate(bigMesh.cells.length*double.sizeof);
		thisE = cast(double[])Mallocator.instance.allocate(bigMesh.cells.length*double.sizeof);
		tmp = cast(double[])Mallocator.instance.allocate(bigMesh.cells.length*double.sizeof);
		R = cast(Vector!4[])Mallocator.instance.allocate(bigMesh.cells.length*Vector!4.sizeof);
		forceBuffer = cast(ubyte[])Mallocator.instance.allocate(buffSize);
		
		scope(exit) Mallocator.instance.deallocate(cast(void[])lastRho);
		scope(exit) Mallocator.instance.deallocate(cast(void[])thisRho);
		scope(exit) Mallocator.instance.deallocate(cast(void[])lastU);
		scope(exit) Mallocator.instance.deallocate(cast(void[])thisU);
		scope(exit) Mallocator.instance.deallocate(cast(void[])lastV);
		scope(exit) Mallocator.instance.deallocate(cast(void[])thisV);
		scope(exit) Mallocator.instance.deallocate(cast(void[])lastE);
		scope(exit) Mallocator.instance.deallocate(cast(void[])thisE);
		scope(exit) Mallocator.instance.deallocate(cast(void[])tmp);
		scope(exit) Mallocator.instance.deallocate(cast(void[])R);
		
		setup(bigMesh, config, lastRho, lastU, lastV, lastE, t, dt, "", ex);
	}

	~this()
	{
		scope(exit) Mallocator.instance.deallocate(cast(void[])forceBuffer);
	}
	
	final override Complex!double Compute(Complex!double[] designVar)
	{
		return doCompute(designVar, 0);
	}

	override void reinitIntegrator(ref UMesh2 mesh)
	{
		integrator.reinit(mesh);
	}

	override void solverIteration(ref UMesh2 mesh)
	{
		double startTime2 = MPI_Wtime;
		double Rmax = 0;
		double newDt = dt;

		//if(mpiRank == 0) logln(iterations, ": Stepping");
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
				//printf("iteration = %d\n", iterations);
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

		if(config.plotIter >= 0)
		{
			if(totalIterations % config.plotIter == 0)
			{
				if(mesh.mpiRank == 0)
				{
					printf("lift force = %f\t drag force = %f\t t = %f\n", ld[1], ld[0], t);
					printf("rho_RMS = %.10e\tu_RMS = %.10e\tv_RMS = %.10e\tE_RMS = %.10e\tFlux_R = %.10e\t dt = %10.10f\n", residRho, residU, residV, residE, Rmax, dt);
				}
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
		totalIterations++;
		lastRmax = Rmax;
	}

	Complex!double doCompute(Complex!double[] designVar, int depth)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		sqp.stop = atomicLoad(interrupted);
		//logln("in compute");
		UMesh2 mesh;

		//logln("paritioning mesh");
		double[] w = new double[p];
		/+
		w[] = 1.0/(cast(double)p);

		double equalTime = 1.0;
		if(depth == 0)21770105
		{
			equalTime = doCompute(w.complex, 1).re;
		}
		+/

		
		auto down = 1.0/p.to!double * 0.85;
		//w[] = (1.0 - 2.0*down)/((p - 2).to!double);
		//w[1] = down;
		//w[3] = down;
		//w[] = 1.0/p.to!double;
		if(mpiRank == 0)
		{
			//if(depth == 0) logln("w = ", w);
		}
		
		//logln("partSum = ", designVar.sum!".re");
		mesh = partitionMesh(bigMesh, p, mpiRank, MPI_COMM_WORLD, w);
		MPI_Barrier(MPI_COMM_WORLD);
		mesh.comm = MPI_COMM_WORLD;
		mesh.mpiRank = mpiRank;

		double equalTime = 1.0;
		
		//logln("building mesh");
		mesh.buildMesh;
		MPI_Barrier(MPI_COMM_WORLD);
		//logln("mesh cells = ", mesh.interiorCells.length);
		size_t[] meshSizes = new size_t[p];
		if(mpiRank == 0)
		{
			meshSizes[0] = mesh.interiorCells.length;
			for(int i = 1; i < p; i++)
			{
				meshSizes[i] = mesh.comm.recv!size_t(i, mesh.meshTag);
			}
			//logln(meshSizes);
		}
		else
		{
			mesh.comm.send(mesh.interiorCells.length, 0, mesh.meshTag);
		}

		if(mesh.interiorCells.length != bigMesh.interiorCells.length)
		{
			w[] = 1.0/(cast(double)p);

			if(depth == 0)
			{
				//equalTime = doCompute(w.complex, 1).re;
			}
		}

		
		foreach(uint i, dv; designVar)
		{
			w[i] = dv.re;
		}
		// let the integrator do any neccessary initialization
		//logln("Re initing integrator");
		integrator.reinit(mesh);

		lastRho = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		thisRho = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		lastU = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		thisU = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		lastV = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		thisV = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		lastE = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		thisE = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		tmp = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
		
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

		// TODO: Redo this so that it doesn't always start the sim again
		setup(mesh, config, lastRho, lastU, lastV, lastE, t, dt, "", ex);
		if(mpiRank == 0)
		{
			foreach(int proc, localToGlobalMap; bigMesh.localToGlobalMaps)
			{
				if(proc != 0)
				{
					//logln("Sending to proc ", proc);
					auto localState = new Vector!4[localToGlobalMap.length];

					foreach(i, globalIdx; localToGlobalMap)
					{
						localState[i] = bigMesh.q[globalIdx];
					} 
					MPI_COMM_WORLD.sendArray(localState, proc, mesh.meshTag);
				}
				else
				{
					foreach(j, i; mesh.interiorCells)
					{
						mesh.q[i] = bigMesh.q[localToGlobalMap[j]];
					}
				}
			}
		}
		else
		{
			auto localState = MPI_COMM_WORLD.recvArray!(Vector!4)(0, mesh.meshTag);

			foreach(j, i; mesh.interiorCells)
			{
				mesh.q[i] = localState[j];
			}
		}

		//logln("Going");
		import core.memory: GC;

		GC.disable;

		MPI_Barrier(mesh.comm);
		double startTime = MPI_Wtime;
		double elapsed;
		uint iterations;

		if(mesh.interiorCells.length != bigMesh.interiorCells.length)
		{
			while((iterations < runIterations))
			{
				solverIteration(mesh);
				iterations++;
				//if(mpiRank == 0) logln("current elapsed: ", MPI_Wtime - startTime2);
			}
			MPI_Barrier(mesh.comm);
			elapsed = MPI_Wtime - startTime;

			GC.enable;
			GC.collect;
			if(mpiRank != 0)
			{
				auto localState = new Vector!4[mesh.interiorCells.length];
				foreach(i; mesh.interiorCells)
				{
					localState[i] = mesh.q[i];
				}

				MPI_COMM_WORLD.sendArray(localState, 0, mesh.meshTag);
			}
			else
			{
				foreach(int proc, localToGlobalMap; bigMesh.localToGlobalMaps)
				{
					if(proc != 0)
					{
						auto localState = MPI_COMM_WORLD.recvArray!(Vector!4)(proc, mesh.meshTag);
						foreach(localIdx, globalIdx; localToGlobalMap)
						{
							bigMesh.q[globalIdx] = localState[localIdx];
						}
					}
					else
					{
						foreach(localIdx, globalIdx; localToGlobalMap)
						{
							bigMesh.q[globalIdx] = mesh.q[localIdx];
						}
					}
				}
			}
		}
		else
		{
			MPI_Barrier(mesh.comm);
			elapsed = 100*equalTime*p;
			GC.enable;
			GC.collect;
		}


		//logln(runIterations, " iterations took ", elapsed, " seconds");
		MPI_Barrier(mesh.comm);
		MPI_Bcast(&elapsed, 1, MPI_DOUBLE, 0, mesh.comm);

		if(((elapsed/runIterations.to!double) < minTime) && (depth == 0))
		{
			minTime = elapsed/runIterations.to!double;
			bestWeights[] = w[];
			//if(mpiRank == 0) logln("new minimum time: ", minTime, " seconds");
			//if(mpiRank == 0) logln("	optimal weights: ", bestWeights);
		}

		if(mpiRank == 0)
		{
			//if(depth == 0) logln(runIterations, " iterations took ", elapsed, " seconds");
			//if(depth == 0) logln(runIterations, " iterations took ", elapsed/equalTime, " normalized seconds; equal partition time: ", equalTime, " seconds; w = ", w);
			if(depth == 0) logln("mesh sizes: ", meshSizes, ";   ", runIterations, " iterations took ", (elapsed), " seconds; average iteration time: ", elapsed/runIterations.to!double);
		}

		MPI_Barrier(mesh.comm);

		if(depth == 1)
		{
			return complex(elapsed, 0.0);
		}
		else
		{
			//return complex((elapsed/runIterations.to!double)/(equalTime/runIterations.to!double), 0.0);
			return complex(elapsed/runIterations.to!double, 0.0);
		}
	}

	alias double delegate(double x) ActiveConstraint;

	import std.container.dlist : DList;
	import std.typecons : tuple, Tuple;

	DList!(Tuple!(int, ActiveConstraint, string)) constraints;

	final override void UpdateActiveSet(double[] designVars)
	{
		import std.algorithm : find;
		foreach(int i, designVar; designVars)
		{
			if((designVar >= 0.0) && (designVar <= 1.0))
			{
				import std.range : array;
				// remove designVar constraint if it exists
				auto c = constraints[].find!("a[0] == b")(i);
				if(!c.empty)
				{
					constraints.remove(c);
				}
			}
			else if(designVar > 1.0)
			{
				// Add 0 = 1 - designVar constraint
				auto foundConstraints = constraints[].find!("a[0] == b")(i);
				if(foundConstraints.empty)
				{
					ActiveConstraint ac = (double x) => 1.0 - x;
					constraints.stableInsertBack(tuple(i, ac, "g1"));
				}
				else
				{
					foundConstraints = constraints[].find!("((a[0] == b) && (a[2] != \"g1\"))")(i);
					if(!foundConstraints.empty)
					{
						constraints.remove(foundConstraints);
						ActiveConstraint ac = (double x) => 1.0 - x;
						constraints.stableInsertBack(tuple(i, ac, "g1"));
					}
				}
			}
			else if(designVar < 0.0)
			{
				// Add 0 = designVar constraint
				auto foundConstraints = constraints[].find!("a[0] == b")(i);
				if(foundConstraints.empty)
				{
					ActiveConstraint ac = (double x) => x;
					constraints.stableInsertBack(tuple(i, ac, "l0"));
				}
				else
				{
					foundConstraints = constraints[].find!("((a[0] == b) && (a[2] != \"l0\"))")(i);
					if(!foundConstraints.empty)
					{
						constraints.remove(foundConstraints);
						ActiveConstraint ac = (double x) => x;
						constraints.stableInsertBack(tuple(i, ac, "l0"));
					}
				}
			}
			else
			{
				if(mpiRank == 0) writeln("Unkown condition, designVar = ", designVar);
			}
		}
	}

	final override Complex!double[] Constraint(Complex!double[] designVar)
	{
		Complex!double[] c = new Complex!double[1];

		double partSum = designVar.sum!".re";
		c[0].re = 1 - partSum;
		c[0].im = 0;

		foreach(constraint; constraints)
		{
			Complex!double co;
			co.re = constraint[1](designVar[constraint[0]].re);
			co.im = 0.0;
			c ~= co;
		}
		//if((mpiRank == 0) && (c.length > 1)) writeln;
		Constraints = cast(int)c.length;
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
