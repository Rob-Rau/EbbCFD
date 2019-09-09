/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.solver;

import core.atomic;
import core.stdc.stdio : fopen, fwrite, fopen, printf, snprintf;
import core.sys.posix.signal;

import std.algorithm : canFind, countUntil, min, map, max, reduce, sum;
import std.conv;
import std.file;
import std.math;
import std.meta;
import std.stdio;
import std.string;

import numd.utility;
import numd.linearalgebra.matrix;

import ebb.config;
import ebb.euler;
import ebb.exception;
import ebb.flux;
import ebb.integrators;
import ebb.limiters;
import ebb.mesh;
import ebb.io;
import ebb.mpid;
import ebb.solve;

static shared bool interrupted = false;

@nogc @system nothrow extern(C) void handle(int sig)
{
	printf("Signal received\n");
	atomicStore(interrupted, true);
}

@nogc uint runIntegrator(alias setup, alias solver, alias integrator)(ref UMesh2 mesh, Config config, string saveFile, uint[] triMap)
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

	FILE* forceFile;
	if(mesh.mpiRank == 0)
	{
		forceFile = fopen("boundaryForces.frc", "wb");
	}

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
	Vector!4[] triQ = cast(Vector!4[])Mallocator.instance.allocate(triMap.length*Vector!4.sizeof);

	lastRho[] = 0.0;
	thisRho[] = 0.0;
	lastU[] = 0.0;
	thisU[] = 0.0;
	lastV[] = 0.0;
	thisV[] = 0.0;
	lastE[] = 0.0;
	thisE[] = 0.0;
	tmp[] = 0.0;

	Vector!4[] q0;
	Vector!4[] q1;
	Vector!4[] q2;
	Vector!4[] dem;

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
	scope(exit) Mallocator.instance.deallocate(triQ);

	// let the integrator do any neccessary initialization
	integrator.init(mesh);

	double lastRmax = 0.0;

	// Setup IC's and BC's
	setup(mesh, config, lastRho, lastU, lastV, lastE, t, dt, saveFile);
	mesh.comm.barrier;
	while((t < config.tEnd) && !atomicLoad(interrupted))
	{
		double Rmax = 0;
		double newDt = dt;

		integrator.step!solver(R, mesh.q, mesh, config, newDt, Rmax);

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
			
			immutable bool haveNan = mesh.q[i][0].isNaN || mesh.q[i][1].isNaN || mesh.q[i][2].isNaN || mesh.q[i][3].isNaN;
			if(haveNan)
			{
				printf("p = %f\n", getPressure(mesh.q[i]));
				printf("q = [%f, %f, %f, %f]\n", mesh.q[i][0], mesh.q[i][1], mesh.q[i][2], mesh.q[i][3]);
				printf("cell = %d\n", i);
				printf("iteration = %d\n", iterations);
				printf("dt = %f\n", dt);
				printf("t = %f\n", t);
				enforce!CellException(!haveNan, "Got nan on cell average value", mesh.cells[i]);
			}
		}

		import std.algorithm : sum;

		@nogc double computeRMSResidual(double[] now, double[] last)
		{
			tmp[] = (now[] - last[])^^2;
			double tmpSum = tmp.sum;
			mesh.comm.allreduce(tmpSum, tmpSum, MPI_SUM);

			uint totLen = cast(uint)mesh.interiorCells.length;
			mesh.comm.allreduce(totLen, totLen, MPI_SUM);

			double valSum = now.sum;
			mesh.comm.allreduce(valSum, valSum, MPI_SUM);
			
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
		if(config.cflAdjust)
		{
			if(residRhoLast < residMax)
			{
				residRhoIncIters++;
			}
			else
			{
				if(residRhoIncIters >= 1)
				{
					residRhoIncIters--;
				}
			}

			if(residRhoIncIters >= 2000)
			{
				config.CFL *= 0.99;
				printf("Max RMS residual hasn't decreased in 2000 iterations, decreasing CFL to %f\n", config.CFL);
				residRhoIncIters = 0;
			}
		}

		residRhoLast = residMax;

		if(config.plotIter != -1)
		{
			if(iterations % config.plotIter == 0)
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
			if(iterations % config.saveIter == 0)
			{
				char[512] filename;
				filename[] = 0;
				snprintf(filename.ptr, 512, "save_%d_%d.esln", saveItr, mesh.mpiRank);
				foreach(i, cell; triMap)
				{
					triQ[i] = mesh.q[cell];
				}
				saveSolution(mesh.q, mesh, filename.ptr, t, dt, cast(uint)config.order);
				saveItr++;
			}
		}
		t += dt;
		iterations++;
		lastRmax = Rmax;
	}

	if(mesh.mpiRank == 0)
	{
		if(buffPos != 0)
		{
			fwrite(forceBuffer.ptr, ubyte.sizeof, buffPos, forceFile);
			buffPos = 0;
		}
		fclose(forceFile);
	}

	char[512] filename;
	filename[] = 0;
	snprintf(filename.ptr, 512, "final_%d.esln", mesh.mpiRank);
	foreach(i, cell; triMap)
	{
		triQ[i] = mesh.q[cell];
	}
	saveSolution(mesh.q, mesh, filename.ptr, t, dt, cast(uint)config.order);
	if(mesh.mpiRank == 0)
	{
		printf("lift force = %f\t drag force = %f\t t = %f\n", ld[1], ld[0], t);
		printf("rho_RMS = %.10e\tu_RMS = %.10e\tv_RMS = %.10e\tE_RMS = %.10e\tFlux_R = %.10e\t dt = %10.10f\n", residRho, residU, residV, residE, lastRmax, dt);
	}

	return iterations;
}

string switchBuilder(int level, string switchVar, Args...)(string statement)
{
	alias list = AliasSeq!Args;

	string fillPlaceHolder(int level)(in string statement, string arg)
	{
		auto strSlice = statement;
		ptrdiff_t searchIdx = 0;
		string newStatement = statement;
		while(searchIdx < statement.length)
		{
			auto idxStart = newStatement.indexOf('{', searchIdx);
			if(idxStart == -1)
			{
				break;
			}
			
			auto idxEnd = newStatement.indexOf('}', idxStart);
			assert(idxEnd > idxStart);
			
			auto sliceLen = idxEnd - (idxStart + 1);
			if(newStatement[idxStart+1..idxStart+1 + sliceLen].isNumeric)
			{
				auto strLevel = newStatement[idxStart+1..idxStart+1 + sliceLen].to!int;
				if(strLevel == level)
				{
					newStatement = newStatement[0..idxStart] ~ arg ~ newStatement[idxEnd + 1..$];
				}
			}

			searchIdx = idxStart + 1;
		}
		
		return newStatement;
	}

	string switchStatement = `final switch(`~switchVar~`)
{`;
	foreach(arg; list)
	{
		auto thisStatement = fillPlaceHolder!level(statement, arg.stringof);
		switchStatement ~= `
case `~arg.stringof~`:
	`~thisStatement~`
	break;`;
	}
	switchStatement ~= `
}`;
	return switchStatement;
}

void startComputation(Config config, string saveFile, uint p, uint id)
{
	try
	{
		auto umesh = UMesh2(MPI_COMM_SELF, id);

		double dt = config.dt;
		double t = 0;

		if(id == 0)
		{
			umesh = importMesh(config.meshFile);
			umesh.comm = MPI_COMM_SELF;
			umesh.mpiRank = id;
		}

		version(Have_mpi)
		{
			if(p > 1)
			{
				umesh = partitionMesh(umesh, p, id, MPI_COMM_WORLD);
				umesh.comm = MPI_COMM_WORLD;
				umesh.mpiRank = id;
			}
			else
			{
				import std.array : array;
				import std.range : iota;
				umesh.localToGlobalElementMap = iota(0, umesh.elements.length).array.to!(uint[]);
			}
		}

		umesh.buildMesh;

		import std.array : split;
		auto triMesh = umesh.triangulate;
		auto triMap = triMesh[1];
		auto tMesh = triMesh[0];

		umesh.comm.barrier;
		tMesh.buildMesh;
		umesh.comm.barrier;
		saveMatlabMesh(tMesh, config.meshFile.split('.')[0]~"_"~id.to!string~".mmsh");
		umesh.comm.barrier;

		mixin(`writeln("Running 2D finite volume solver");
				writeln("-limited: ", config.limited);
				writeln("-order: ", config.order);
				writeln("-flux: "~{1});
				writeln("-integrator: "~{0});
				double startTime = wtime;
				auto iterations = runIntegrator!(ufvmSetup, ufvmSolver!(mixin({1}), 4), mixin({0}))(umesh, config, saveFile, triMap);
				MPI_COMM_WORLD.barrier;
				double elapsed = wtime - startTime;
				if(id == 0)
				{
					writeln("total solver time: ", elapsed, " s");
					auto timePerIteration = elapsed/cast(double)iterations;
					writeln("average time per iteration: ", timePerIteration, " s");
					auto timePerCell = timePerIteration/cast(double)umesh.q.length;
					writeln("average time per cell: ", timePerCell*1_000_000.0, " us");
				}`.
				switchBuilder!(1, "config.flux", fluxList).
				switchBuilder!(0, "config.integrator", integratorList));
	}
	catch(CellException ce)
	{
		writeln("Solver encountered an error: ", ce.msg);
		writeln("	In file ", ce.file);
		writeln("	On line ", ce.line);
		MPI_COMM_WORLD.abort(1);
	}
	catch(EdgeException ex)
	{
		writeln("Solver encountered an error: ", ex.msg);
		writeln("	In file ", ex.file);
		writeln("	On line ", ex.line);
		writeln("	pL = ", getPressure(ex.edge.q[0]));
		writeln("	pR = ", getPressure(ex.edge.q[1]));
		writeln("	Flux = ", ex.edge.flux);
		writeln("	qL = ", ex.edge.q[0]);
		writeln("	qR = ", ex.edge.q[1]);
		writeln("	cell L = ", ex.edge.cellIdx[0]);
		writeln("	cell R = ", ex.edge.cellIdx[1]);
		writeln("	normal = ", ex.edge.normal);
		MPI_COMM_WORLD.abort(1);
	}
	catch(Exception ex)
	{
		writeln("Caught unknown exception. Message: ", ex.msg);
		writeln("	In file ", ex.file);
		writeln("	On line ", ex.line);
		MPI_COMM_WORLD.abort(1);
	}
	finally
	{
		writeln("exiting");
	}
}

int main(string[] args)
{
	import std.getopt : getopt;

	string configFile;
	string saveFile = "";

	init(args);

	int p = MPI_COMM_WORLD.size;
	int id = MPI_COMM_WORLD.rank;

	vec4dataType = toMPIType!(Vector!4);

	double startTime = wtime;

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

	MPI_COMM_WORLD.barrier;
	double elapsed = wtime - startTime;
	if(id == 0)
	{
		writeln("total run time: ", elapsed);
	}

	return shutdown;
}
