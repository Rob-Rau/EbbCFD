module ebb.solve;

import core.atomic;
import core.stdc.stdio : fopen, fwrite, fopen, printf, snprintf;
import core.sys.posix.signal;

import std.algorithm : canFind, countUntil, min, max, reduce, sum;
import std.complex;
import std.conv;
import std.file;
import std.math;
import std.stdio;

import ebb.config;
import ebb.mpid;
import ebb.solver;

import numd.linearalgebra.matrix;

import mpi;
import mpi.util;

int main(string[] args)
{
	import std.getopt;

	string configFile;
	string saveFile = "";

	int argc = cast(int)args.length;
	auto argv = args.toArgv;

	int p;
	int id;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

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

	startComputation(config, saveFile, p, -1, id);

	MPI_Barrier(MPI_COMM_WORLD);
	double elapsed = MPI_Wtime() - startTime;
	if(id == 0)
	{
		writeln("total time: ", elapsed);
	}
	
	return MPI_Shutdown;
}