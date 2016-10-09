/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.mesher;

import std.algorithm;
import std.array;
import std.conv;
import std.exception;
import std.getopt;
import std.math;
import std.stdio;
import std.string;

import numd.linearalgebra.matrix;

import ebb.mesh;

void main(string[] args)
{
	bool convert;
	bool generate;
	string inFile;
	string outFile;
	
	auto result = getopt(args, "c|convert", "Convert from su2/xflow mesh format to EbbCFD format (normal or extended)", &convert,
								"g|generate", "generate mesh", &generate,
								"f|file", "File to convert", &inFile,
								"o|ofile", "File to output", &outFile);

	if(result.helpWanted)
	{
		writeln("ebb-mesh options:");
		foreach(opt; result.options)
		{
			writeln(opt.optShort, " | ", opt.optLong, "\t\t", opt.help);
		}
		return;
	}

	if(convert)
	{
		if(inFile == "")
		{
			writeln("No file supplied, exiting");
			return;
		}

		if(outFile == "")
		{
			writeln("No output file supplied, exiting");
			return;
		}

		UMesh2 mesh;
		if(inFile.canFind("gri"))
		{
			writeln("Reading xflow mesh file");
			mesh = parseXflowMesh(inFile);
		}
		else if(inFile.canFind("su2"))
		{
			//mesh = parseSu2Mesh(inFile);
			writeln("SU2 mesh not yet supported, exiting");
			return;
		}
		else if(inFile.canFind("emsh"))
		{
			writeln("EbbCFD mesh not yet supported lol, exiting");
			return;
		}
		else
		{
			writeln("Unsupported mesh type, exiting");
			return;
		}

		if(outFile.canFind("emsh"))
		{
			writeln("EbbCFD mesh not yet supported lol, exiting");
			return;
		}
		else if(outFile.canFind("mmsh"))
		{
			writeln("Save mesh in matlab format");
			saveMatlabMesh(mesh, outFile.ptr);
		}
		else
		{
			writeln("Unsupported mesh output type, exiting");
			return;
		}
	}
	else if(generate)
	{
		writeln("Unsupported, exiting");
		return;
	}
}