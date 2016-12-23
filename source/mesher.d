/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.mesher;

import std.algorithm;
import std.array;
import std.conv;
import std.getopt;
import std.math;
import std.stdio;
import std.string;

import numd.linearalgebra.matrix;

import ebb.mesh;
import ebb.integrators;
import ebb.io;

void main(string[] args)
{
	bool convert;
	bool generate;
	bool info;
	string inFile;
	string outFile;
	
	auto result = getopt(args, std.getopt.config.caseSensitive, std.getopt.config.bundling,
								"c|convert", "Convert from su2/xflow mesh format to EbbCFD format (normal or extended)", &convert,
								"g|generate", "generate mesh", &generate,
								"f|file", "File to convert", &inFile,
								"o|ofile", "File to output", &outFile,
								"i|info", "Get mesh information", &info);

	if(result.helpWanted)
	{
		writeln("ebb-mesh options:");
		foreach(opt; result.options)
		{
			writeln(opt.optShort, " | ", opt.optLong, "\t\t", opt.help);
		}
		return;
	}

	if(info)
	{
		if(inFile == "")
		{
			writeln("No mesh file supplied, exiting");
			return;
		}

		UMesh2 mesh; 

		if(inFile.canFind("gri"))
		{
			writeln("Reading xflow mesh file");
			mesh = parseXflowMesh(inFile, false);
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
		else if(inFile.canFind("mmsh"))
		{
			writeln("Loading matlab mesh");
			mesh = loadMatlabMesh(inFile);
			return;
		}
		else
		{
			writeln("Unsupported mesh type, exiting");
			return;
		}

		mesh.buildMesh;

		auto dAve = mesh.cells.sum!".d"/mesh.cells.length;

		double x1 = -1.5;
		double y1 = 1.5;
		double x2 = 10;
		double y2 = -1.5;

		double dxSum = 0;
		uint dxCnt = 0;
		foreach(cell; mesh.cells)
		{
			if((cell.centroid[0] >= x1) && (cell.centroid[1] <= y1) && (cell.centroid[0] <= x2) && (cell.centroid[1] >= y2))
			{
				dxSum += cell.d;
				dxCnt++;
			}
		}

		double dxAve = dxSum/cast(double)dxCnt;

		writeln("----------------------------Mesh Info----------------------------");
		writeln("mesh: ", inFile);
		writeln("Dimensions: ", 2);
		writeln("cells: ", mesh.cells.length);
		writeln("D_ave: ", dAve);
		writeln("D_ave (bounded): ", dxAve);
		writeln("nodes: ", mesh.nodes.length);
		writeln("boundaries: ", mesh.bTags.length);
		for(uint i = 0; i < mesh.bTags.length; i++)
		{
			writeln("    boundary ", i, ": ", mesh.bTags[i]);
			writeln("        faces: ", mesh.bGroups[i].length);
		}
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
			saveMatlabMesh(mesh, outFile);
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