/+ Copyright (c) 2016 Robert F. Rau II +/
module mesher;

import std.array;
import std.conv;
import std.getopt;
import std.math : abs;
import std.stdio;
import std.string;

import mesh;

void initMesh(ref Mesh mesh, double width, double height, double dx, double dy)
{
	import std.math : round;
	import std.conv : to;
	int N = (width/dx).round.to!int;
	int M = (height/dy).round.to!int;
	
	mesh = Mesh(N,M);
	
	for(int i = 0; i < mesh.N; i++)
	{
		for(int j = 0; j < mesh.M; j++)
		{
			mesh[i, j].x = i*dx;
			mesh[i, j].y = j*dy;
			mesh[i, j].dx = dx;
			mesh[i, j].dy = dy;
		}
	}	
}

void addBox(ref Mesh mesh, double xt, double yt, double xb, double yb, Vec q, CellType bottom, CellType top, CellType right, CellType left)
{
	int xtIdx, xbIdx;
	int ytIdx, ybIdx;
	double ytDelta, ybDelta;
	double xtDelta, xbDelta;
	for(int i = 0; i < mesh.N; i++)
	{
		for(int j = 0; j < mesh.M; j++)
		{
			if(ytDelta > abs(yt - mesh[i,j].y))
			{
				ytIdx = j;
			}
			if(ybDelta > abs(yb - mesh[i,j].y))
			{
				ybIdx = j;
			}
			
			ytDelta = abs(yt - mesh[i,j].y);
			ybDelta = abs(yb - mesh[i,j].y);
		}
		
		if(xtDelta > abs(xt - mesh[i,0].x))
		{
			xtIdx = i;
		}
			
		if(xbDelta > abs(xb - mesh[i,0].x))
		{
			xbIdx = i;
		}
		xtDelta = abs(xt - mesh[i,0].x);
		xbDelta = abs(xb - mesh[i,0].x);
	}
	
	for(int i = xtIdx; i <= xbIdx; i++)
	{
		for(int j = ybIdx; j <= ytIdx; j++)
		{
			mesh[i,j].cellType = CellType.Normal;
			mesh[i,j].q = q;
			mesh[i,j].qR = q;
			mesh[i,j].qL = q;
			mesh[i,j].qT = q;
			mesh[i,j].qB = q;
		}
	}
	
	for(int i = xtIdx-1; i <= xbIdx+1; i++)
	{
		if(mesh[i,ytIdx+1].cellType == CellType.Solid)
		{
			mesh[i,ytIdx+1].cellType = top;
			mesh[i,ytIdx+1].q = q;
			mesh[i,ytIdx+1].qR = q;
			mesh[i,ytIdx+1].qL = q;
			mesh[i,ytIdx+1].qT = q;
			mesh[i,ytIdx+1].qB = q;
		}
		else if(mesh[i,ytIdx-1].cellType != CellType.Normal)
		{
			mesh[i, ytIdx-1].corner = true;
			mesh[i, ytIdx-1].cornerType = mesh[i,ybIdx-1].cellType;
			mesh[i, ytIdx-1].cellType = top;
		}
		
		if(mesh[i,ybIdx-1].cellType == CellType.Solid)
		{
			mesh[i,ybIdx-1].cellType = bottom;
			mesh[i,ybIdx-1].q = q;
			mesh[i,ybIdx-1].qR = q;
			mesh[i,ybIdx-1].qL = q;
			mesh[i,ybIdx-1].qT = q;
			mesh[i,ybIdx-1].qB = q;
		}
		else if((mesh[i,ybIdx-1].cellType != CellType.Normal) && ((mesh[i,ybIdx-1].cellType == CellType.GhostMirrorXL) || (mesh[i,ybIdx-1].cellType == CellType.GhostMirrorXR)
			|| (mesh[i,ybIdx-1].cellType == CellType.GhostNoGradXR) || (mesh[i,ybIdx-1].cellType == CellType.GhostNoGradXL))) 
		{
			mesh[i, ybIdx-1].corner = true;
			mesh[i, ybIdx-1].cornerType = mesh[i,ybIdx-1].cellType;
			mesh[i, ybIdx-1].cellType = bottom; 
		}
	}

	for(int j = ybIdx; j <= ytIdx; j++)
	{
		if(mesh[xtIdx-1,j].cellType == CellType.Solid)
		{
			mesh[xtIdx-1,j].cellType = left;
			mesh[xtIdx-1,j].q = q;
			mesh[xtIdx-1,j].qR = q;
			mesh[xtIdx-1,j].qL = q;
			mesh[xtIdx-1,j].qT = q;
			mesh[xtIdx-1,j].qB = q;
		}
		else if(mesh[xtIdx-1,j].cellType != CellType.Normal)
		{
			mesh[xtIdx-1,j].corner = true;
			mesh[xtIdx-1,j].cornerType = left;
		}
		
		if(mesh[xbIdx+1,j].cellType == CellType.Solid)
		{
			mesh[xbIdx+1,j].cellType = right;
			mesh[xbIdx+1,j].q = q;
			mesh[xbIdx+1,j].qR = q;
			mesh[xbIdx+1,j].qL = q;
			mesh[xbIdx+1,j].qT = q;
			mesh[xbIdx+1,j].qB = q;
		}
		else if(mesh[xbIdx+1,j].cellType != CellType.Normal)
		{
			mesh[xbIdx+1,j].corner = true;
			mesh[xbIdx+1,j].cornerType = right;
		}
	}
}

struct SU2Cell
{
	
}

enum GeomType
{
	Triangle,
	Quad
}

void parseSu2Mesh(string meshFile)
{
	auto file = File(meshFile);
	bool readingCells = false;
	bool readingNodes = false;
	uint numCells = 0;
	uint cellIdx = 0;
	uint numNodes = 0;
	uint nodeIdx = 0;
	uint numDims = 0;

	uint[] elements;

	foreach(line; file.byLine)
	{
		auto cleanLine = line.strip.chomp;
		if(cleanLine.indexOf("NDIME") > -1)
		{
			numDims = cleanLine.split("=")[$-1].strip.to!uint;
			writeln("Number of dimensions = "~numDims.to!string);
		}
		else if(cleanLine.indexOf("NELEM") > -1)
		{
			numCells = cleanLine.split("=")[$-1].strip.to!uint;
			readingCells = true;
			writeln("Number of cells = "~numCells.to!string);
		}
		else if(readingCells)
		{

		}
	}
}

void main(string[] args)
{
	Mesh mesh;
	bool convert;
	bool generate;
	string inFile;
	string outFile;
	
	auto result = getopt(args, "c|convert", "Convert from su2 mesh format", &convert,
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
		parseSu2Mesh(inFile);
	}
	else if(generate)
	{
		/+
		initMesh(mesh, 50, 40, 0.5, 0.5);
		//addBox(ref Mesh mesh, double xt, double yt, double xb, double yb, Vec q, CellType bottom, CellType top, CellType right, CellType left)
		addBox(mesh, 0.5, 22.0, 3.0, 18.0, buildQ(1.6367, 0.178, 0, 1.245), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostConst, CellType.GhostConst);
		addBox(mesh, 3.5, 22.0, 15.0, 18.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostConst, CellType.GhostConst);
		addBox(mesh, 15.0, 39.0, 49.0, 0.5, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostNoGradYT, CellType.GhostNoGradYB, CellType.GhostNoGradXL, CellType.GhostMirrorXR);
		
		addBox(mesh, 1.0, 39.0, 13.0, 30.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 6.0, 39.0, 10.0, 22.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		
		addBox(mesh, 1.0, 10.0, 13.0, 0.5, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 6.0, 19.0, 10.0, 10.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		mesh.updateGhosts();
		printMesh(mesh);
		saveMesh(mesh, "dissipationTubes.mesh", 0.01, 0);
		+/
		/+
		initMesh(mesh, 50, 20, 0.5, 0.5);
		//addBox(ref Mesh mesh, double xt, double yt, double xb, double yb, Vec q, CellType bottom, CellType top, CellType right, CellType left)
		addBox(mesh, 0.5, 2.0, 3.0, 0.5, buildQ(1.6367, 0.178, 0, 1.245), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostConst, CellType.GhostConst);
		addBox(mesh, 3.5, 2.0, 15.0, 0.5, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostConst, CellType.GhostConst);
		addBox(mesh, 15.0, 19.0, 49.0, 0.5, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostNoGradYB, CellType.GhostNoGradXL, CellType.GhostMirrorXR);
		
		addBox(mesh, 1.0, 19.0, 13.0, 10.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 6.0, 10.0, 10.0, 2.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		
		
		//addBox(mesh, 1.0, 10.0, 13.0, 0.5, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		//addBox(mesh, 6.0, 18.0, 10.0, 10.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		mesh.updateGhosts();
		printMesh(mesh);
		saveMesh(mesh, "mirrorMesh.mesh", 0.01, 0);
		+/
		
		/+
		initMesh(mesh, 60, 40, 0.5, 0.5);
		
		// initial tunnel
		addBox(mesh, 0.5, 22.0, 3.0, 18.0, buildQ(1.6367, 0.178, 0, 1.245), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostConst, CellType.GhostConst);
		addBox(mesh, 3.5, 22.0, 25.0, 18.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostConst, CellType.GhostConst);
		addBox(mesh, 25.0, 39.0, 59.0, 0.5, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostNoGradYT, CellType.GhostNoGradYB, CellType.GhostNoGradXL, CellType.GhostMirrorXR);
		
		addBox(mesh, 19.0, 35.0, 23.0, 22.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 6.0, 35.0, 10.0, 22.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 6.5, 35.0, 22.5, 30.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		
		addBox(mesh, 6.5, 10.0, 22.5, 5.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 6.0, 19.0, 10.0, 5.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 19.0, 19.0, 23.0, 5.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		
		mesh.updateGhosts();
		printMesh(mesh);
		saveMesh(mesh, "dissipationTunnel.mesh", 0.01, 0);
		+/
		/+
		initMesh(mesh, 80, 40, 0.5, 0.5);
		
		// initial tunnel
		addBox(mesh, 0.5, 22.0, 3.0, 18.0, buildQ(1.6367, 0.178, 0, 1.245), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostConst, CellType.GhostConst);
		addBox(mesh, 3.5, 22.0, 38.0, 18.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostConst, CellType.GhostConst);
		addBox(mesh, 38.0, 39.0, 79.0, 0.5, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostNoGradYT, CellType.GhostNoGradYB, CellType.GhostNoGradXL, CellType.GhostMirrorXR);
		
		addBox(mesh, 19.0, 35.0, 23.0, 22.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 6.0, 35.0, 10.0, 22.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 32.0, 35.0, 36.0, 22.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 6.5, 35.0, 35.5, 30.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		
		addBox(mesh, 6.5, 10.0, 35.5, 5.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 6.0, 19.0, 10.0, 5.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 19.0, 19.0, 23.0, 5.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 32.0, 19.0, 36.0, 5.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		
		mesh.updateGhosts();
		printMesh(mesh);
		saveMesh(mesh, "dissipationTunnel2.mesh", 0.01, 0);
		+/
		/+
		initMesh(mesh, 100, 150, 0.5, 0.5);
		
		addBox(mesh, 0.5, 149.0, 4.0, 0.5, buildQ(1.4, 2, 0, 1), CellType.GhostNoGradYT, CellType.GhostNoGradYB, CellType.GhostConst, CellType.GhostConst);

		addBox(mesh, 44.0, 70.0, 56.0, 0.5, buildQ(1.4, 0, 0, 1), CellType.GhostNoGradYT, CellType.GhostMirrorYB, CellType.GhostNoGradXL, CellType.GhostNoGradXR);
		addBox(mesh, 4.0, 149.0, 45.0, 0.5, buildQ(1.4, 0, 0, 1), CellType.GhostNoGradYT, CellType.GhostNoGradYB, CellType.GhostMirrorXL, CellType.GhostNoGradXR);
		addBox(mesh, 45.0, 149.0, 55.0, 80.0, buildQ(1.4, 0, 0, 1), CellType.GhostMirrorYT, CellType.GhostNoGradYB, CellType.GhostNoGradXL, CellType.GhostNoGradXR);
		addBox(mesh, 55.0, 149.0, 99.0, 0.5, buildQ(1.4, 0, 0, 1), CellType.GhostNoGradYT, CellType.GhostNoGradYB, CellType.GhostNoGradXL, CellType.GhostMirrorXR);
		
		mesh.updateGhosts();
		printMesh(mesh);
		saveMesh(mesh, "box.mesh", 0.01, 0);
		+/
		/+
		initMesh(mesh, 100, 100, 0.5, 0.5);
		import std.math : sqrt;
		
		addBox(mesh, 45.0, 98.0, 55.0, 55.0, buildQ(1.4, 0, 0, 1), CellType.GhostMirrorYT, CellType.GhostNoGradYB, CellType.GhostNoGradXL, CellType.GhostNoGradXR);
		addBox(mesh, 3.0, 99.0, 98.5, 98.0, buildQ(1.4, sqrt(100.0), -sqrt(100.0), 1), CellType.GhostNoGradYT, CellType.GhostConst, CellType.GhostNoGradXL, CellType.GhostConst);

		addBox(mesh, 44.0, 45.0, 56.0, 0.5, buildQ(1.4, 0, 0, 1), CellType.GhostNoGradYT, CellType.GhostMirrorYB, CellType.GhostNoGradXL, CellType.GhostNoGradXR);
		addBox(mesh, 2.0, 98.0, 45.0, 0.5, buildQ(1.4, 0, 0, 1), CellType.GhostNoGradYT, CellType.GhostNoGradYB, CellType.GhostMirrorXL, CellType.GhostNoGradXR);
		addBox(mesh, 0.5, 99.0, 4.0, 0.5, buildQ(1.4, sqrt(100.0), -sqrt(100.0), 1), CellType.GhostNoGradYT, CellType.GhostConst, CellType.GhostConst, CellType.GhostConst);
		addBox(mesh, 55.0, 98.0, 99.0, 0.5, buildQ(1.4, 0, 0, 1), CellType.GhostNoGradYT, CellType.GhostNoGradYB, CellType.GhostNoGradXL, CellType.GhostMirrorXR);
		
		mesh[89, 195].corner = false;
		mesh[89, 195].cornerType = CellType.Normal;
		mesh[111, 195].corner = false;
		mesh[111, 195].cornerType = CellType.Normal;
		
		mesh.updateGhosts();
		printMesh(mesh);
		saveMesh(mesh, "box2.mesh", 0.01, 0);
		+/

		initMesh(mesh, 1.0 + 0.640, 0.026, 0.00005, 0.00005);
		import std.math : sqrt;
		
		// void addBox(ref Mesh mesh, double xt, double yt, double xb, double yb, Vec q, CellType bottom, CellType top, CellType right, CellType left)
		// Vec buildQ(double rho, double u, double v, double p)

		addBox(mesh, 1.0 + 0.566, 0.02225, 1.0 + 0.636, 0.01975, buildQ(1.2, 0, 0, 101325), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 1.0 + 0.566, 0.00625, 1.0 + 0.616, 0.00375, buildQ(1.2, 0, 0, 101325), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);

		addBox(mesh, 1.0 + 0.509, 0.025, 1.0 + 0.566, 0.017, buildQ(1.2, 0, 0, 101325), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 1.0 + 0.509, 0.009, 1.0 + 0.566, 0.001, buildQ(1.2, 0, 0, 101325), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);

		addBox(mesh, 0.001, 0.017, 1.0 + 0.0011, 0.009, buildQ(8.195, 0, 0, 689476), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		addBox(mesh, 1.0 + 0.0011, 0.017, 1.0 + 0.501, 0.009, buildQ(1.2, 0, 0, 101325), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostNoGradXL, CellType.GhostNoGradXR);
		addBox(mesh, 1.0 + 0.501, 0.025, 1.0 + 0.509, 0.001, buildQ(1.2, 0, 0, 101325), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostMirrorXL, CellType.GhostMirrorXR);
		
/+
		addBox(mesh, 2.0, 98.0, 45.0, 0.5, buildQ(1.2, 0, 0, 101325), CellType.GhostNoGradYT, CellType.GhostNoGradYB, CellType.GhostMirrorXL, CellType.GhostNoGradXR);
		addBox(mesh, 0.5, 99.0, 4.0, 0.5, buildQ(1.2, 0, 0, 101325), CellType.GhostNoGradYT, CellType.GhostConst, CellType.GhostConst, CellType.GhostConst);
		addBox(mesh, 55.0, 98.0, 99.0, 0.5, buildQ(1.2, 0, 0, 101325), CellType.GhostNoGradYT, CellType.GhostNoGradYB, CellType.GhostNoGradXL, CellType.GhostMirrorXR);
+/
/+
		uint i = 11;
		for(uint j = 0; j < mesh.M; j++)
		{
			if(mesh[i,j].corner)
			{
				mesh[i,j].corner = false;
				writeln("corner");
			}
		}
+/
		mesh.updateGhosts();
		printMesh(mesh, "pressureSensorSytem.txt");
		writeln("saving mesh");
		saveMesh(mesh, "pressureSensorSytem.mesh", 0.01, 0);

		/+
		initMesh(mesh, 20, 20, 0.5, 0.5);
		//addBox(ref Mesh mesh, double xt, double yt, double xb, double yb, Vec q, CellType bottom, CellType top, CellType right, CellType left)
		addBox(mesh, 0.5, 12.0, 3.0, 8.0, buildQ(1.6367, 0.178, 0, 1.245), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostConst, CellType.GhostConst);
		addBox(mesh, 3.5, 12.0, 5.0, 8.0, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostMirrorYT, CellType.GhostMirrorYB, CellType.GhostConst, CellType.GhostConst);
		addBox(mesh, 5.0, 19.0, 19.0, 0.5, buildQ(1.4, 0.0, 0, 1.0), CellType.GhostNoGradYT, CellType.GhostNoGradYB, CellType.GhostNoGradXL, CellType.GhostMirrorXR);
		
		mesh.updateGhosts();
		printMesh(mesh);
		saveMesh(mesh, "standard.mesh", 0.01, 0);
		+/
	}
}