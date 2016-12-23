module ebb.io;

import core.stdc.stdio;

import std.algorithm;
import std.array;
import std.conv;
import std.file;
import std.meta;
import std.stdio : File, writeln;
import std.string;

import numd.linearalgebra.matrix;

import ebb.exception;
import ebb.mesh;

/+
@nogc writeln(T...)(T args)
{
	//string printfCall = "printf(\"\");";
	char[512] printfCall;
	foreach(arg; AliasSeq!args)
	{
		//printfCall ~= arg.stringof;
	}

	//mixin(printfCall);
}
+/
struct MeshHeader
{
	const uint meshMagic = 0xB1371AC7;
	uint meshVersion;
	uint dims;
	uint nNodes;
	uint nElems;
	uint nBGroups;
}

/+
@nogc void saveMesh(ref UMesh2 mesh, char* filename)
{
	MeshHeader header = {meshVersion: 1, dims: 2, nNodes: cast(uint)mesh.nodes.length, nElems: cast(uint)mesh.elements.length, nBGroups: cast(uint)mesh.bTags.length};

	ulong totSize = MeshHeader.sizeof;
	totSize += header.nNodes*header.dims*double.sizeof;
	for(uint i = 0; i < mesh.elements.length; i++)
	{
		totSize += 
	}
}
+/

struct SlnHeader
{
	static const uint slnMagic = 0xEA98E1F5;
	uint slnVersion;
	uint dataPoints;
	double t;
	double dt;
}

@nogc void saveSolution(ref UMesh2 mesh, char* filename, double t, double dt)
{
	import std.experimental.allocator.mallocator : Mallocator;
	import std.bitmanip : write;

	SlnHeader header = {slnVersion: 1, dataPoints: cast(uint)mesh.interiorCells.length, t: t, dt: dt};

	ulong totSize = SlnHeader.sizeof + header.dataPoints*4*double.sizeof + uint.sizeof + uint.sizeof;

	ubyte[] buffer = cast(ubyte[])Mallocator.instance.allocate(totSize);
	scope(exit) Mallocator.instance.deallocate(buffer);

	size_t offset = 0;
	buffer.write!uint(header.slnMagic, &offset);
	buffer.write!uint(header.slnVersion, &offset);
	buffer.write!uint(header.dataPoints, &offset);
	buffer.write!double(header.t, &offset);
	buffer.write!double(header.dt, &offset);

	//for(uint i = 0; i < mesh.cells.length; i++)
	foreach(i; mesh.interiorCells)
	{
		buffer.write!double(mesh.q[i][0], &offset);
		buffer.write!double(mesh.q[i][1], &offset);
		buffer.write!double(mesh.q[i][2], &offset);
		buffer.write!double(mesh.q[i][3], &offset);
	}

	import std.digest.crc : CRC32;

	CRC32 crc;
	crc.start();
	crc.put(buffer);
	auto crc32 = crc.finish();

	buffer.write!ubyte(crc32[0], &offset);
	buffer.write!ubyte(crc32[1], &offset);
	buffer.write!ubyte(crc32[2], &offset);
	buffer.write!ubyte(crc32[3], &offset);

	auto file = fopen(filename, "wb");
	ulong writeOffset = 0;
	while(writeOffset < buffer.length)
	{
		ulong chunkSize = 1024*1024*1024;
		if(buffer.length - writeOffset < chunkSize)
		{
			chunkSize = buffer.length - writeOffset;
		}
		fwrite(buffer[writeOffset..writeOffset+chunkSize].ptr, ubyte.sizeof, chunkSize, file);
		writeOffset += chunkSize;
	}
	fclose(file);
}

@nogc void saveLimits(ref UMesh2 mesh, char* filename, double t, double dt)
{
	import std.experimental.allocator.mallocator : Mallocator;
	import std.bitmanip : write;

	SlnHeader header = {slnVersion: 1, dataPoints: cast(uint)mesh.interiorCells.length, t: t, dt: dt};

	ulong totSize = SlnHeader.sizeof + header.dataPoints*4*double.sizeof + uint.sizeof + uint.sizeof;

	ubyte[] buffer = cast(ubyte[])Mallocator.instance.allocate(totSize);
	scope(exit) Mallocator.instance.deallocate(buffer);

	size_t offset = 0;
	buffer.write!uint(header.slnMagic, &offset);
	buffer.write!uint(header.slnVersion, &offset);
	buffer.write!uint(header.dataPoints, &offset);
	buffer.write!double(header.t, &offset);
	buffer.write!double(header.dt, &offset);

	//for(uint i = 0; i < mesh.cells.length; i++)
	foreach(i; mesh.interiorCells)
	{
		buffer.write!double(mesh.cells[i].gradErr[0], &offset);
		buffer.write!double(mesh.cells[i].gradErr[1], &offset);
		buffer.write!double(mesh.cells[i].gradErr[2], &offset);
		buffer.write!double(mesh.cells[i].gradErr[3], &offset);
	}

	import std.digest.crc : CRC32;

	CRC32 crc;
	crc.start();
	crc.put(buffer);
	auto crc32 = crc.finish();

	buffer.write!ubyte(crc32[0], &offset);
	buffer.write!ubyte(crc32[1], &offset);
	buffer.write!ubyte(crc32[2], &offset);
	buffer.write!ubyte(crc32[3], &offset);

	auto file = fopen(filename, "wb");
	ulong writeOffset = 0;
	while(writeOffset < buffer.length)
	{
		ulong chunkSize = 1024*1024*1024;
		if(buffer.length - writeOffset < chunkSize)
		{
			chunkSize = buffer.length - writeOffset;
		}
		fwrite(buffer[writeOffset..writeOffset+chunkSize].ptr, ubyte.sizeof, chunkSize, file);
		writeOffset += chunkSize;
	}
	fclose(file);
}

@nogc bool loadSolution(ref UMesh2 mesh, ref double t, ref double dt, string filename)
{
	import std.algorithm : canFind;
	import std.bitmanip : peek;
	import std.experimental.allocator.mallocator : Mallocator;
	import std.digest.crc : CRC32;
	import std.string : toStringz;
	char[1024] filenamePtr;
	filenamePtr[] = 0;
	filenamePtr[0..filename.length] = filename[];
	auto file = fopen(filenamePtr.ptr, "rb");
	assert(file != null);

	ubyte[] buffer = cast(ubyte[])Mallocator.instance.allocate(8);
	scope(exit) Mallocator.instance.deallocate(buffer);

	CRC32 crc;
	crc.start();

	fread(buffer.ptr, 1, 4, file);
	crc.put(buffer[0..4]);
	uint slnMagic = buffer.peek!uint;

	if(slnMagic != SlnHeader.slnMagic)
	{
		return false;
	}

	fread(buffer.ptr, 1, 4, file);
	crc.put(buffer[0..4]);
	uint slnVersion = buffer.peek!uint;
	fread(buffer.ptr, 1, 4, file);
	crc.put(buffer[0..4]);
	uint dataPoints = buffer.peek!uint;

	if(dataPoints != mesh.interiorCells.length)
	{
		return false;
	}

	fread(buffer.ptr, 1, 8, file);
	crc.put(buffer[]);
	t = buffer.peek!double;
	fread(buffer.ptr, 1, 8, file);
	crc.put(buffer[]);
	dt = buffer.peek!double;

	//for(uint i = 0; i < mesh.q.length; i++)
	foreach(i; mesh.interiorCells)
	{
		fread(buffer.ptr, 1, 8, file);
		crc.put(buffer[]);
		mesh.q[i][0] = buffer.peek!double;
		fread(buffer.ptr, 1, 8, file);
		crc.put(buffer[]);
		mesh.q[i][1] = buffer.peek!double;
		fread(buffer.ptr, 1, 8, file);
		crc.put(buffer[]);
		mesh.q[i][2] = buffer.peek!double;
		fread(buffer.ptr, 1, 8, file);
		crc.put(buffer[]);
		mesh.q[i][3] = buffer.peek!double;
	}
	//ubyte[4] readCrc32;
	fread(buffer.ptr, 1, 4, file);
	auto crc32 = crc.finish();
	bool crcGood = true;
	/+
	for(uint i = 0; i < 4; i++)
	{
		crcGood &= crc32[i] == buffer[i]; 
	}
	
	printf("%x %x %x %x\n", crc32[0], crc32[1], crc32[2], crc32[3]);
	printf("%x %x %x %x\n", buffer[0], buffer[1], buffer[2], buffer[3]);
	+/
	fclose(file);

	return crcGood;
}

void loadMatlabSolution(ref UMesh2 mesh, string filename)
{
	//UMesh2 mesh;
	import std.algorithm : canFind;
	import std.bitmanip : read;
	import std.conv : to;

	//writeln("Reading file ", filename);
	auto slnFile = File(filename);

	auto fileSize = slnFile.size;

	auto buffer = slnFile.rawRead(new ubyte[fileSize]);
	auto nNodes = buffer.read!ulong;
	enforce(nNodes == mesh.nodes.length, "Mesh file has different number of nodes than solution file.");

	//writeln("nNodes = ", nNodes);

	for(uint i = 0; i < nNodes; i++)
	{
		buffer.read!(double);
		buffer.read!(double);
	}

	auto nEls = buffer.read!ulong;
	enforce(nEls == mesh.cells.length, "Mesh file has different number of cells than solution file.");

	//writeln("nEdges = ", nEls);
	for(uint i = 0; i < nEls; i++)
	{
		buffer.read!(double);
		buffer.read!(double);
		buffer.read!(double);
	}

	auto nIe = buffer.read!ulong;
	//writeln("nIe = ", nIe);
	for(uint i = 0; i < nIe; i++)
	{
		buffer.read!(double);
		buffer.read!(double);
		buffer.read!(double);
		buffer.read!(double);
	}

	auto nBe = buffer.read!ulong;
	//writeln("nBe = ", nBe);
	for(uint i = 0; i < nBe; i++)
	{
		buffer.read!(double);
		buffer.read!(double);
		buffer.read!(double);
		buffer.read!(double);
	}

	auto nTags = buffer.read!ulong;
	//writeln("nTags = ", nTags);
	for(uint i = 0; i < nTags; i++)
	{
		auto strLen = buffer.read!uint;
		for(uint j = 0; j < strLen; j++)
		{
			auto str = buffer.read!char;
		}
	}

	auto nCells = buffer.read!ulong;
	//writeln("nCells = ", nCells);
	enforce(nEls == mesh.cells.length, "Mesh file has different number of cells than solution file.");

	for(uint i = 0; i < nCells; i++)
	{
		mesh.q[i][0] = buffer.read!(double);
		mesh.q[i][1] = buffer.read!(double);
		mesh.q[i][2] = buffer.read!(double);
		mesh.q[i][3] = buffer.read!(double);
	}

}

UMesh2 loadMatlabMesh(string filename)
{
	UMesh2 mesh;
	import std.algorithm : canFind;
	import std.bitmanip : read;
	import std.conv : to;

	//writeln("Reading file ", filename);
	auto meshFile = File(filename);

	auto fileSize = meshFile.size;

	auto buffer = meshFile.rawRead(new ubyte[fileSize]);
	auto nNodes = buffer.read!ulong;

	for(uint i = 0; i < nNodes; i++)
	{
		mesh.nodes ~= [buffer.read!(double), buffer.read!(double)];
	}

	auto nEls = buffer.read!ulong;

	for(uint i = 0; i < nEls; i++)
	{
		mesh.elements ~= [cast(uint)buffer.read!(double), cast(uint)buffer.read!(double), cast(uint)buffer.read!(double)];
	}

	auto nIe = buffer.read!ulong;
	//writeln("nIe = ", nIe);
	for(uint i = 0; i < nIe; i++)
	{
		buffer.read!(double);
		buffer.read!(double);
		buffer.read!(double);
		buffer.read!(double);
	}

	auto nBe = buffer.read!ulong;
	uint[] bGroup;
	//writeln("nBe = ", nBe);
	for(uint i = 0; i < nBe; i++)
	{
		mesh.bNodes ~= [cast(uint)buffer.read!(double) - 1, cast(uint)buffer.read!(double) - 1];
		buffer.read!(double);
		bGroup ~= cast(uint)buffer.read!(double);
	}

	writeln(bGroup);
	auto nTags = buffer.read!ulong;
	//writeln("nTags = ", nTags);
	for(uint i = 0; i < nTags; i++)
	{
		auto strLen = buffer.read!uint;
		string str;
		for(uint j = 0; j < strLen; j++)
		{
			str ~= buffer.read!char;
		}
		mesh.bTags ~= str;
	}

	writeln(mesh.bTags);
	/+
	f = fopen(meshFile, 'rb');
	
	nNodes = fread(f, 1, 'uint64', 'b');
	nodes = fread(f, [2, nNodes], 'double', 'b')';
	
	nEls = fread(f, 1, 'uint64', 'b');
	e2n = fread(f, [3, nEls], 'double', 'b')';
	
	ieSize = fread(f, 1, 'uint64', 'b');
	ie = fread(f, [4, ieSize], 'double', 'b')';
	
	beSize = fread(f, 1, 'uint64', 'b');
	be = fread(f, [4, beSize], 'double', 'b')';
	
	nTags = fread(f, 1, 'uint64', 'b');
	for i = 1:nTags
		tagLen = fread(f, 1, 'uint32', 'b');
		tags{i} = fread(f, tagLen, '*char');
	end
	fclose(f);
	+/
	return mesh;
}
@nogc void saveMatlabMesh(ref UMesh2 mesh, immutable (string) filename)
{
	import std.experimental.allocator.mallocator : Mallocator;
	import std.bitmanip : write;

	ulong nodesSize = cast(ulong)(2*mesh.nodes.length*double.sizeof);
	ulong e2nSize = cast(ulong)(3*mesh.elements.length*double.sizeof);
	ulong ieSize = 0;
	ulong beSize = 0;
	for(uint i = 0; i < mesh.edges.length; i++)
	{
		if(!mesh.edges[i].isBoundary)
		{
			ieSize += 4;
		}
		else
		{
			beSize += 4;
		}
	}
	
	ieSize *= double.sizeof;
	beSize *= double.sizeof;

	ulong bNameSize = 0;
	for(uint i = 0; i < mesh.bTags.length; i++)
	{
		bNameSize += (uint.sizeof + mesh.bTags[i].length);  
	}

	immutable ulong nodeHeaderSize = ulong.sizeof;
	immutable ulong e2nHeaderSize = ulong.sizeof;
	immutable ulong ieHeaderSize = ulong.sizeof;
	immutable ulong beHeaderSize = ulong.sizeof;
	immutable ulong tagHeaderSize = ulong.sizeof;

	ulong totSize = nodeHeaderSize + e2nHeaderSize + ieHeaderSize + beHeaderSize + tagHeaderSize;
	totSize += nodesSize + e2nSize + ieSize + beSize + bNameSize;
	ubyte[] buffer = cast(ubyte[])Mallocator.instance.allocate(totSize);
	scope(exit) Mallocator.instance.deallocate(buffer);

	size_t offset = 0;

	buffer.write!ulong((nodesSize/(2*double.sizeof)), &offset);
	for(uint i = 0; i < mesh.nodes.length; i++)
	{
		buffer.write!(double)(mesh.nodes[i][0], &offset);
		buffer.write!(double)(mesh.nodes[i][1], &offset);
	}

	buffer.write!ulong((e2nSize/(3*double.sizeof)), &offset);
	for(uint i = 0; i < mesh.elements.length; i++)
	{
		buffer.write!(double)(mesh.elements[i][0], &offset);
		buffer.write!(double)(mesh.elements[i][1], &offset);
		buffer.write!(double)(mesh.elements[i][2], &offset);
	}

	buffer.write!ulong((ieSize/(4*double.sizeof)), &offset);
	for(uint i = 0; i < mesh.edges.length; i++)
	{
		if(!mesh.edges[i].isBoundary)
		{
			buffer.write!(double)(mesh.edges[i].nodeIdx[0] + 1, &offset);
			buffer.write!(double)(mesh.edges[i].nodeIdx[1] + 1, &offset);
			buffer.write!(double)(mesh.edges[i].cellIdx[0] + 1, &offset);
			buffer.write!(double)(mesh.edges[i].cellIdx[1] + 1, &offset);
		}
	}

	buffer.write!ulong((beSize/(4*double.sizeof)), &offset);
	for(uint i = 0; i < mesh.edges.length; i++)
	{
		if(mesh.edges[i].isBoundary)
		{
			buffer.write!(double)(mesh.edges[i].nodeIdx[0] + 1, &offset);
			buffer.write!(double)(mesh.edges[i].nodeIdx[1] + 1, &offset);
			buffer.write!(double)(mesh.edges[i].cellIdx[0] + 1, &offset);
			auto bgIdx = mesh.bTags.countUntil(mesh.edges[i].boundaryTag) + 1;
			buffer.write!(double)(bgIdx, &offset);
		}
	}

	buffer.write!ulong(mesh.bTags.length, &offset);
	for(uint i = 0; i < mesh.bTags.length; i++)
	{
		buffer.write!uint(cast(uint)mesh.bTags[i].length, &offset);
		for(uint j = 0; j < mesh.bTags[i].length; j++)
		{
			buffer.write!char(mesh.bTags[i][j], &offset);
		}
	}

	char[1024] filenamePtr;
	filenamePtr[] = 0;
	filenamePtr[0..filename.length] = filename[];
	auto file = fopen(filenamePtr.ptr, "wb");
	ulong writeOffset = 0;
	while(writeOffset < buffer.length)
	{
		ulong chunkSize = 1024*1024*1024;
		if(buffer.length - writeOffset < chunkSize)
		{
			chunkSize = buffer.length - writeOffset;
		}
		fwrite(buffer[writeOffset..writeOffset+chunkSize].ptr, ubyte.sizeof, chunkSize, file);
		writeOffset += chunkSize;
	}
	fclose(file);
}

char[][] readCleanLine(ref File file)
{
	return file.readln.strip.chomp.detab.split(' ').filter!(a => a != "").array;
}

UMesh2 parseXflowMesh(string meshFile, bool chatty = true)
{
	UMesh2 mesh;
	auto file = File(meshFile);

	auto headerLine = file.readCleanLine;
	immutable uint nNodes = headerLine[0].to!uint;
	immutable uint nElems = headerLine[1].to!uint;
	immutable uint nDims = headerLine[2].to!uint;

	enforce(nDims == 2, nDims.to!string~" dimensional meshes not supported");

	if(chatty)
	{
		writeln("Importing XFlow mesh "~meshFile);
		writeln("    nNodes = ", nNodes);
		writeln("    nElems = ", nElems);
		writeln("    nDims = ", nDims);
	}

	for(uint i = 0; i < nNodes; i++)
	{
		mesh.nodes ~= file.readCleanLine.to!(double[]);
	}

	immutable uint nBoundaryGroups = file.readCleanLine[0].to!uint;
	if(chatty) writeln("    nBoundaryGroups = ", nBoundaryGroups);

	uint[][] bNodes;
	size_t[] bGroupStart;
	string[] bTags;
	for(uint i = 0; i < nBoundaryGroups; i++)
	{
		auto boundaryHeader = file.readCleanLine;
		immutable uint bFaces = boundaryHeader[0].to!uint;
		immutable uint nodesPerFace = boundaryHeader[1].to!uint;
		string faceTag;
		if(boundaryHeader.length == 3)
		{
			faceTag = boundaryHeader[2].to!string;
			bTags ~= faceTag;
		}

		if(chatty) writeln("        Boundary group ", i, ": faces = ", bFaces, ", nodes per face = ", nodesPerFace, ", tag = ", faceTag);

		bGroupStart ~= bNodes.length;
		for(uint j = 0; j < bFaces; j++)
		{
			auto bn = file.readCleanLine.to!(uint[]);
			bNodes ~= [bn[0] - 1, bn[$-1] - 1];
		}
	}

	uint[][] elements;
	uint foundElements;
	uint faces;
	uint eGroup;
	while(foundElements < nElems)
	{
		char[][] elementLine = file.readCleanLine;
		uint q = elementLine[1].to!uint;
		uint subElements = elementLine[0].to!uint;

		enforce((q < 4) && (q != 0), "Unsuported q");

		if(elementLine[2].canFind("Tri"))
		{
			faces = 3;
		}
		else if(elementLine[2].canFind("Quad"))
		{
			faces = 4;
		}
		else
		{
			enforce(false, "Unsuported cell type");
		}

		if(chatty) writeln("    Element group ", eGroup, ": faces = ", faces, ", q = ", q, ", subElements = ", subElements);

		for(uint i = 0; i < subElements; i++)
		{
			auto els = file.readCleanLine.to!(uint[]);

			if(q == 1)
			{
				elements ~= els;
			}
			else if(q == 2)
			{
				if(faces == 3)
				{
					elements ~= [els[0], els[2], els[5]];
				}
				else
				{
					elements ~= [els[0], els[2], els[8], els[6]];
				}
			}
			else
			{
				if(faces == 3)
				{
					elements ~= [els[0], els[3], els[9]];
				}
				else
				{
					elements ~= [els[0], els[3], els[15], els[12]];
				}
			}
		}

		foundElements += subElements;
		eGroup++;
	}
	
	mesh.cells = new UCell2[nElems];
	mesh.q = new Vector!4[nElems];

	for(uint i = 0; i < nElems; i++)
	{
		mesh.cells[i].nEdges = faces;
	}
	mesh.elements = elements;
	mesh.bNodes = bNodes;
	mesh.bGroupStart = bGroupStart;
	mesh.bTags = bTags;
	//mesh.buildMesh;

	return mesh;
}

UMesh2 parseSu2Mesh(string meshFile)
{
	UMesh2 mesh;
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

	return mesh;
}