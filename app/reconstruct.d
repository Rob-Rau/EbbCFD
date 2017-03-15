/+ Copyright (c) 2017 Robert F. Rau II +/
module ebb.reconstruct;

import std.algorithm : canFind, filter, sort;
import std.array : array, split;
import std.conv;
import std.file : dirEntries, DirEntry, mkdir, readText, SpanMode;
import std.getopt;
import std.math;
import std.stdio;
import std.string;

import numd.linearalgebra.matrix;

import ebb.config;
import ebb.euler;
import ebb.exception;
import ebb.flux;
import ebb.integrators;
import ebb.limiters;
import ebb.mesh;
import ebb.manufacturedsolution;
import ebb.io;
import ebb.mpid;
import ebb.solve;

void stepMesh(ref UMesh2 mesh, Config config, double t, double dt)
{
	try
	{
		mesh.comm = MPI_COMM_SELF;
		mesh.mpiRank = 0;
		// Setup bc's
		double M = config.ic[0];
		double aoa = config.ic[1] * (PI/180);
		double rho = config.ic[3];
		double U = 1.0;
		double a = U/M;
		double p = (rho*a^^2.0)/gamma;
		double u = U*cos(aoa);
		double v = U*sin(aoa);
		if(config.viscosity)
		{
			config.physicalConfig.mu = (rho*U*config.physicalConfig.L)/config.physicalConfig.Re;
			writeln("mu = ", config.physicalConfig.mu);
		}

		for(uint i = 0; i < mesh.bGroups.length; i++)
		{
			@nogc uint findBcIndex(string tag)
			{
				auto tagIdx = config.boundaries.countUntil!"a.bTag == b"(tag);
				if(tagIdx < 0)
				{
					char[64] str;
					str[] = '\0';
					str[0..tag.length] = tag[];
					printf("Could not find tag %s\n", str.ptr);
					enforce(false, "Could not find matching boundary condition tag");
				}
				return cast(uint)tagIdx;
			}

			uint bcIdx = findBcIndex(mesh.bTags[i]);

			for(uint j = 0; j < mesh.bGroups[i].length; j++)
			{
				enforce(mesh.edges[mesh.bGroups[i][j]].isBoundary, "Edge not boundary edge but should be");
				enforce(mesh.edges[mesh.bGroups[i][j]].boundaryTag == config.boundaries[bcIdx].bTag, "Incorrect boundary tag");

				mesh.edges[mesh.bGroups[i][j]].bIdx = bcIdx;

				M = config.boundaries[bcIdx].boundaryData[0];
				aoa = config.boundaries[bcIdx].boundaryData[1] * (PI/180);
				rho = config.boundaries[bcIdx].boundaryData[3];

				U = 1.0;
				a = U/M;
				p = (rho*a^^2.0)/gamma;
				u = U*cos(aoa);
				v = U*sin(aoa);

				mesh.edges[mesh.bGroups[i][j]].boundaryType = config.boundaries[bcIdx].type;

				if(mesh.edges[mesh.bGroups[i][j]].boundaryType == BoundaryType.FullState)
				{
					mesh.edges[mesh.bGroups[i][j]].q[1][0] = rho;
					mesh.edges[mesh.bGroups[i][j]].q[1][1] = rho*u;
					mesh.edges[mesh.bGroups[i][j]].q[1][2] = rho*v;
					mesh.edges[mesh.bGroups[i][j]].q[1][3] = p/(gamma - 1.0) + 0.5*rho*(u^^2.0 + v^^2.0);
				}
				else if(mesh.edges[mesh.bGroups[i][j]].boundaryType == BoundaryType.ConstPressure)
				{
					enforce(config.boundaries[bcIdx].boundaryData.length == 1, "Constant pressure boundary data should only have one element; the constant pressure");
					mesh.edges[mesh.bGroups[i][j]].bData = config.boundaries[bcIdx].boundaryData;
				}
				else if(mesh.edges[mesh.bGroups[i][j]].boundaryType == BoundaryType.Dirichlet)
				{
					auto mid = mesh.edges[mesh.bGroups[i][j]].mid;
					mesh.edges[mesh.bGroups[i][j]].q[1] = solution(mid[0], mid[1]);
					config.boundaries[bcIdx].dFunc = &solution;
				}
			}
		}

		auto R = new Vector!4[mesh.interiorCells.length];
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
											double Rmax;
											ufvmSolver!(mixin(lim), mixin(fl), 4)(R, mesh.q, mesh, config, t, Rmax, config.limited, false);
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

	}
}

int main(string[] args)
{
	if(args.length < 3)
	{
		writeln("Not enough input arguments");
		writeln("Usage: ebb-reconstruct configFile save_#_*");
		return -1;
	}

	init(args);

	string configFile = args[1];
	string[] slnFiles = args[2..$];
	slnFiles = slnFiles.sort!`a.chomp(".esln")[a.lastIndexOf('_')+1..$].to!int < b.chomp(".esln")[b.lastIndexOf('_')+1..$].to!int`.array;

	auto saveNum = slnFiles[0].split("_")[$-2];
	writeln("saveNum: ", saveNum);
	string savePrefix = "recon_save_"~saveNum;

	import std.path : dirName, asAbsolutePath;
	auto configStr = readText(configFile);
	auto config = loadConfig(configStr);
	string meshFile = "/";
	writeln(meshFile);
	meshFile = configFile.dirName.asAbsolutePath.to!string ~ meshFile;
	writeln(meshFile);
	meshFile ~= config.meshFile;
	writeln(meshFile);
	UMesh2 mesh = importMesh(meshFile);
	mesh.buildMesh;

	double dt, t;
	
	foreach(i, slnFile; slnFiles)
	{
		enforce(loadSolution(mesh, t, dt, slnFile, true), "Failed to load solution file: "~slnFile);
	}

	immutable size_t dims = 4;
	stepMesh(mesh, config, t, dt);
	Vector!dims[] nodeVals;
	uint[] edgeNodeIdx;
	double[][] edgeNodes;
	Vector!dims[] edgeVals;
	
	uint[] centroidNodeIdx;
	double[][] centroidNodes;

	int[] newNodeMap = new int[mesh.nodes.length];
	newNodeMap[] = -1;
	double[][] newNodes;
	size_t[][] nodeCells;
	bool[] isBnode;
	foreach(i, node; mesh.nodes)
	{
		size_t[] cells;
		foreach(j, element; mesh.elements)
		{
			enforce(!element.canFind(0), "found 0 element");
			auto elIdx = element.countUntil(i+1);
			if(elIdx >= 0)
			{
				cells ~= j;
			}
		}
		if(cells.length > 0)
		{
			// Not all nodes will neccessarily be associated with an element
			// so we need to remove orphaned nodes.
			newNodes ~= node;
			newNodeMap[i] = newNodes.length.to!int;
			nodeCells ~= cells;
			isBnode ~= false;
			foreach(bNode; mesh.bNodes)
			{
				if(bNode.canFind(i))
				{
					isBnode[$-1] = true;
					break;
				}
			}
		}
	}

	// Remap element list to new node list
	foreach(i; 0..mesh.elements.length)
	{
		foreach(j; 0..mesh.elements[i].length)
		{
			enforce(newNodeMap[mesh.elements[i][j]-1] != -1, "negative mapping at i: "~i.to!string~" j: "~j.to!string);
			mesh.elements[i][j] = newNodeMap[mesh.elements[i][j]-1].to!uint;
		}
	}

	// Remap bNode list to new node list
	foreach(i; 0..mesh.bNodes.length)
	{
		foreach(j; 0..mesh.bNodes[i].length)
		{
			enforce(newNodeMap[mesh.bNodes[i][j]] != -1, "negative mapping at i: "~i.to!string~" j: "~j.to!string);
			mesh.bNodes[i][j] = newNodeMap[mesh.bNodes[i][j]].to!uint - 1;
		}
	}

	foreach(i; 0..mesh.edges.length)
	{
		edgeNodes ~= [mesh.edges[i].mid[0], mesh.edges[i].mid[1]];
		edgeVals ~= 0.5*(mesh.edges[i].q[0] + mesh.edges[i].q[1]);
	}

	writeln("mesh.bNodes.length = ", mesh.bNodes.length);
	writeln("mesh.boundaryEdges.length = ", mesh.boundaryEdges.length);
	// run through the nodes and compute node averaged values using
	// second order gradient reconstruction.
	foreach(i, node; newNodes)
	{
		auto nodeVal = Vector!dims(0);
		if(!isBnode[i])
		{
			foreach(cell; nodeCells[i])
			{
				auto qM = mesh.q[cell];
				auto grad = mesh.cells[cell].gradient;
				auto centroid = mesh.cells[cell].centroid;
				auto nodeCoords = Vector!2(node[0], node[1]);

				auto dx = nodeCoords[0] - centroid[0];
				auto dy = nodeCoords[1] - centroid[1];

				for(uint j = 0; j < dims; j++)
				{
					nodeVal[j] += qM[j] + mesh.cells[cell].lim[j][0]*grad[j,0]*dx + 
										  mesh.cells[cell].lim[j][1]*grad[j,1]*dy;
				}
			}
			nodeVal *= (1.0/nodeCells[i].length);
			immutable bool haveNan = (nodeVal[0].isNaN || nodeVal[1].isNaN || nodeVal[2].isNaN || nodeVal[3].isNaN);
			enforce(!haveNan, "New node value is NaN");
			nodeVals ~= nodeVal;
		}
		else
		{
			uint numEdges = 0;
			foreach(cell; nodeCells[i])
			{
				foreach(edgeIdx; 0..mesh.cells[cell].nEdges)
				{
					auto n1 = mesh.elements[cell][edgeIdx];
					auto n2 = mesh.elements[cell][(edgeIdx+1)%mesh.cells[cell].nEdges];
					if(mesh.edges[mesh.cells[cell].edges[edgeIdx]].isBoundary && ((n1 == i+1) || (n2 == i+1)))
					{
						nodeVal += edgeVals[mesh.cells[cell].edges[edgeIdx]];
						numEdges++;
					}
				}
				
			}
			enforce(numEdges == 2, "node is connected to more than 2 boundary edges somehow");
			nodeVal *= (1.0/numEdges);
			immutable bool haveNan = (nodeVal[0].isNaN || nodeVal[1].isNaN || nodeVal[2].isNaN || nodeVal[3].isNaN);
			enforce(!haveNan, "New boundary node value is NaN");
			nodeVals ~= nodeVal;
		}
	}

	foreach(i; 0..mesh.edges.length)
	{
		edgeNodeIdx ~= (newNodes.length + i + 1).to!uint;
	}
	enforce(edgeNodeIdx.length == mesh.edges.length, "edge node index map and edges lengths don't match");

	foreach(i; mesh.interiorCells)
	{
		centroidNodes ~= [mesh.cells[i].centroid[0], mesh.cells[i].centroid[1]];
		centroidNodeIdx ~= (newNodes.length + edgeNodes.length + i + 1).to!uint;
	}
	enforce(centroidNodeIdx.length == mesh.interiorCells.length, "centroid node index map and edges lengths don't match");

	auto reconNodes = newNodes ~ edgeNodes ~ centroidNodes;
	auto reconQ = nodeVals ~ edgeVals ~ mesh.q[0..mesh.interiorCells.length];

	enforce(reconNodes.length == reconQ.length, "different number of nodes and solution points");
	uint[][] reconElements;

	UMesh2 reconMesh;
	foreach(i, element; mesh.elements)
	{
		auto cell = mesh.cells[i];
		enforce(element.length == cell.nEdges, "Element and edge count don't match");
		foreach(j; 0..cell.nEdges)
		{
			auto n1 = element[j];
			auto n2 = element[(j + 1)%cell.nEdges];
			auto edgeIdx = cell.edges[j];
			auto en = edgeNodeIdx[edgeIdx];
			auto cn = centroidNodeIdx[i];
			enforce(![n1, en, cn].canFind(-1), "found negative element 1 indicie: "~[n1, en, cn].to!string);
			enforce(![en, n2, cn].canFind(-1), "found negative element 2 indicie: "~[n1, en, cn].to!string);
			enforce(![n1, en, cn].canFind(0), "found zero element 1 indicie: "~[n1, en, cn].to!string);
			enforce(![en, n2, cn].canFind(0), "found zero element 2 indicie: "~[n1, en, cn].to!string);
			reconElements ~= [n1, en, cn];
			reconElements ~= [en, n2, cn];
		}
	}

	reconMesh.nodes = reconNodes;
	reconMesh.elements = reconElements;
	reconMesh.q = reconQ;

	reconMesh.cells = new UCell2[reconMesh.nodes.length];
	reconMesh.interiorCells = std.range.iota(0, reconMesh.nodes.length).array.to!(uint[]);
	reconMesh.localToGlobalElementMap = std.range.iota(0, reconMesh.nodes.length).array.to!(uint[]);
	foreach(i; 0..reconMesh.cells.length)
	{
		reconMesh.cells[i].nEdges = 0;
	}

	saveMatlabMesh(reconMesh, "recon_"~config.meshFile.split(".")[0]~".mmsh");
	saveSolution(reconMesh.q, reconMesh, cast(char*)(savePrefix~".esln").toStringz, t, dt, cast(uint)config.order);

	writeln("exiting");
	return shutdown;
}
