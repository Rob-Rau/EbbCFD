/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.mesh;

import std.algorithm;
import std.array;
import std.conv;
import std.exception;
import std.math;
import std.stdio;
import std.string;
import std.typecons : Tuple, tuple;

import numd.linearalgebra.matrix;
import numd.utility;

import ebb.euler;
import ebb.mpid;

import mpi;
import mpi.util;

import parmetis;

alias Vec = Vector!4;
alias Mat = Matrix!(4, 4);

enum BoundaryType
{
	FullState,
	ConstPressure,
	InviscidWall
}

struct Edge
{
	uint[2] nodeIdx;
	double len;
	double sMax;
	Vector!2 normal;
	Vector!2 tangent;
	Vector!2 mid;

	Matrix!(2, 2) rotMat;
	// neighboring cells index
	uint[2] cellIdx;

	// Edge values for 2nd order
	Vector!4[2] q;
	
	// Flux on this edge
	Vector!4 flux;

	// This is a boundary edg1e
	bool isBoundary;
	string boundaryTag;
	BoundaryType boundaryType;
	Vector!2 bNormal;
}

struct CommEdgeNodes
{
	uint n1;
	uint n2;
}

struct UCell2
{
	uint[6] edges;
	double[6] fluxMultiplier;
	uint[6] neighborCells;
	Matrix!(2, 6) gradMat;
	Vector!2[4] gradient;
	Vector!2[4] lim;
	double[4] gradErr;
	bool[4] useLP;
	Vector!4 minQ;
	Vector!4 maxQ;

	uint nNeighborCells;
	uint nEdges;
	double area = 0;
	double d = 0;
	double perim = 0;
	Vector!2 centroid;
	double dt;
}

@nogc double determinant(Matrix!(2, 2) mat)
{
	return mat[0]*mat[3] - mat[1]*mat[2];
}

struct UMesh2
{
	MPI_Comm comm;
	uint mpiRank;
	static const int meshTag = 2000;
	// Raw mesh data
	double[][] nodes;
	uint[][] elements;
	uint[][] bNodes;
	string[] bTags;
	size_t[] bGroupStart;
	uint[][] bGroups;

	uint[] commProc;
	CommEdgeNodes[][] commEdgeLists;

	uint[][] commEdgeIdx;
	uint[][] commCellIdx;

	// Computed mesh data
	uint[] interiorCells;
	uint[] ghostCells;
	UCell2[] cells;
	Vector!4[] q;

	uint[] interiorEdges;
	uint[] boundaryEdges;
	Edge[] edges;

	MPI_Datatype vec2Type;
	MPI_Datatype vec4Type;

	Vector!4[][] stateBuffers;

	this(uint nCells)
	{
		vec2Type = toMPIType!(Vector!2);
		vec4Type = toMPIType!(Vector!4);
		cells = new UCell2[nCells];
	}

	this(uint nCells, MPI_Comm comm, uint rank)
	{
		vec2Type = toMPIType!(Vector!2);
		vec4Type = toMPIType!(Vector!4);
		cells = new UCell2[nCells];
		this.comm = comm;
		this.mpiRank = rank;
	}

	this(MPI_Comm, uint rank)
	{
		vec2Type = toMPIType!(Vector!2);
		vec4Type = toMPIType!(Vector!4);
		this.comm = comm;
		this.mpiRank = rank;
	}

	ref UCell2 opIndex(size_t i)
	{
		return cells[i];
	}

	private bool isBoundaryEdge(ref Edge edge, ref uint bGroup, ref uint bNodeIdx)
	{
		for(uint i = 0; i < bNodes.length; i++)
		{
			if(((bNodes[i][0] == edge.nodeIdx[0]) && (bNodes[i][1] == edge.nodeIdx[1])) ||
			   ((bNodes[i][1] == edge.nodeIdx[0]) && (bNodes[i][0] == edge.nodeIdx[1])))
			{
				for(uint j = 0; j < bGroupStart.length-1; j++)
				{
					if((i >= bGroupStart[j]) && (i < bGroupStart[j+1]))
					{
						bGroup = j;
					}
				}
				if(i >= bGroupStart[$-1])
				{
					bGroup = bGroupStart.length.to!uint - 1;
				}
				bNodeIdx = i;
				return true;
			}
		}
		return false;
	}

	private bool isCommEdge(ref Edge edge)
	{
		foreach(commEdgeList; commEdgeLists)
		{
			auto pred = (CommEdgeNodes a, Edge b) => (((b.nodeIdx[0] == a.n1) && (b.nodeIdx[1] == a.n2)) || ((b.nodeIdx[0] == a.n2) && (b.nodeIdx[1] == a.n1)));
			if(commEdgeList.canFind!pred(edge))
			{
				return true;
			}

		}
		return false;
	}

	private bool isExistingEdge(ref Edge edge, ref UMesh2 mesh, ref uint idx)
	{
		for(uint k = 0; k < mesh.edges.length; k++)
		{
			if(((mesh.edges[k].nodeIdx[0] == edge.nodeIdx[0]) && (mesh.edges[k].nodeIdx[1] == edge.nodeIdx[1])) ||
			   ((mesh.edges[k].nodeIdx[1] == edge.nodeIdx[0]) && (mesh.edges[k].nodeIdx[0] == edge.nodeIdx[1])))
			{
				idx = k;
				return true;
			}
		}
		return false;
	}

	/++
		Compute cell and edge interconnections
	+/
	void buildMesh()
	{
		bGroups = new uint[][](bTags.length);

		for(uint i = 0; i < elements.length; i++)
		{
			q[i] = Vector!4(0);
			cells[i].edges[] = 0;
			cells[i].fluxMultiplier[] = 1.0;
			cells[i].area = 0.0;
			cells[i].d = 0.0;
			cells[i].perim = 0;
			cells[i].centroid = Vector!2(0);
			cells[i].nNeighborCells = 0;
			cells[i].lim[] = Vector!2(1.0);

			for(uint j = 0; j < cells[i].nEdges; j++)
			{
				Edge edge;
				uint edgeidx;
				edge.nodeIdx = [elements[i][j]-1, elements[i][(j + 1)%cells[i].nEdges]-1];
				uint[] ni = edge.nodeIdx;

				uint bGroup;
				uint bNodeIdx;
				edge.isBoundary = isBoundaryEdge(edge, bGroup, bNodeIdx);
				if(!edge.isBoundary)
				{
					uint k = 0;
					if(isExistingEdge(edge, this, k))
					{
						edges[k].cellIdx[1] = i;
						edgeidx = k;
						cells[i].fluxMultiplier[j] = -1.0;
					}
					else
					{
						edge.cellIdx[0] = i;
						edge.len = sqrt((nodes[ni[0]][0] - nodes[ni[1]][0])^^2.0 + (nodes[ni[0]][1] - nodes[ni[1]][1])^^2.0);
						Vector!2 normal = Vector!2(nodes[ni[1]][1] - nodes[ni[0]][1], nodes[ni[0]][0] - nodes[ni[1]][0]);
						normal = normal.normalize();
						Vector!2 tangent = Vector!2(-normal[1], normal[0]);
						tangent = tangent.normalize();
						edge.normal = normal;
						edge.tangent = tangent;
						edge.rotMat = Matrix!(2, 2)(normal[0], tangent[0], normal[1], tangent[1]).Inverse;
						double x = 0.5*(nodes[ni[0]][0] + nodes[ni[1]][0]);
						double y = 0.5*(nodes[ni[0]][1] + nodes[ni[1]][1]);
						edge.mid = Vector!2(x, y);
						cells[i].fluxMultiplier[j] = 1.0;
						edges ~= edge;
						edgeidx = edges.length.to!uint - 1;
						if(!isCommEdge(edge))
						{
							interiorEdges ~= edgeidx;
						}
					}
				}
				else
				{
					edge.cellIdx[0] = i;
					edge.len = sqrt((nodes[ni[0]][0] - nodes[ni[1]][0])^^2 + (nodes[ni[0]][1] - nodes[ni[1]][1])^^2);
					Vector!2 normal = (1.0/edge.len)*Vector!2(nodes[ni[1]][1] - nodes[ni[0]][1], nodes[ni[0]][0] - nodes[ni[1]][0]);
					Vector!2 tangent = Vector!2(-normal[1], normal[0]);
					edge.normal = normal;
					edge.tangent = tangent;
					edge.rotMat = Matrix!(2, 2)(normal[0], tangent[0], normal[1], tangent[1]).Inverse;
					edge.boundaryTag = bTags[bGroup];
					edge.bNormal = (1.0/edge.len)*Vector!2(nodes[bNodes[bNodeIdx][1]][1] - nodes[bNodes[bNodeIdx][0]][1], nodes[bNodes[bNodeIdx][0]][0] - nodes[bNodes[bNodeIdx][1]][0]);
					double x = 0.5*(nodes[ni[0]][0] + nodes[ni[1]][0]);
					double y = 0.5*(nodes[ni[0]][1] + nodes[ni[1]][1]);
					edge.mid = Vector!2(x, y);
					bGroups[bGroup] ~= edges.length.to!uint;
					cells[i].fluxMultiplier[j] = 1.0;
					edges ~= edge;
					edgeidx = edges.length.to!uint - 1;
					boundaryEdges ~= edgeidx;
				}

				cells[i].edges[j] = edgeidx;
				cells[i].perim += edges[edgeidx].len;
				auto mat = Matrix!(2, 2)(nodes[ni[0]][0], nodes[ni[1]][0], nodes[ni[0]][1], nodes[ni[1]][1]);
				//mat.writeln;
				cells[i].area += mat.determinant;
				//cells[i].area.writeln;
				cells[i].centroid[0] += (nodes[ni[0]][0] + nodes[ni[1]][0])*(nodes[ni[0]][0]*nodes[ni[1]][1] - nodes[ni[1]][0]*nodes[ni[0]][1]);
				cells[i].centroid[1] += (nodes[ni[0]][1] + nodes[ni[1]][1])*(nodes[ni[0]][0]*nodes[ni[1]][1] - nodes[ni[1]][0]*nodes[ni[0]][1]);
			}

			interiorCells ~= i;
			cells[i].area *= 0.5;
			//cells[i].area.writeln;
			cells[i].d = (2*cells[i].area)/cells[i].perim;
			cells[i].centroid[0] *= 1/(6*cells[i].area);
			cells[i].centroid[1] *= 1/(6*cells[i].area);
		}

		foreach(bEdge; boundaryEdges)
		{
			UCell2 cell;
			cell.edges[] = 0;
			cell.fluxMultiplier[] = -1.0;
			cell.area = 0.0;
			cell.d = 0.0;
			cell.perim = 0;
			
			// reflect centroid across edge
			auto edge = edges[bEdge];
			auto otherCell = cells[edge.cellIdx[0]];
			double xr = 0;
			double yr = 0;
			if(nodes[edge.nodeIdx[0]][0] != nodes[edge.nodeIdx[1]][0])
			{
				double x1 = otherCell.centroid[0];
				double y1 = otherCell.centroid[1];
				double m = (nodes[edge.nodeIdx[0]][1] - nodes[edge.nodeIdx[1]][1])/(nodes[edge.nodeIdx[0]][0] - nodes[edge.nodeIdx[1]][0]);
				
				double x = nodes[edge.nodeIdx[0]][0];
				double y = nodes[edge.nodeIdx[0]][1];
				double c = y - m*x;
				double d = (x1 + (y1 - c)*m)/(1.0 + m^^2.0);
				xr = 2.0*d - x1;
				yr = 2.0*d*m - y1 + 2.0*c;
			}
			else
			{
				// edge is vertical
				double x = nodes[edge.nodeIdx[0]][0];
				double y = nodes[edge.nodeIdx[0]][1];
				double d = otherCell.centroid[0] - x;
				double x1 = otherCell.centroid[0];
				xr = -(2.0*d - x1);
				yr = otherCell.centroid[1];
				//logln("xc = ", otherCell.centroid[0], ", x = ", nodes[edge.nodeIdx[0]][0], ", xr = ", xr);
			}
			cell.centroid = Vector!2(xr, yr);

			//logln("new cell centroid = ", cell.centroid);
			//logln("other cell centroid = ", otherCell.centroid, "\n");
			cell.nNeighborCells = 1;
			cell.lim[] = Vector!2(1.0);
			cell.neighborCells[0] = edges[bEdge].cellIdx[0];
			cell.edges[0] = bEdge;

			auto newq = Vector!4(0);
			q ~= newq;

			cells ~= cell; 
			edges[bEdge].cellIdx[1] = cast(uint)(cells.length - 1);

			//logln("Adding ghost cell. Boundary type = ", edges[bEdge].boundaryTag, "\t cell index = ", edges[bEdge].cellIdx[1]);
			ghostCells ~= edges[bEdge].cellIdx[1];
		}

		foreach(commEdgeList; commEdgeLists)
		{
			commEdgeIdx.length++;
			commCellIdx.length++;
			stateBuffers.length++;
			foreach(commEdge; commEdgeList)
			{
				foreach(uint edgeIdx, edge; edges)
				{
					// is edge a communication edge
					if(((edge.nodeIdx[0] == commEdge.n1) && (edge.nodeIdx[1] == commEdge.n2)) ||
						((edge.nodeIdx[0] == commEdge.n2) && (edge.nodeIdx[1] == commEdge.n1)))
					{
						commEdgeIdx[$-1] ~= edgeIdx;
						UCell2 cell;
						cell.edges[] = 0;
						cell.fluxMultiplier[] = -1.0;
						cell.area = 0.0;
						cell.d = 0.0;
						cell.perim = 0;
						cell.neighborCells[0] = edge.cellIdx[0];
						cell.centroid = Vector!2(0.0);
						cell.nNeighborCells = 1;
						cell.lim[] = Vector!2(1.0);
						cell.edges[0] = edgeIdx;

						auto newq = Vector!4(0);
						q ~= newq;
						cells ~= cell;

						commCellIdx[$-1] ~= cast(uint)(cells.length - 1);
						edges[edgeIdx].cellIdx[1] = commCellIdx[$-1][$-1];
					}
				}

				stateBuffers[$-1] = new Vector!4[commCellIdx[$-1].length];
			}
		}

		//MPI_Barrier(comm);
		// We now need to distribute cell centroids to neighboring 
		// processors so they can set up their ghost cells
		foreach(commIdx, commEdges; commEdgeIdx)
		{
			auto centroids = new Vector!2[commEdges.length];
			foreach(i, edge; commEdges)
			{
				centroids[i] = cells[edges[edge].cellIdx[0]].centroid; 
			}
			MPI_Send(centroids.ptr, cast(uint)centroids.length, toMPIType!(Vector!2), commProc[commIdx], meshTag, comm);
		}

		// recieve cell centroids
		for(uint c = 0; c < commProc.length; c++)
		{
			MPI_Status status;
			MPI_Probe(MPI_ANY_SOURCE, meshTag, comm, &status);

			// determine which proc this is comming from, will be different order
			// than commProc
			if(commProc.canFind(status.MPI_SOURCE))
			{
				auto commIdx = commProc.countUntil(status.MPI_SOURCE);
				auto commEdges = commEdgeIdx[commIdx];
				auto centroids = new Vector!2[commEdges.length];

				MPI_Recv(centroids.ptr, cast(int)centroids.length, toMPIType!(Vector!2), commProc[commIdx], meshTag, comm, &status);

				foreach(i, edge; commEdges)
				{
					cells[edges[edge].cellIdx[1]].centroid = centroids[i];
				}
			}
			else
			{
				logln("Unexpected source message from ", status.MPI_SOURCE);
			}
		}

		foreach(i; interiorCells)
		{
			cells[i].gradMat = Matrix!(2, 6)(0);
			auto tmpMat = Matrix!(6, 2)(0);

			// Run through the edges and find indecies
			// of cell neighbors
			for(uint j = 0; j < cells[i].nEdges; j++)
			{
				auto edge = edges[cells[i].edges[j]];

				auto v1 = edge.mid - cells[i].centroid;
				
				/+
				if(edge.isBoundary)
				{
					auto bDot = v1.dot(edge.bNormal);
					if(bDot > 0.0)
					{
						edges[cells[i].edges[j]].bNormal *= -1.0;
					}
				}
				+/
				auto eDot = v1.dot(edge.normal);
				cells[i].fluxMultiplier[j] = eDot/abs(eDot);

				if(edge.cellIdx[0] == i)
				{
					cells[i].neighborCells[cells[i].nNeighborCells] = edge.cellIdx[1];
				}
				else
				{
					cells[i].neighborCells[cells[i].nNeighborCells] = edge.cellIdx[0];
				}
				uint idx = cells[i].neighborCells[cells[i].nNeighborCells];
				tmpMat[cells[i].nNeighborCells, 0] = cells[idx].centroid[0] - cells[i].centroid[0];
				tmpMat[cells[i].nNeighborCells, 1] = cells[idx].centroid[1] - cells[i].centroid[1];

				cells[i].nNeighborCells++;
			}

			// This cell has a boundary edge and we need to find another
			// cell close by (sharing a node with) to reconstruct the
			// cell gradient
			// TODO: This
			if(cells[i].nEdges == 3)
			{
				if(cells[i].nNeighborCells != cells[i].nEdges)
				{
					logln("Not equal");
					Edge edge;
					// find non-boundary edge
					for(uint j = 0; j < cells[i].nEdges; j++)
					{
						if(!edges[cells[i].edges[j]].isBoundary)
						{
							edge = edges[cells[i].edges[j]];
							break;
						}
					}

					// find shared node
					uint node = 0;
					for(uint j = 0; j < cells[i].nEdges; j++)
					{
						if(edges[cells[i].edges[j]].isBoundary)
						{
							if((edge.nodeIdx[0] != edges[cells[i].edges[j]].nodeIdx[0]) &&
								(edge.nodeIdx[0] != edges[cells[i].edges[j]].nodeIdx[1]))
							{
								node = edge.nodeIdx[0];
							}
							else
							{
								node = edge.nodeIdx[1];
							}
							break;
						}
					}
					
					// find edges with shared node
					uint nIdx = 0;
					double r = double.infinity;
					foreach(uint j, ref e; edges)
					{
						if(!e.isBoundary)
						{
							// does it share the node
							if((e.nodeIdx[0] == node) || (e.nodeIdx[1] == node))
							{
								// is not not our edge
								if((e.cellIdx[0] != i) || (e.cellIdx[1] != i))
								{
									if(cells[i].neighborCells[].countUntil(e.cellIdx[0]) == -1)
									{
										// not a neighboring cell
										auto vec = Vector!2(cells[e.cellIdx[0]].centroid[0] - cells[i].centroid[0], cells[e.cellIdx[0]].centroid[1] - cells[i].centroid[1]);
										auto vec1 = Vector!2(tmpMat[0,0], tmpMat[0,1]);
										auto vec2 = Vector!2(tmpMat[1,0], tmpMat[1,1]);

										// find magnitude of angles of this vec and other grad vecs
										double angle1 = abs(acos(vec.dot(vec1)/(vec.magnitude*vec1.magnitude)));
										double angle2 = abs(acos(vec.dot(vec2)/(vec.magnitude*vec2.magnitude)));

										// 
										double newR = abs(angle1 - angle2);
										if(newR < r)
										{
											r = newR;
											nIdx = e.cellIdx[0];
										}
									}

									if(cells[i].neighborCells[].countUntil(e.cellIdx[1]) == -1)
									{
										// not a neighboring cell
										auto vec = Vector!2(cells[e.cellIdx[1]].centroid[0] - cells[i].centroid[0], cells[e.cellIdx[1]].centroid[1] - cells[i].centroid[1]);
										auto vec1 = Vector!2(tmpMat[0,0], tmpMat[0,1]);
										auto vec2 = Vector!2(tmpMat[1,0], tmpMat[1,1]);

										// find magnitude of angles of this vec and other grad vecs
										double angle1 = abs(acos(vec.dot(vec1)/(vec.magnitude*vec1.magnitude)));
										double angle2 = abs(acos(vec.dot(vec2)/(vec.magnitude*vec2.magnitude)));

										// 
										double newR = abs(angle1 - angle2);
										if(newR < r)
										{
											r = newR;
											nIdx = e.cellIdx[1];
										}
									}
								}
							}
						}
					}

					cells[i].neighborCells[cells[i].nNeighborCells] = nIdx;
					tmpMat[cells[i].nNeighborCells, 0] = cells[nIdx].centroid[0] - cells[i].centroid[0];
					tmpMat[cells[i].nNeighborCells, 1] = cells[nIdx].centroid[1] - cells[i].centroid[1];
					cells[i].nNeighborCells++;
				}
			}
			
			auto a = tmpMat.transpose;
			auto b = a*tmpMat;
			cells[i].gradMat = b.Inverse*tmpMat.transpose;
		}
	}

	@nogc Vector!2 computeBoundaryForces(string[] tags)
	{
		auto f = Vector!2(0);

		foreach(ref bTag; tags)
		{
			auto bgIdx = bTags.countUntil(bTag);
			if(bgIdx > -1)
			{
				for(uint i = 0; i < bGroups[bgIdx].length; i++)
				{
					//double p = getPressure(edges[bGroups[bgIdx][i]].q[0]);
					double p = getPressure(q[edges[bGroups[bgIdx][i]].cellIdx[0]]);
					auto len = edges[bGroups[bgIdx][i]].len;
					f += p*len*edges[bGroups[bgIdx][i]].bNormal;
				}
			}
		}

		auto globalF = Vector!2(0);

		MPI_Reduce(f.mData.ptr, globalF.mData.ptr, 2, MPI_DOUBLE, MPI_SUM, 0, comm);

		return globalF;
	}
}

idxtype[] buildPartitionMap(ref UMesh2 bigMesh, uint p, uint id, MPI_Comm comm)
{
	idxtype[2] elmdist = [0, bigMesh.elements.length.to!idxtype];
	idxtype[] eptr;
	idxtype[] eind;
	int wgtflag = 0; // un-weighted
	int numflag = 0; // c style indexing
	int ncon = 1; // one constraint. still not totally sure what this is.
	int ncommonnodes = 2; // 2 shared nodes between cells
	int nparts = p;
	float[] tpwgts = new float[ncon*nparts];
	float[] ubvec = new float[ncon];
	int[4] options;
	int edgecut;
	idxtype[] part;

	tpwgts[] = 1.0/cast(double)nparts;
	ubvec[] = 1.05; // recommended value from the manual
	options[] = 0; // default options

	uint currIdx = 0;
	foreach(el; bigMesh.elements)
	{
		eptr~= currIdx;
		foreach(ind; el)
		{
			eind ~= (ind - 1);
		}
		currIdx += el.length;
	}
	eptr~= currIdx;

	logln("eind.length = ", eind.length);
	logln("eptr.length = ", eptr.length);
	part = new idxtype[bigMesh.elements.length];

	logln("Partitioning mesh");
	MPI_Comm tmpComm = comm;
	comm = MPI_COMM_SELF;
	ParMETIS_V3_PartMeshKway(elmdist.ptr, eptr.ptr, eind.ptr, null, &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts, tpwgts.ptr, ubvec.ptr, options.ptr, &edgecut, part.ptr, &comm);
	comm = tmpComm;

	logln("edgecut = ", edgecut);

	return part;
}

uint localElementMap(uint elNode, ref double[][] localNodes, double[][] globalNodes)
{
	auto idx = localNodes.countUntil(globalNodes[elNode]);

	if(idx < 0)
	{
		idx = localNodes.length;
		localNodes ~= globalNodes[elNode];
	}

	return cast(uint)(idx + 1);
}

UMesh2 partitionMesh(ref UMesh2 bigMesh, uint p, uint id, MPI_Comm comm)
{
	int partTag = 2001;
	auto smallMesh = UMesh2(comm, id);

	if(id == 0)
	{
		bigMesh.buildMesh;
		
		auto part = bigMesh.buildPartitionMap(p, id, comm);

		for(uint i = 0; i < p; i++)
		{
			uint nElems = 0;
			uint[][] localElements;
			uint[] nodesPerElement;
			double[][] localNodes;
			uint[][] localbNodes;
			uint[][] localbNodesUnsorted;
			uint[] localbGroupStart;
			string[] localbTag;

			uint[][] commEdges;
			//logln(bigMesh.bGroupStart);
			uint[] commP;

			// holds the localy mapped edge nodes for comm boundaries.
			// first dim matches up with commP, aka the proc to send edge data to
			CommEdgeNodes[][] commEdgeList;
			Tuple!(size_t, uint)[] bNodeMap;

			foreach(uint j, pa; part)
			{
				if(pa == i)
				{
					uint[] localEl;
					nodesPerElement ~= cast(uint)bigMesh.elements[j].length; 
					foreach(uint k, el; bigMesh.elements[j])
					{
						// map global node index to new proc local node index
						localEl ~= (el - 1).localElementMap(localNodes, bigMesh.nodes);

						// Is this edge a boundary edge
						size_t bIdx = 0;
						uint b1 = 0;
						uint b2 = 0;
						auto bIdx1 = bigMesh.bNodes.countUntil([el - 1, bigMesh.elements[j][(k + 1)%bigMesh.elements[j].length] - 1]);
						auto bIdx2 = bigMesh.bNodes.countUntil([bigMesh.elements[j][(k + 1)%bigMesh.elements[j].length] - 1, el - 1]);

						// map global boundary node indexes to proc local boundary node indexes
						if(bIdx1 >= 0)
						{
							bIdx = bIdx1;
							b1 = localEl[$-1] - 1;
							b2 = (bigMesh.elements[j][(k + 1)%bigMesh.elements[j].length] - 1).localElementMap(localNodes, bigMesh.nodes) - 1;
						}
						else if(bIdx2 >= 0)
						{
							bIdx = bIdx2;
							b2 = localEl[$-1] - 1;
							b1 = (bigMesh.elements[j][(k + 1)%bigMesh.elements[j].length] - 1).localElementMap(localNodes, bigMesh.nodes) - 1;
						}

						if((bIdx1 >= 0) || (bIdx2 >= 0))
						{
							// If boundary edge compute local boundary group
							localbNodesUnsorted ~= [b1, b2];

							bNodeMap ~= tuple(bIdx, cast(uint)localbNodesUnsorted.length - 1);
						}
						nElems++;
					}
					localElements ~= localEl;
				}
			}

			auto sortedMap = bNodeMap.sort!"a[0] < b[0]";

			foreach(kvPair; sortedMap)
			{
				long bGroup = bigMesh.bGroupStart.countUntil!"b < a"(kvPair[0]) - 1;

				localbNodes ~= localbNodesUnsorted[kvPair[1]];
				if(bGroup < 0)
				{
					bGroup = bigMesh.bGroups.length - 1;
				}

				//logln("bIdx = ", kvPair[0], ", unsortedIdx = ", kvPair[1], ", group = ", bigMesh.bTags[bGroup]);

				if(!localbTag.canFind(bigMesh.bTags[bGroup]))
				{
					localbTag ~= bigMesh.bTags[bGroup];
					localbGroupStart ~= cast(uint)localbNodes.length - 1;
				}
			}

			// compute communication edges
			foreach(uint j, edge; bigMesh.interiorEdges.map!(a => bigMesh.edges[a]).array)
			{
				// Is edge comm boundary
				if(((part[edge.cellIdx[0]] == i) && (part[edge.cellIdx[1]] != i)) ||
					((part[edge.cellIdx[1]] == i) && (part[edge.cellIdx[0]] != i)))
				{
					// map edge nodes
					uint n1 = edge.nodeIdx[0].localElementMap(localNodes, bigMesh.nodes) - 1;
					uint n2 = edge.nodeIdx[1].localElementMap(localNodes, bigMesh.nodes) - 1;
					commEdges ~= [n1, n2];

					uint partIdx = 0;
					if(part[edge.cellIdx[0]] != i)
					{
						partIdx = edge.cellIdx[0];
					}
					else if(part[edge.cellIdx[1]] != i)
					{
						partIdx = edge.cellIdx[1];
					}

					// add comm proc to commP if not already in there.
					auto cIdx = commP.countUntil(part[partIdx]);
					if(cIdx < 0)
					{
						if(!edge.isBoundary)
						{
							commP ~= part[partIdx];
							commEdgeList.length++;
							cIdx = commP.length - 1;
						}
					}

					if(!edge.isBoundary)
					{
						commEdgeList[cIdx] ~= CommEdgeNodes(n1, n2);
					}
				}
			}

			logln("Proc ", i, " comm: ", commP);

			uint[] flatElementMap = localElements.joiner.array;
			double[] flatNodes = localNodes.joiner.array;

			//localbGroupStart ~= cast(uint)localbNodes.length;
			//localbNodes ~= commEdges;
			//localbTag ~= "Comm";

			if(i != 0)
			{
				comm.send!uint(2, i, partTag);
				comm.sendArray(flatElementMap, i, partTag);
				comm.sendArray(nodesPerElement, i, partTag);
				comm.sendArray(flatNodes, i, partTag);
				comm.sendArray(localbGroupStart.to!(uint[]), i, partTag);

				uint[] flatBnodes = localbNodes.joiner.array;
				comm.sendArray(flatBnodes, i, partTag);

				comm.send!uint(cast(uint)localbTag.length, i, partTag);
				foreach(tag; localbTag)
				{
					comm.sendArray!char(tag.to!(char[]), i, partTag);
				}

				comm.sendArray(commP, i, partTag);

				comm.send!uint(cast(uint)commEdgeList.length, i, partTag);
				foreach(commEdge; commEdgeList)
				{
					comm.sendArray(commEdge, i, partTag);
				}
			}
			else
			{
				smallMesh = UMesh2(cast(uint)localElements.length);
				smallMesh.elements = localElements;
				foreach(uint j, el; smallMesh.elements)
				{
					smallMesh.cells[j].nEdges = cast(uint)el.length;
				}
				smallMesh.q = new Vector!4[smallMesh.cells.length];

				smallMesh.nodes = localNodes[];
				smallMesh.bGroupStart = localbGroupStart.to!(size_t[]);
				smallMesh.bNodes = localbNodes[];
				smallMesh.bTags = localbTag[];
				smallMesh.commProc = commP;
				smallMesh.commEdgeLists = commEdgeList;
				logln("Local elements: ", smallMesh.elements.length);
			}

		}
		logln("Finished partitioning");
	}
	else
	{
		immutable auto dims = comm.recv!uint(0, partTag);
		auto flatLocalElements = comm.recvArray!uint(0, partTag);
		auto nodesPerElement = comm.recvArray!uint(0, partTag);
		auto flatLocalNodes = comm.recvArray!double(0, partTag);
		auto bGroupStart = comm.recvArray!uint(0, partTag);

		smallMesh.bGroupStart = bGroupStart.to!(size_t[]);

		auto flatBnodes = comm.recvArray!uint(0, partTag);

		auto nBTags = comm.recv!uint(0, partTag);
		for(uint i = 0; i < nBTags; i++)
		{
			smallMesh.bTags ~= comm.recvArray!(char)(0, partTag).to!string;
		}

		auto commP = comm.recvArray!uint(0, partTag);

		auto nCommEdges = comm.recv!uint(0, partTag);
		CommEdgeNodes[][] commEdgeList;
		for(uint i = 0; i < nCommEdges; i++)
		{
			commEdgeList ~= comm.recvArray!CommEdgeNodes(0, partTag);
		}

		smallMesh.commProc = commP;
		smallMesh.commEdgeLists = commEdgeList;

		auto nNodes = flatLocalNodes.length/dims;

		// unpack node coords
		for(uint i = 0; i < flatLocalNodes.length; i += dims)
		{
			double[] node;
			for(uint j = 0; j < dims; j++)
			{
				node ~= flatLocalNodes[i + j];
			}
			smallMesh.nodes ~= node;
		}
		
		// unpack elements
		uint npeIdx = 0;
		for(uint i = 0; i < flatLocalElements.length; i += nodesPerElement[npeIdx])
		{
			uint[] element;
			for(uint j = 0; j < nodesPerElement[npeIdx]; j++)
			{
				element ~= flatLocalElements[i + j];
			}
			UCell2 cell;
			cell.nEdges = nodesPerElement[npeIdx];
			smallMesh.cells ~= cell;

			smallMesh.elements ~= element;
			npeIdx++;
			if(npeIdx >= nodesPerElement.length)
			{
				npeIdx = cast(uint)nodesPerElement.length - 1;
			}
		}
	
		logln("Local elements: ", smallMesh.elements.length);

		smallMesh.q = new Vector!4[smallMesh.cells.length];

		for(uint i = 0; i < flatBnodes.length; i += dims)
		{
			uint[] node;
			for(uint j = 0; j < dims; j++)
			{
				node ~= flatBnodes[i + j];
			}
			smallMesh.bNodes ~= node;
		}
	}

	MPI_Barrier(comm);

	return smallMesh;
}

@nogc Vec buildQ(double rho, double u, double v, double p)
{
	return Vec([rho, rho*u, rho*v, p/(gamma - 1) + 0.5*rho*(u^^2.0 + v^^2.0)]);
}
