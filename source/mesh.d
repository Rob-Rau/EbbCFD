/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.mesh;

import std.algorithm;
import std.array;
import std.conv;
import std.math;
import std.stdio;
import std.string;
import std.typecons : Tuple, tuple;

import numd.linearalgebra.matrix;
import numd.utility;

import ebb.euler;
import ebb.exception;
import ebb.mpid;

import parmetis;

import mir.sparse;

alias Vec = Vector!4;
alias Mat = Matrix!(4, 4);

immutable uint MAX_EDGES = 6;
enum BoundaryType
{
	FullState,
	ConstPressure,
	InviscidWall,
	ViscousWall,
	Symmetry,
	TempPresInflow,
	Dirichlet
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
	uint bIdx;
	BoundaryType boundaryType;
	Vector!2 bNormal;
	// data that may be needed for this boundary
	double[] bData;
}

struct CommEdgeNodes
{
	uint n1;
	uint n2;
}

struct UCell2
{
	uint[MAX_EDGES] edges;
	double[MAX_EDGES] fluxMultiplier;
	uint[MAX_EDGES] neighborCells;
	Matrix!(2,MAX_EDGES) gradMat;
	//Vector!2[4] gradient;
	Matrix!(4,2) gradient;
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
	Comm comm;
	uint mpiRank;
	static const int meshTag = 2000;
	// Raw mesh data
	double[][] nodes;

	// elements are 1 indexed
	uint[][] elements;

	/// bNodes are 0 indexed
	uint[][] bNodes;
	string[] bTags;
	size_t[] bGroupStart;
	uint[][] bGroups;

	uint[] commProc;
	Request[] sendRequests;
	Request[] sendGradRequests;
	Request[] recvRequests;
	Request[] recvGradRequests;
	Status[] statuses;
	CommEdgeNodes[][] commEdgeLists;

	uint[] localToGlobalElementMap;

	uint[][] commEdgeIdx;
	uint[][] commCellSendIdx;
	uint[][] commCellRecvIdx;

	// Computed mesh data
	uint[] interiorCells;
	uint[] nonCommCells; // interior cells that DON'T have a comm ghost neighbor
	uint[] needCommCells; // interior cells that DO have a comm ghost neighbor
	uint[] ghostCells;
	UCell2[] cells;
	Vector!4[] q;

	uint[] interiorEdges;
	uint[] boundaryEdges;
	Edge[] edges;

	Datatype vec2Type;
	Datatype vec4Type;

	Vector!4[][] sendStateBuffers;
	Vector!4[][] recvStateBuffers;

	Matrix!(4,2)[][] sendGradBuffers;
	Matrix!(4,2)[][] recvGradBuffers;

	this(uint nCells)
	{
		vec2Type = toMPIType!(Vector!2);
		vec4Type = toMPIType!(Vector!4);
		cells = new UCell2[nCells];
	}

	this(uint nCells, Comm comm, uint rank)
	{
		vec2Type = toMPIType!(Vector!2);
		vec4Type = toMPIType!(Vector!4);
		cells = new UCell2[nCells];
		this.comm = comm;
		this.mpiRank = rank;
	}

	this(Comm, uint rank)
	{
		vec2Type = toMPIType!(Vector!2);
		vec4Type = toMPIType!(Vector!4);
		this.comm = comm;
		this.mpiRank = rank;
	}

	private bool isBoundaryEdge(H)(H boundaryHash, ref Edge edge, ref uint bGroup, ref uint bNodeIdx)
	{
		immutable auto n1 = edge.nodeIdx[0] > edge.nodeIdx[1] ? edge.nodeIdx[1] : edge.nodeIdx[0];
		immutable auto n2 = edge.nodeIdx[0] > edge.nodeIdx[1] ? edge.nodeIdx[0] : edge.nodeIdx[1];

		if(boundaryHash[n1, n2] != 0)
		{
			for(uint j = 0; j < bGroupStart.length-1; j++)
			{
				if(((boundaryHash[n1, n2] - 1) >= bGroupStart[j]) && ((boundaryHash[n1, n2] - 1) < bGroupStart[j+1]))
				{
					bGroup = j;
				}
			}
			if((boundaryHash[n1, n2] - 1) >= bGroupStart[$-1])
			{
				bGroup = bGroupStart.length.to!uint - 1;
			}
			bNodeIdx = boundaryHash[n1, n2] - 1;
			return true;
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

	/++
		Compute cell and edge interconnections
	+/
	void buildMesh()
	{
		bGroups = new uint[][](bTags.length);
		auto edgeHash = sparse!uint(nodes.length, nodes.length);
		auto boundaryHash = sparse!uint(nodes.length, nodes.length);

		foreach(uint i, bNode; bNodes)
		{
			immutable auto n1 = bNode[0] > bNode[1] ? bNode[1] : bNode[0];
			immutable auto n2 = bNode[0] > bNode[1] ? bNode[0] : bNode[1];
			boundaryHash[n1, n2] = i + 1;
		}

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
				edge.isBoundary = isBoundaryEdge(boundaryHash, edge, bGroup, bNodeIdx);
				if(!edge.isBoundary)
				{
					immutable auto n1 = ni[0] > ni[1] ? ni[1] : ni[0];
					immutable auto n2 = ni[0] > ni[1] ? ni[0] : ni[1];

					if(edgeHash[n1, n2] != 0)
					{
						auto idx = edgeHash[n1, n2].to!size_t - 1;
						edges[idx].cellIdx[1] = i;
						edgeidx = edgeHash[n1, n2] - 1;
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
						auto rotMat = Matrix!(2, 2)(normal[0], tangent[0], normal[1], tangent[1]).inverse;
						enforce(!rotMat.isNull, "Failed to build rotation matrix for edge "~j.to!string);
						edge.rotMat = rotMat.get;
						immutable double x = 0.5*(nodes[ni[0]][0] + nodes[ni[1]][0]);
						immutable double y = 0.5*(nodes[ni[0]][1] + nodes[ni[1]][1]);
						edge.mid = Vector!2(x, y);
						cells[i].fluxMultiplier[j] = 1.0;
						edges ~= edge;
						edgeidx = edges.length.to!uint - 1;
						edgeHash[n1, n2] = edgeidx + 1;
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
					auto rotMat = Matrix!(2, 2)(normal[0], tangent[0], normal[1], tangent[1]).inverse;
					enforce(!rotMat.isNull, "Failed to build rotation matrix for edge "~j.to!string);
					edge.rotMat = rotMat.get;
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
				cells[i].area += mat.determinant;
				cells[i].centroid[0] += (nodes[ni[0]][0] + nodes[ni[1]][0])*(nodes[ni[0]][0]*nodes[ni[1]][1] - nodes[ni[1]][0]*nodes[ni[0]][1]);
				cells[i].centroid[1] += (nodes[ni[0]][1] + nodes[ni[1]][1])*(nodes[ni[0]][0]*nodes[ni[1]][1] - nodes[ni[1]][0]*nodes[ni[0]][1]);
			}

			interiorCells ~= i;
			cells[i].area *= 0.5;
			enforce(cells[i].area > 0, "Computed cell "~i.to!string~" area negative");
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
			}
			cell.centroid = Vector!2(xr, yr);
			//cell.gradient[] = Vector!2(0);
			cell.gradient = Matrix!(4, 2)(0);
			cell.nNeighborCells = 1;
			cell.lim[] = Vector!2(1.0);
			cell.neighborCells[0] = edges[bEdge].cellIdx[0];
			cell.edges[0] = bEdge;

			auto newq = Vector!4(0);
			q ~= newq;

			cells ~= cell; 
			edges[bEdge].cellIdx[1] = cast(uint)(cells.length - 1);

			ghostCells ~= edges[bEdge].cellIdx[1];
		}

		foreach(commEdgeList; commEdgeLists)
		{
			commEdgeIdx.length++;
			commCellRecvIdx.length++;
			commCellSendIdx.length++;
			sendStateBuffers.length++;
			recvStateBuffers.length++;
			sendGradBuffers.length++;
			recvGradBuffers.length++;
			sendRequests.length++;
			recvRequests.length++;
			sendGradRequests.length++;
			recvGradRequests.length++;
			statuses.length++;
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

						commCellSendIdx[$-1] ~= cell.neighborCells[0];
						commCellRecvIdx[$-1] ~= cast(uint)(cells.length - 1);
						edges[edgeIdx].cellIdx[1] = commCellRecvIdx[$-1][$-1];
					}
				}

				sendStateBuffers[$-1] = new Vector!4[commCellRecvIdx[$-1].length];
				recvStateBuffers[$-1] = new Vector!4[commCellRecvIdx[$-1].length];

				sendGradBuffers[$-1] = new Matrix!(4, 2)[commCellRecvIdx[$-1].length];
				recvGradBuffers[$-1] = new Matrix!(4, 2)[commCellRecvIdx[$-1].length];
			}
		}

		for(uint commIdx = 0; commIdx < sendRequests.length; commIdx++)
		{
			comm.sendInit(sendStateBuffers[commIdx], commProc[commIdx], meshTag, sendRequests[commIdx]);
			comm.sendInit(sendGradBuffers[commIdx], commProc[commIdx], meshTag, sendGradRequests[commIdx]);
		}

		for(uint commIdx = 0; commIdx < recvRequests.length; commIdx++)
		{
			comm.recvInit(recvStateBuffers[commIdx], commProc[commIdx], meshTag, recvRequests[commIdx]);
			comm.recvInit(recvGradBuffers[commIdx], commProc[commIdx], meshTag, recvGradRequests[commIdx]);
		}

		auto tmpRecvRequests = new Request[commProc.length];
		auto tmpRecvStatuses = new Status[commProc.length];
		Vector!2[][] centroidsRcv;
		centroidsRcv.length = commProc.length;

		foreach(commIdx, commEdges; commEdgeIdx)
		{
			auto commEdgesRcv = commEdgeIdx[commIdx];
			centroidsRcv[commIdx].length = commEdgesRcv.length;
			tmpRecvRequests[commIdx] = comm.irecv(centroidsRcv[commIdx], commProc[commIdx], meshTag);
		}

		// We now need to distribute cell centroids to neighboring 
		// processors so they can set up their ghost cells
		foreach(commIdx, commEdges; commEdgeIdx)
		{
			auto centroids = new Vector!2[commEdges.length];
			foreach(i, edge; commEdges)
			{
				centroids[i] = cells[edges[edge].cellIdx[0]].centroid; 
			}
			comm.send(centroids, commProc[commIdx], meshTag);
		}

		tmpRecvRequests.waitall(tmpRecvStatuses);
		foreach(commIdx, commEdges; commEdgeIdx)
		{
			foreach(i, edge; commEdgeIdx[commIdx])
			{
				cells[edges[edge].cellIdx[1]].centroid = centroidsRcv[commIdx][i];
			}
		}

		foreach(i; interiorCells)
		{
			cells[i].gradMat = Matrix!(2, MAX_EDGES)(0);
			auto tmpMat = Matrix!(MAX_EDGES, 2)(0);
			cells[i].nNeighborCells = 0;
			// Run through the edges and find indecies
			// of cell neighbors
			for(uint j = 0; j < cells[i].nEdges; j++)
			{
				auto edge = edges[cells[i].edges[j]];

				auto v1 = edge.mid - cells[i].centroid;

				auto eDot = v1.dot(edge.normal);
				//cells[i].fluxMultiplier[j] = eDot/abs(eDot);

				if(edge.cellIdx[0] == i)
				{
					cells[i].neighborCells[cells[i].nNeighborCells] = edge.cellIdx[1];
				}
				else
				{
					cells[i].neighborCells[cells[i].nNeighborCells] = edge.cellIdx[0];
				}
				immutable uint idx = cells[i].neighborCells[cells[i].nNeighborCells];

				tmpMat[cells[i].nNeighborCells, 0] = cells[idx].centroid[0] - cells[i].centroid[0];
				tmpMat[cells[i].nNeighborCells, 1] = cells[idx].centroid[1] - cells[i].centroid[1];

				cells[i].nNeighborCells++;
			}

			bool hasCommNeighbor = false;
			for(uint j = 0; j < cells[i].nNeighborCells; j++)
			{
				uint idx = cells[i].neighborCells[j];
				foreach(commCells; commCellRecvIdx)
				{
					if(commCells.canFind(idx))
					{
						hasCommNeighbor = true;
					}
				}
			}

			if(hasCommNeighbor)
			{
				needCommCells ~= i;
			}
			else
			{
				nonCommCells ~= i;
			}

			// This cell has a boundary edge and we need to find another
			// cell close by (sharing a node with) to reconstruct the
			// cell gradient
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
			auto c = b.inverse;
			//cells[i].gradMat = b.inverse*tmpMat.transpose;
			cells[i].gradMat = c*a;
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
					//double p = getPressure(q[edges[bGroups[bgIdx][i]].cellIdx[0]]);
					//auto len = edges[bGroups[bgIdx][i]].len;
					//f += p*len*edges[bGroups[bgIdx][i]].bNormal;
					auto pn = Vector!2(edges[bGroups[bgIdx][i]].flux[1], edges[bGroups[bgIdx][i]].flux[2]);
					f += pn*edges[bGroups[bgIdx][i]].len;
				}
			}
		}

		auto globalF = Vector!2(0);

		comm.reduce(f[], globalF[], MPI_SUM, 0);

		return globalF;
	}
}

/++
	Takes a mesh of arbitrary cell types and converts all
	cells to triangles. Purely for display purposes only
+/
Tuple!(UMesh2, uint[]) triangulate(ref UMesh2 inMesh)
{
	UMesh2 tMesh;

	tMesh.nodes = inMesh.nodes;
	uint[] triMap;

	foreach(uint i, element; inMesh.elements)
	{
		if(element.length == 3)
		{
			tMesh.elements ~= element;
			triMap ~= i;
		}
		else
		{
			foreach(j; 0..(element.length - 2))
			{
				tMesh.elements ~= [element[0], element[j + 1], element[j + 2]];
				triMap ~= i;
			}
		}
	}

	tMesh.cells = new UCell2[tMesh.elements.length];
	for(uint i = 0; i < tMesh.elements.length; i++)
	{
		tMesh.cells[i].nEdges = 3;
	}
	tMesh.q = new Vector!4[tMesh.elements.length];
	tMesh.bNodes = inMesh.bNodes;
	tMesh.bGroupStart = inMesh.bGroupStart;
	tMesh.bTags = inMesh.bTags;

	//tMesh.buildMesh;

	return tuple(tMesh, triMap);
}

version(Have_mpi)
{
	idxtype[] buildPartitionMap(ref UMesh2 bigMesh, uint p, uint id, Comm comm)
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
		Comm tmpComm = comm;
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
}

UMesh2 partitionMesh(ref UMesh2 bigMesh, uint p, uint id, Comm comm)
{
	version(Have_mpi)
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
				uint[] localToGlobalElementMap;
				uint[] nodesPerElement;
				double[][] localNodes;
				uint[][] localbNodes;
				uint[][] localbNodesUnsorted;
				uint[] localbGroupStart;
				string[] localbTag;
				uint[] commP;

				// holds the locally mapped edge nodes for comm boundaries.
				// first dim matches up with commP, aka the proc to send edge data to
				CommEdgeNodes[][] commEdgeList;
				Tuple!(size_t, uint)[] bNodeMap;

				foreach(uint j, pa; part)
				{
					if(pa == i)
					{
						uint[] localEl;
						nodesPerElement ~= bigMesh.elements[j].length.to!uint;
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
						localToGlobalElementMap ~= j;
					}
				}

				auto sortedMap = bNodeMap.sort!"a[0] < b[0]";

				foreach(kvPair; sortedMap)
				{
					auto bGroup = bigMesh.bGroupStart.countUntil!"b < a"(kvPair[0]) - 1;

					localbNodes ~= localbNodesUnsorted[kvPair[1]];
					if(bGroup < 0)
					{
						bGroup = bigMesh.bGroups.length - 1;
					}

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

					comm.sendArray(localToGlobalElementMap, i, partTag);
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
					smallMesh.localToGlobalElementMap = localToGlobalElementMap;
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

			smallMesh.localToGlobalElementMap ~= comm.recvArray!uint(0, partTag);

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
			for(uint i = 0; i < flatLocalElements.length; i += nodesPerElement[npeIdx-1])
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

		comm.barrier;

		return smallMesh;
	}
	else
	{
		bigMesh.localToGlobalElementMap = std.range.iota(0, bigMesh.elements.length).array.to!(uint[]);
		return bigMesh;
	}
}

@nogc Vec buildQ(double rho, double u, double v, double p)
{
	return Vec([rho, rho*u, rho*v, p/(gamma - 1) + 0.5*rho*(u^^2.0 + v^^2.0)]);
}
