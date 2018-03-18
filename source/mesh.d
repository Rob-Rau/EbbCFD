/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.mesh;

import std.algorithm;
import std.array;
import std.conv;
import std.math;
import std.range;
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
	size_t[2] nodeIdx;
	double len;
	Vector!2 normal;
	Vector!2 tangent;
	Vector!2 mid;

	Matrix!(2, 2) rotMat;
	// neighboring cells index
	size_t cellIdxL;
	size_t cellIdxR;

	// This is a boundary edg1e
	bool isBoundary;
	string boundaryTag;
	uint bTag;
	size_t bIdx;
	BoundaryType boundaryType;
	Vector!2 bNormal;
	// data that may be needed for this boundary
	double[] bData;
}

struct CommEdgeNodes
{
	size_t n1;
	size_t n2;
}

struct Cell
{
	size_t[MAX_EDGES] edges;
	double[MAX_EDGES] fluxMultiplier;
	size_t[MAX_EDGES] neighborCells;
	Matrix!(2,MAX_EDGES) gradMat;

	uint nNeighborCells;
	uint nEdges;
	double area = 0;
	double d = 0;
	double perim = 0;
	Vector!2 centroid;
}

@nogc double determinant(Matrix!(2, 2) mat)
{
	return mat[0]*mat[3] - mat[1]*mat[2];
}

struct Mesh
{
	Comm comm;
	uint mpiRank;
	static const int meshTag = 2000;
	// Raw mesh data
	double[][] nodes;

	// elements are 1 indexed
	size_t[][] elements;

	/// bNodes are 0 indexed
	size_t[][] bNodes;
	string[] bTags;
	size_t[] bGroupStart;
	size_t[][] bGroups;

	size_t[] commProc;
	Request[] sendRequests;
	Request[] sendGradRequests;
	Request[] recvRequests;
	Request[] recvGradRequests;
	Status[] statuses;
	CommEdgeNodes[][] commEdgeLists;

	size_t[] localToGlobalElementMap;

	size_t[][] commEdgeIdx;
	size_t[][] commCellSendIdx;
	size_t[][] commCellRecvIdx;

	// Computed mesh data
	size_t[] interiorCells;
	size_t[] nonCommCells; // interior cells that DON'T have a comm ghost neighbor
	size_t[] needCommCells; // interior cells that DO have a comm ghost neighbor
	size_t[] ghostCells;
	Cell[] cells;

	size_t[] interiorEdges;
	size_t[] boundaryEdges;
	Edge[] edges;
	Vector!2[] normals;

	this(uint nCells)
	{
		cells = new Cell[nCells];
	}

	this(uint nCells, Comm comm, uint rank)
	{
		cells = new Cell[nCells];
		this.comm = comm;
		this.mpiRank = rank;
	}

	this(Comm, uint rank)
	{
		this.comm = comm;
		this.mpiRank = rank;
	}

	alias BoundaryResult = Tuple!(bool, "isBoundary", size_t, "bGroup", size_t, "bNodeIdx");
	alias boundaryResult = tuple!("isBoundary", "bGroup", "bNodeIdx");
	private BoundaryResult isBoundaryEdge(H)(H boundaryHash, ref Edge edge)//, ref uint bGroup, ref uint bNodeIdx)
	{
		immutable auto n1 = edge.nodeIdx[0] > edge.nodeIdx[1] ? edge.nodeIdx[1] : edge.nodeIdx[0];
		immutable auto n2 = edge.nodeIdx[0] > edge.nodeIdx[1] ? edge.nodeIdx[0] : edge.nodeIdx[1];

		if(boundaryHash[n1, n2] != 0)
		{
			size_t bGroup;
			size_t bNodeIdx;
			for(uint j = 0; j < bGroupStart.length-1; j++)
			{
				if(((boundaryHash[n1, n2] - 1) >= bGroupStart[j]) && ((boundaryHash[n1, n2] - 1) < bGroupStart[j+1]))
				{
					bGroup = j;
				}
			}
			if((boundaryHash[n1, n2] - 1) >= bGroupStart[$-1])
			{
				bGroup = bGroupStart.length - 1;
			}
			bNodeIdx = boundaryHash[n1, n2] - 1;
			//return true;
			return boundaryResult(true, bGroup, bNodeIdx);
		}
		//return false;
		return boundaryResult(false, 0UL, 0UL);
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
		bGroups = new size_t[][](bTags.length);
		// Hash tables are 1 indexed so that 0 can identify an
		// edge that has not been computed yet.
		auto edgeHash = sparse!size_t(nodes.length, nodes.length);
		auto commHash = sparse!size_t(nodes.length, nodes.length);
		auto boundaryHash = sparse!uint(nodes.length, nodes.length);

		foreach(uint i, bNode; bNodes)
		{
			immutable auto n1 = bNode[0] > bNode[1] ? bNode[1] : bNode[0];
			immutable auto n2 = bNode[0] > bNode[1] ? bNode[0] : bNode[1];
			boundaryHash[n1, n2] = i + 1;
		}

		Edge[] iEdges;
		Edge[] bEdges;
		Edge[] cEdges;

		foreach(i, element; elements)
		{
			cells[i].edges[] = 0;
			cells[i].fluxMultiplier[] = 1.0;
			cells[i].area = 0.0;
			cells[i].d = 0.0;
			cells[i].perim = 0;
			cells[i].centroid = Vector!2(0);
			cells[i].nNeighborCells = 0;
			cells[i].nEdges = 0;

			foreach(j, nodeIdx; element)
			{
				Edge edge;
				edge.nodeIdx = [elements[i][j]-1, elements[i][(j + 1)%element.length]-1];
				size_t[] ni = edge.nodeIdx;

				edge.len = sqrt((nodes[ni[0]][0] - nodes[ni[1]][0])^^2 + (nodes[ni[0]][1] - nodes[ni[1]][1])^^2);
				Vector!2 normal = (1.0/edge.len)*Vector!2(nodes[ni[1]][1] - nodes[ni[0]][1], nodes[ni[0]][0] - nodes[ni[1]][0]);
				Vector!2 tangent = Vector!2(-normal[1], normal[0]);
				edge.normal = normal;
				edge.tangent = tangent;
				auto rotMat = Matrix!(2, 2)(normal[0], tangent[0], normal[1], tangent[1]).inverse;
				enforce(!rotMat.isNull, "Failed to build rotation matrix for edge "~j.to!string);
				edge.rotMat = rotMat.get;

				auto boundaryRes = isBoundaryEdge(boundaryHash, edge);
				edge.isBoundary = boundaryRes.isBoundary;
				if(!edge.isBoundary)
				{
					immutable auto n1 = ni[0] > ni[1] ? ni[1] : ni[0];
					immutable auto n2 = ni[0] > ni[1] ? ni[0] : ni[1];

					enforce(
						(edgeHash[n1, n2] == 0) || (commHash[n1, n2] == 0),
						"Edge should not be in both comm list end interior list"
					);

					if((edgeHash[n1, n2] == 0) && (commHash[n1, n2] == 0))
					{
						edge.cellIdxL = i;
						if(isCommEdge(edge))
						{
							cEdges ~= edge;
							commHash[n1, n2] = cEdges.length;
						}
						else
						{
							iEdges ~= edge;
							edgeHash[n1, n2] = iEdges.length;
						}
					}
					else if((edgeHash[n1, n2] == 0) && (commHash[n1, n2] != 0))
					{
						enforce(false, "We hit an already hashed comm edge. This should not happen");
					}
					else if((edgeHash[n1, n2] != 0) && (commHash[n1, n2] == 0))
					{
						iEdges[edgeHash[n1, n2]-1].cellIdxR = i;
					}
				}
				else
				{
					edge.cellIdxL = i;
					edge.boundaryTag = bTags[boundaryRes.bGroup];
					auto bNorm1 = nodes[bNodes[boundaryRes.bNodeIdx][1]][1] - nodes[bNodes[boundaryRes.bNodeIdx][0]][1];
					auto bNorm2 = nodes[bNodes[boundaryRes.bNodeIdx][0]][0] - nodes[bNodes[boundaryRes.bNodeIdx][1]][0];
					edge.bNormal = (1.0/edge.len)*Vector!2(bNorm1, bNorm2);
					bEdges ~= edge;
				}

				cells[i].perim += edge.len;
				auto mat = Matrix!(2, 2)(nodes[ni[0]][0], nodes[ni[1]][0], nodes[ni[0]][1], nodes[ni[1]][1]);
				cells[i].area += mat.determinant;
				cells[i].centroid[0] += (nodes[ni[0]][0] + nodes[ni[1]][0])*(nodes[ni[0]][0]*nodes[ni[1]][1] - nodes[ni[1]][0]*nodes[ni[0]][1]);
				cells[i].centroid[1] += (nodes[ni[0]][1] + nodes[ni[1]][1])*(nodes[ni[0]][0]*nodes[ni[1]][1] - nodes[ni[1]][0]*nodes[ni[0]][1]);
			}

			cells[i].area *= 0.5;
			enforce(cells[i].area > 0, "Computed cell "~i.to!string~" area negative");
			cells[i].d = (2*cells[i].area)/cells[i].perim;
			cells[i].centroid[0] *= 1/(6*cells[i].area);
			cells[i].centroid[1] *= 1/(6*cells[i].area);
		}

		edges = iEdges ~ cEdges ~ bEdges;
		interiorEdges = iota(0, iEdges.length).array;
		auto commEdges = iota(iEdges.length, iEdges.length + cEdges.length).array;
		boundaryEdges = iota(cEdges.length, cEdges.length + bEdges.length).array;

		foreach(i; interiorEdges)
		{
			auto edge = edges[i];
			auto lIdx = edge.cellIdxL;
			auto rIdx = edge.cellIdxR;
			cells[lIdx].edges[cells[lIdx].nEdges] = i;
			cells[rIdx].edges[cells[rIdx].nEdges] = i;
			// Cells on the right of the edge have flux 
			// multipliers of -1.0. Left cells have 1.0
			cells[lIdx].fluxMultiplier[cells[lIdx].nEdges] = 1.0;
			cells[rIdx].fluxMultiplier[cells[rIdx].nEdges] = -1.0;
			cells[lIdx].nEdges++;
			cells[rIdx].nEdges++;
		}

		foreach(i; commEdges)
		{
			auto edge = edges[i];
			auto lIdx = edge.cellIdxL;
			cells[lIdx].edges[cells[lIdx].nEdges] = i;
			cells[lIdx].fluxMultiplier[cells[lIdx].nEdges] = 1.0;
			cells[lIdx].nEdges++;
		}

		foreach(i; boundaryEdges)
		{
			auto edge = edges[i];
			auto lIdx = edge.cellIdxL;
			cells[lIdx].edges[cells[lIdx].nEdges] = i;
			cells[lIdx].fluxMultiplier[cells[lIdx].nEdges] = 1.0;
			cells[lIdx].nEdges++;
		}

		foreach(commEdgeList; commEdgeLists)
		{
			commEdgeIdx.length++;
			commCellRecvIdx.length++;
			commCellSendIdx.length++;

			foreach(commEdge; commEdgeList)
			{
				foreach(edgeIdx; commEdges)
				{
					auto edge = edges[edgeIdx];
					// is edge a communication edge
					if(((edge.nodeIdx[0] == commEdge.n1) && (edge.nodeIdx[1] == commEdge.n2)) ||
						((edge.nodeIdx[0] == commEdge.n2) && (edge.nodeIdx[1] == commEdge.n1)))
					{
						commEdgeIdx[$-1] ~= edgeIdx;
						Cell cell;
						cell.edges[] = 0;
						cell.fluxMultiplier[] = -1.0;
						cell.area = 0.0;
						cell.d = 0.0;
						cell.perim = 0;
						cell.neighborCells[0] = edge.cellIdxL;
						cell.centroid = Vector!2(0.0);
						cell.nNeighborCells = 1;
						cell.edges[0] = edgeIdx;

						commCellSendIdx[$-1] ~= cell.neighborCells[0];
						commCellRecvIdx[$-1] ~= cells.length;
						cells ~= cell;
						edges[edgeIdx].cellIdxR = commCellRecvIdx[$-1][$-1];
					}
				}
			}
		}

		foreach(bEdge; boundaryEdges)
		{
			Cell cell;
			cell.edges[] = 0;
			cell.fluxMultiplier[] = -1.0;
			cell.area = 0.0;
			cell.d = 0.0;
			cell.perim = 0;
			
			// reflect centroid across edge
			auto edge = edges[bEdge];
			auto otherCell = cells[edge.cellIdxL];
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
			cell.nNeighborCells = 1;
			cell.neighborCells[0] = edges[bEdge].cellIdxL;
			cell.edges[0] = bEdge;

			edges[bEdge].cellIdxR = cells.length;
			cells ~= cell; 

			ghostCells ~= edges[bEdge].cellIdxR;
		}
		assert(false);
		/+
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
		}+/
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
					assert(false);
					//auto pn = Vector!2(edges[bGroups[bgIdx][i]].flux[1], edges[bGroups[bgIdx][i]].flux[2]);
					//f += pn*edges[bGroups[bgIdx][i]].len;
				}
			}
		}

		auto globalF = Vector!2(0);

		//comm.reduce(f[], globalF[], MPI_SUM, 0);

		return globalF;
	}
}

/++
	Takes a mesh of arbitrary cell types and converts all
	cells to triangles. Purely for display purposes only
+/
Tuple!(Mesh, uint[]) triangulate(ref Mesh inMesh)
{
	Mesh tMesh;

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

	tMesh.cells = new Cell[tMesh.elements.length];
	for(uint i = 0; i < tMesh.elements.length; i++)
	{
		tMesh.cells[i].nEdges = 3;
	}
	//tMesh.q = new Vector!4[tMesh.elements.length];
	assert(false);
	tMesh.bNodes = inMesh.bNodes;
	tMesh.bGroupStart = inMesh.bGroupStart;
	tMesh.bTags = inMesh.bTags;

	//tMesh.buildMesh;

	return tuple(tMesh, triMap);
}

version(Have_mpi)
{
	idxtype[] buildPartitionMap(ref Mesh bigMesh, uint p, uint id, Comm comm)
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

		idxtype currIdx = 0;
		foreach(el; bigMesh.elements)
		{
			eptr~= currIdx;
			foreach(ind; el)
			{
				eind ~= (ind - 1).to!idxtype;
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

	size_t localElementMap(size_t elNode, ref double[][] localNodes, double[][] globalNodes)
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

Mesh partitionMesh(Communication)(ref Mesh bigMesh, uint p, uint id, Communication comm)
{
	version(Have_mpi)
	{
		int partTag = 2001;
		auto smallMesh = Mesh(comm, id);

		if(id == 0)
		{
			bigMesh.buildMesh;
			
			auto part = bigMesh.buildPartitionMap(p, id, comm);

			for(size_t i = 0; i < p; i++)
			{
				size_t nElems = 0;
				size_t[][] localElements;
				size_t[] localToGlobalElementMap;
				size_t[] nodesPerElement;
				double[][] localNodes;
				size_t[][] localbNodes;
				size_t[][] localbNodesUnsorted;
				size_t[] localbGroupStart;
				string[] localbTag;
				size_t[] commP;

				// holds the locally mapped edge nodes for comm boundaries.
				// first dim matches up with commP, aka the proc to send edge data to
				CommEdgeNodes[][] commEdgeList;
				Tuple!(size_t, size_t)[] bNodeMap;

				foreach(size_t j, pa; part)
				{
					if(pa == i)
					{
						size_t[] localEl;
						nodesPerElement ~= bigMesh.elements[j].length.to!size_t;
						foreach(size_t k, el; bigMesh.elements[j])
						{
							// map global node index to new proc local node index
							localEl ~= (el - 1).localElementMap(localNodes, bigMesh.nodes);

							// Is this edge a boundary edge
							size_t bIdx = 0;
							size_t b1 = 0;
							size_t b2 = 0;
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

								bNodeMap ~= tuple(bIdx, localbNodesUnsorted.length - 1);
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
						localbGroupStart ~= cast(size_t)localbNodes.length - 1;
					}
				}

				// compute communication edges
				foreach(size_t j, edge; bigMesh.interiorEdges.map!(a => bigMesh.edges[a]).array)
				{
					// Is edge comm boundary
					if(((part[edge.cellIdxL] == i) && (part[edge.cellIdxR] != i)) ||
						((part[edge.cellIdxR] == i) && (part[edge.cellIdxL] != i)))
					{
						// map edge nodes
						size_t n1 = edge.nodeIdx[0].localElementMap(localNodes, bigMesh.nodes) - 1;
						size_t n2 = edge.nodeIdx[1].localElementMap(localNodes, bigMesh.nodes) - 1;

						size_t partIdx = 0;
						if(part[edge.cellIdxL] != i)
						{
							partIdx = edge.cellIdxL;
						}
						else if(part[edge.cellIdxR] != i)
						{
							partIdx = edge.cellIdxR;
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

				size_t[] flatElementMap = localElements.joiner.array;
				double[] flatNodes = localNodes.joiner.array;

				if(i != 0)
				{
					comm.send!size_t(2, i, partTag);
					comm.sendArray(flatElementMap, i, partTag);
					comm.sendArray(nodesPerElement, i, partTag);
					comm.sendArray(flatNodes, i, partTag);
					comm.sendArray(localbGroupStart.to!(size_t[]), i, partTag);

					size_t[] flatBnodes = localbNodes.joiner.array;
					comm.sendArray(flatBnodes, i, partTag);

					comm.send!size_t(localbTag.length, i, partTag);
					foreach(tag; localbTag)
					{
						comm.sendArray!char(tag.to!(char[]), i, partTag);
					}

					comm.sendArray(commP, i, partTag);

					comm.send!size_t(commEdgeList.length, i, partTag);
					foreach(commEdge; commEdgeList)
					{
						comm.sendArray(commEdge, i, partTag);
					}

					comm.sendArray(localToGlobalElementMap, i, partTag);
				}
				else
				{
					smallMesh = Mesh(localElements.length.to!uint);
					smallMesh.elements = localElements;
					foreach(uint j, el; smallMesh.elements)
					{
						smallMesh.cells[j].nEdges = el.length.to!uint;
					}
					//smallMesh.q = new Vector!4[smallMesh.cells.length];
					assert(false);

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
			immutable auto dims = comm.recv!size_t(0, partTag);
			auto flatLocalElements = comm.recvArray!size_t(0, partTag);
			auto nodesPerElement = comm.recvArray!size_t(0, partTag);
			auto flatLocalNodes = comm.recvArray!double(0, partTag);
			auto bGroupStart = comm.recvArray!size_t(0, partTag);

			smallMesh.bGroupStart = bGroupStart.to!(size_t[]);

			auto flatBnodes = comm.recvArray!size_t(0, partTag);

			auto nBTags = comm.recv!size_t(0, partTag);
			for(uint i = 0; i < nBTags; i++)
			{
				smallMesh.bTags ~= comm.recvArray!(char)(0, partTag).to!string;
			}

			auto commP = comm.recvArray!size_t(0, partTag);

			auto nCommEdges = comm.recv!size_t(0, partTag);
			CommEdgeNodes[][] commEdgeList;
			for(uint i = 0; i < nCommEdges; i++)
			{
				commEdgeList ~= comm.recvArray!CommEdgeNodes(0, partTag);
			}

			smallMesh.localToGlobalElementMap ~= comm.recvArray!size_t(0, partTag);

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
				size_t[] element;
				for(uint j = 0; j < nodesPerElement[npeIdx]; j++)
				{
					element ~= flatLocalElements[i + j];
				}
				Cell cell;
				cell.nEdges = nodesPerElement[npeIdx].to!uint;
				smallMesh.cells ~= cell;

				smallMesh.elements ~= element;
				npeIdx++;
			}
		
			logln("Local elements: ", smallMesh.elements.length);

			//smallMesh.q = new Vector!4[smallMesh.cells.length];
			assert(false);

			for(uint i = 0; i < flatBnodes.length; i += dims)
			{
				size_t[] node;
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
/+
@nogc Vec buildQ(double rho, double u, double v, double p)
{
	return Vec([rho, rho*u, rho*v, p/(gamma - 1) + 0.5*rho*(u^^2.0 + v^^2.0)]);
}
+/