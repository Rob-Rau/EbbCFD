/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.mesh;

import std.algorithm;
import std.array;
import std.conv;
import std.exception;
import std.file;
import std.math;
import std.stdio;
import std.string;

import numd.linearalgebra.matrix;
import numd.utility;

import ebb.euler;

import mpi;
import mpi.util;

import parmetis;

alias Vec = Vector!4;
alias Mat = Matrix!(4, 4);

private MPI_Datatype toMPIType(T)()
{
	static if(!is(T : char))
	{
		static if(is(T : int) && !is(T: uint))
		{
			pragma(msg, "here int");
			return MPI_INT32_T;
		}
		else static if(is(T : uint))
		{
			pragma(msg, "here uint");
			return MPI_UINT32_T;
		}
		else static if(is(T : long) && !is(T : ulong))
		{
			pragma(msg, "here long");
			return MPI_INT64_T;
		}
		else static if(is(T : ulong))
		{
			pragma(msg, "here ulong");
			return MPI_UINT64_T;
		}
		else static if(is(T : float) && !is(T : double))
		{
			pragma(msg, "here float");
			return MPI_FLOAT;
		}
		else static if(is(T : double))
		{
			pragma(msg, "here double");
			return MPI_DOUBLE;
		}
		else
		{
			pragma(msg, "here");
			static assert(false, "type not supported");
		}
	}
	else
	{
		pragma(msg, "here char");
		return MPI_CHAR;
	}
}

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
	// Raw mesh data
	double[][] nodes;
	uint[][] elements;
	uint[][] bNodes;
	string[] bTags;
	size_t[] bGroupStart;
	uint[][] bGroups;

	// Computed mesh data
	UCell2[] cells;
	Vector!4[] q;
	Edge[] edges;

	this(uint nCells)
	{
		cells = new UCell2[nCells];
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
					edge.bNormal = (1/edge.len)*Vector!2(nodes[bNodes[bNodeIdx][1]][1] - nodes[bNodes[bNodeIdx][0]][1], nodes[bNodes[bNodeIdx][0]][0] - nodes[bNodes[bNodeIdx][1]][0]);
					double x = 0.5*(nodes[ni[0]][0] + nodes[ni[1]][0]);
					double y = 0.5*(nodes[ni[0]][1] + nodes[ni[1]][1]);
					edge.mid = Vector!2(x, y);
					bGroups[bGroup] ~= edges.length.to!uint;
					cells[i].fluxMultiplier[j] = 1.0;
					edges ~= edge;
					edgeidx = edges.length.to!uint - 1;
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

			cells[i].area *= 0.5;
			//cells[i].area.writeln;
			cells[i].d = (2*cells[i].area)/cells[i].perim;
			cells[i].centroid[0] *= 1/(6*cells[i].area);
			cells[i].centroid[1] *= 1/(6*cells[i].area);
		}

		for(uint i = 0; i < cells.length; i++)
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

				if(!edge.isBoundary)
				{
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
			}

			// This cell has a boundary edge and we need to find another
			// cell close by (sharing a node with) to reconstruct the
			// cell gradient
			// TODO: This
			if(cells[i].nEdges == 3)
			{
				if(cells[i].nNeighborCells != cells[i].nEdges)
				{
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
		return f;
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

void logln(S...)(S args)
{
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	writeln("Proc ", id, ": ", args);
}

void send(T)(T data, uint proc, int tag, MPI_Comm comm)
{
	MPI_Send(&data, 1, toMPIType!T, proc, tag, comm);
}

T recv(T)(uint from, int tag, MPI_Comm comm)
{
	T data;
	MPI_Status status;
	MPI_Recv(&data, 1, toMPIType!T, from, tag, comm, &status);
	return data;
}

void sendArray(T)(T[] data, uint proc, int tag, MPI_Comm comm)
{
	uint len = cast(uint)data.length;
	MPI_Send(&len, 1, MPI_UINT32_T, proc, tag, comm);
	MPI_Send(data.ptr, cast(uint)data.length, toMPIType!T, proc, tag, comm);
}

T[] recvArray(T)(uint from, int tag, MPI_Comm comm)
{
	uint nElems;
	T[] data;

	MPI_Status status;
	MPI_Recv(&nElems, 1, MPI_UINT32_T, from, tag, comm, &status);
	data = new T[nElems];
	
	enforce(status.MPI_ERROR == MPI_SUCCESS, "Error receiving array length. MPI Error: "~status.MPI_ERROR.to!string);
	
	MPI_Recv(data.ptr, nElems, toMPIType!T, 0, tag, comm, &status);

	return data;
}

UMesh2 partitionMesh(ref UMesh2 bigMesh, uint p, uint id, MPI_Comm comm)
{
	int partTag = 2001;
	UMesh2 smallMesh;

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
			uint[] localbGroupStart;
			string[] localbTag;

			logln(bigMesh.bGroupStart);
			foreach(uint j, pa; part)
			{
				if(pa == i)
				{
					uint[] localEl;
					nodesPerElement ~= cast(uint)bigMesh.elements[j].length; 
					foreach(uint k, el; bigMesh.elements[j])
					{
						localEl ~= (el - 1).localElementMap(localNodes, bigMesh.nodes);

						auto bIdx1 = bigMesh.bNodes.countUntil([el - 1, bigMesh.elements[j][(k + 1)%bigMesh.elements[j].length] - 1]);
						auto bIdx2 = bigMesh.bNodes.countUntil([bigMesh.elements[j][(k + 1)%bigMesh.elements[j].length] - 1, el - 1]);

						if(bIdx1 > 0)
						{
							uint b1 = localEl[$-1] - 1;
							uint b2 = (bigMesh.elements[j][(k + 1)%bigMesh.elements[j].length] - 1).localElementMap(localNodes, bigMesh.nodes) - 1;
							localbNodes ~= [b1, b2];

							auto bGroup = bigMesh.bGroupStart.countUntil!"b < a"(bIdx1);
							 
							if(bGroup < 0)
							{
								bGroup = bigMesh.bGroups.length - 1;
							}
							
							//logln("bGroup = ", bGroup);
							//logln("bIdx1 = ", bIdx1);
							if(!localbTag.canFind(bigMesh.bTags[bGroup]))
							{
								localbTag ~= bigMesh.bTags[bGroup];
								localbGroupStart ~= cast(uint)localbNodes.length - 1;
							}
						}
						if(bIdx2 > 0)
						{
							uint b2 = localEl[$-1] - 1;
							uint b1 = (bigMesh.elements[j][(k + 1)%bigMesh.elements[j].length] - 1).localElementMap(localNodes, bigMesh.nodes) - 1;
							localbNodes ~= [b1, b2];

							//auto bGroup = bigMesh.bGroupStart.countUntil!"a > b"(bIdx2) - 1;
							auto bGroup = bigMesh.bGroupStart.countUntil!"b < a"(bIdx2);

							if(bGroup < 0)
							{
								bGroup = bigMesh.bGroups.length - 1;
							}

							//logln("bGroup = ", bGroup);
							if(!localbTag.canFind(bigMesh.bTags[bGroup]))
							{
								localbTag ~= bigMesh.bTags[bGroup];
								localbGroupStart ~= cast(uint)localbNodes.length - 1;
							}
						}
						nElems++;
					}
					localElements ~= localEl;
				}
			}

			logln("localNode length = ", localNodes.length);
			uint[] flatElementMap = localElements.joiner.array;
			double[] flatNodes = localNodes.joiner.array;

			logln("localElements length = ", localElements.length);
			logln("flat localElements length = ", flatElementMap.length);

			if(i != 0)
			{
				send!uint(2, i, partTag, comm);
				sendArray(flatElementMap, i, partTag, comm);
				sendArray(nodesPerElement, i, partTag, comm);
				sendArray(flatNodes, i, partTag, comm);
				sendArray(localbGroupStart.to!(uint[]), i, partTag, comm);

				uint[] flatBnodes = localbNodes.joiner.array;
				sendArray(flatBnodes, i, partTag, comm);

				send!uint(cast(uint)localbTag.length, i, partTag, comm);
				foreach(tag; localbTag)
				{
					sendArray!char(tag.to!(char[]), i, partTag, comm);
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
				logln("Number of local elements: ", nElems);
			}
		}

		logln("Finished partitioning");
	}
	else
	{
		auto dims = recv!uint(0, partTag, comm);
		auto flatLocalElements = recvArray!uint(0, partTag, comm);
		auto nodesPerElement = recvArray!uint(0, partTag, comm);
		auto flatLocalNodes = recvArray!double(0, partTag, comm);
		auto bGroupStart = recvArray!uint(0, partTag, comm);

		smallMesh.bGroupStart = bGroupStart.to!(size_t[]);

		auto flatBnodes = recvArray!uint(0, partTag, comm);

		auto nBTags = recv!uint(0, partTag, comm);
		logln("Recieving ", nBTags, " boundary tags");
		for(uint i = 0; i < nBTags; i++)
		{
			smallMesh.bTags ~= recvArray!(char)(0, partTag, comm).to!string;
		}

		logln(dims);
		auto nNodes = flatLocalNodes.length/dims;

		logln("flatLocalElements.length = ", flatLocalElements.length);
		logln("flatLocalNodes.length = ", flatLocalNodes.length);

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

	logln("waiting");
	MPI_Barrier(comm);

/+
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
+/
	return smallMesh;
}

import core.stdc.stdio;

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

	SlnHeader header = {slnVersion: 1, dataPoints: cast(uint)mesh.cells.length, t: t, dt: dt};

	ulong totSize = SlnHeader.sizeof + header.dataPoints*4*double.sizeof + uint.sizeof + uint.sizeof;

	ubyte[] buffer = cast(ubyte[])Mallocator.instance.allocate(totSize);
	scope(exit) Mallocator.instance.deallocate(buffer);

	size_t offset = 0;
	buffer.write!uint(header.slnMagic, &offset);
	buffer.write!uint(header.slnVersion, &offset);
	buffer.write!uint(header.dataPoints, &offset);
	buffer.write!double(header.t, &offset);
	buffer.write!double(header.dt, &offset);

	for(uint i = 0; i < mesh.cells.length; i++)
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

	SlnHeader header = {slnVersion: 1, dataPoints: cast(uint)mesh.cells.length, t: t, dt: dt};

	ulong totSize = SlnHeader.sizeof + header.dataPoints*4*double.sizeof + uint.sizeof + uint.sizeof;

	ubyte[] buffer = cast(ubyte[])Mallocator.instance.allocate(totSize);
	scope(exit) Mallocator.instance.deallocate(buffer);

	size_t offset = 0;
	buffer.write!uint(header.slnMagic, &offset);
	buffer.write!uint(header.slnVersion, &offset);
	buffer.write!uint(header.dataPoints, &offset);
	buffer.write!double(header.t, &offset);
	buffer.write!double(header.dt, &offset);

	for(uint i = 0; i < mesh.cells.length; i++)
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

	if(dataPoints != mesh.q.length)
	{
		return false;
	}

	fread(buffer.ptr, 1, 8, file);
	crc.put(buffer[]);
	t = buffer.peek!double;
	fread(buffer.ptr, 1, 8, file);
	crc.put(buffer[]);
	dt = buffer.peek!double;

	for(uint i = 0; i < mesh.q.length; i++)
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
	if(nNodes != mesh.nodes.length)
	{
		throw new Exception("Mesh file has different number of nodes than solution file.");
	}

	//writeln("nNodes = ", nNodes);

	for(uint i = 0; i < nNodes; i++)
	{
		buffer.read!(double);
		buffer.read!(double);
	}

	auto nEls = buffer.read!ulong;
	if(nEls != mesh.cells.length)
	{
		throw new Exception("Mesh file has different number of cells than solution file.");
	}
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
	if(nEls != mesh.cells.length)
	{
		throw new Exception("Mesh file has different number of cells than solution file.");
	}

	for(uint i = 0; i < nCells; i++)
	{
		mesh.q[i][0] = buffer.read!(double);
		mesh.q[i][1] = buffer.read!(double);
		mesh.q[i][2] = buffer.read!(double);
		mesh.q[i][3] = buffer.read!(double);
	}

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

	enforce(nDims == 2, new Exception(nDims.to!string~" dimensional meshes not supported"));

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

		enforce((q < 4) && (q != 0), new Exception("Unsuported q"));

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
			throw new Exception("Unsuported cell type");
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

@nogc Vec buildQ(double rho, double u, double v, double p)
{
	return Vec([rho, rho*u, rho*v, p/(gamma - 1) + 0.5*rho*(u^^2.0 + v^^2.0)]);
}
