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

alias Vec = Vector!4;
alias Mat = Matrix!(4, 4);

enum CellType
{
	Normal,
	Solid,
	GhostConst,
	GhostNoGradXL,
	GhostNoGradXR,
	GhostNoGradYT,
	GhostNoGradYB,
	GhostMirrorXL,
	GhostMirrorXR,
	GhostMirrorYT,
	GhostMirrorYB,
	GhostConstPressureXL,
	GhostConstPressureXR,
	GhostConstPressureYT,
	GhostConstPressureYB
}

struct Cell
{
	// x position of cell center
	double x;
	// y position of cell center
	double y;
	// x width of current cell
	double dx;
	// y width of current cell
	double dy;
	// conserved variables
	Vec q;
	Vec qL;
	Vec qR;
	Vec qT;
	Vec qB;
	
	// flux in the x direction
	Vec xFlux;
	// flux in the y direction
	Vec yFlux;
	// x direction slopes for 2nd order.
	Vec xSp;
	Vec xSm;
	// y direction slopes for 2nd order.
	Vec ySp;
	Vec ySm;
	// Type of cell.
	CellType cellType;

	bool corner;
	CellType cornerType;
}

struct Mesh
{
	Cell[] cells;
	ulong N;
	ulong M;
	this(size_t N, size_t M)
	{
		this.N = N;
		this.M = M;
		cells = new Cell[N*M];
		/+
		foreach(ref cell; cells)+/
		for(int i = 0; i < N*M; i++)
			with(CellType)
		{
			cells[i].cellType = Solid;
			cells[i].q = double.nan;
			cells[i].qL = double.nan;
			cells[i].qR = double.nan;
			cells[i].qT = double.nan;
			cells[i].qB = double.nan;
		}
	}
	
	void updateGhosts()
	{
		import std.conv : to;
		for(int i = 0; i < N; i++)
		{
			for(int j = 0; j < M; j++)
			{
				// actual array index
				size_t index = M*i + j;
				
				size_t indexT = index + 1;
				if(j == M-1)
				{
					indexT = index;
				}
				size_t indexB = 0;
				if(j > 0)
				{
					indexB = index - 1;
				}
				size_t indexL = 0;
				if(i > 0)
				{
					indexL = index - M;
				}
				size_t indexR = index + M;
				if(i == N - 1)
				{
					indexR = index;
				}
				
				size_t fixupCorner(size_t curCell, Vec cornerCell1)
				{
					size_t cornerCell2;
					switch(cells[curCell].cornerType)
						with(CellType)
					{
						case GhostMirrorXR:
						 	cornerCell2 = curCell + M;
							//writeln("GMXR fixup");
							for(int i; i < cells[curCell].q.rows; i++)
							{
								if(i == 1)
								{
									//cells[curCell].q[i] = (-cells[cornerCell1].q[i] + cells[cornerCell2].q[i])/2.0;
									cells[curCell].q[i] = (-cells[cornerCell2].qL[i]/cells[cornerCell2].qL[0])*cells[curCell].q[0];
									//cells[curCell].q[i] = (cornerCell1[i] - cells[cornerCell2].qL[0])/2.0;
									//cells[curCell].q[i] = -cells[cornerCell2].qL[i];
								}
								else if(i == 2)
								{
									cells[curCell].q[i] = (-cornerCell1[i]/cornerCell1[0])*cells[curCell].q[0];
									//cells[curCell].q[i] = (-cornerCell1[i] + cells[cornerCell2].qL[0])/2.0;
									//cells[curCell].q[i] = -cornerCell1[i];
								}
								else if(i != 2)
								{
									//cells[curCell].q[i] = (cells[cornerCell1].q[i] + cells[cornerCell2].q[i])/2.0;
									cells[curCell].q[i] = (cornerCell1[i] + cells[cornerCell2].qL[i])/2.0; 
								}
							}
							cells[curCell].qR = cells[cornerCell2].q;
							break;
							
						case GhostMirrorXL:
							cornerCell2 = curCell - M;
							for(int i; i < cells[curCell].q.rows; i++)
							{
								if(i == 1)
								{
									//cells[curCell].q[i] = (-cells[cornerCell1].q[i] + cells[cornerCell2].q[i])/2.0;
									cells[curCell].q[i] = (-cells[cornerCell2].qR[i]/cells[cornerCell2].qR[0])*cells[curCell].q[0];
									//cells[curCell].q[i] = -cells[cornerCell2].q[i];
								}
								else if(i == 2)
								{
									cells[curCell].q[i] = (-cornerCell1[i]/cornerCell1[0])*cells[curCell].q[0];
									//cells[curCell].q[i] = -cornerCell1[i];
								}
								else if(i != 2)
								{
									//cells[curCell].q[i] = (cells[cornerCell1].q[i] + cells[cornerCell2].q[i])/2.0;
									cells[curCell].q[i] = (cornerCell1[i] + cells[cornerCell2].q[i])/2.0; 
								}
							}
							/+
							foreach(int i, ref el; cells[curCell].q)
							{
								if(i == 2)
								{
									el = (-cornerCell1[i] + cells[cornerCell2].q[i])/2.0;
								}
								else
								{
									el = (cornerCell1[i] + cells[cornerCell2].q[i])/2.0;
								}
							}
							+/
							cells[curCell].qL = cells[curCell].q;
							break;
							
						case GhostMirrorYT:
							cornerCell2 = curCell + 1;
							foreach(int i, ref el; cells[curCell].q)
							{
								if(i == 2)
								{
									el = (-cornerCell1[i] + cells[cornerCell2].q[i])/2.0;
								}
								else
								{
									el = (cornerCell1[i] + cells[cornerCell2].q[i])/2.0;
								}
							}
							break;
						case GhostMirrorYB:
							cornerCell2 = curCell - 1;
							foreach(int i, ref el; cells[curCell].q)
							{
								if(i == 2)
								{
									el = (-cornerCell1[i] + cells[cornerCell2].q[i])/2.0;
								}
								else
								{
									el = (cornerCell1[i] + cells[cornerCell2].q[i])/2.0;
								}
							}
							break;
							
						case GhostNoGradXL:
							cornerCell2 = curCell - M;
							cells[curCell].q = (cornerCell1 + cells[cornerCell2].q)/2.0;
							break;
							
						case GhostNoGradXR:
							cornerCell2 = curCell + M;
							cells[curCell].q = (cornerCell1 + cells[cornerCell2].q)/2.0;
							break;
							
						case GhostNoGradYT:
							cornerCell2 = curCell + 1;
							cells[curCell].q = (cornerCell1 + cells[cornerCell2].q)/2.0;
							break;
							
						case GhostNoGradYB:
							cornerCell2 = curCell - 1;
							cells[curCell].q = (cornerCell1 + cells[cornerCell2].q)/2.0;
							break;
							
						default:
							break;
					}
					return cornerCell2;
				}
				
				if(cells[index].cellType == CellType.GhostMirrorXL)
				{
					// This cell is a mirror of the cell on its left. flip u velocity component.
					cells[index].q = cells[indexL].q;
					cells[index].q[1] = -cells[index].q[1];
					
					cells[index].qL = cells[indexL].qR;
					cells[index].qL[1] = -cells[index].qL[1];
					if(cells[index].corner)
					{
						fixupCorner(index, cells[indexL].qR);
					}
				}
				else if(cells[index].cellType == CellType.GhostMirrorXR)
				{
					// This cell is a mirror of the cell on its right. flip u velocity component.
					cells[index].q = cells[indexR].qL;
					cells[index].q[1] = -cells[index].q[1];
					
					cells[index].qR = cells[indexR].qL;
					cells[index].qR[1] = -cells[index].qR[1];
					
					if(cells[index].corner)
					{
						fixupCorner(index, cells[indexR].qL);
					}
				}
				else if(cells[index].cellType == CellType.GhostMirrorYT)
				{
					// This cell is a mirror of the cell on above it. flip v velocity component.
					/+
					for(int k = 0; k < 4; k++)
					{
						assert(cells[indexT].q[k].to!string != "nan");
					}
					+/
					//cells[indexT].cellType.to!string.writeln;
					cells[index].q = cells[indexT].q;
					//cells[index].q[2].writeln;
					cells[index].q[2] = -cells[index].q[2];
					//cells[index].q[2].writeln;
					
					//cells[index].qT = cells[indexT].qB;
					//cells[index].qT[2] = -cells[index].qT[2];
					if(cells[index].corner)
					{
						size_t cornerCell2 = fixupCorner(index, cells[indexT].q);
						//cells[index].qR = cells[cornerCell2].q;
						//cells[index].qR[1] = -cells[cornerCell2].q[1];
						//cells[index].qR = cells[index].q;
					}
					cells[index].qT = cells[index].q;
					//cells[index].qR = cells[index].q;
				}
				else if(cells[index].cellType == CellType.GhostMirrorYB)
				{
					/+
					cells[indexB].cellType.to!string.writeln;
					i.written;
					j.writeln;
					+/
					/+
					for(int k = 0; k < 4; k++)
					{
						assert(cells[indexB].q[k].to!string != "nan");
					}
					+/
					// This cell is a mirror of the cell on above it. flip v velocity component.
					cells[index].q = cells[indexB].qT;
					cells[index].q[2] = -cells[index].q[2];
					//cells[index].qB = cells[indexB].qT;
					//cells[index].qB[2] = -cells[index].qB[2];
					if(cells[index].corner)
					{
						size_t cornerCell2 = fixupCorner(index, cells[indexB].q);
						//cells[index].qR = cells[cornerCell2].q;
						//cells[index].qR[1] = -cells[cornerCell2].q[1];
					}
					//cells[index].q[2] = -cells[indexB].qT[2];
					cells[index].qB = cells[index].q;
					//cells[index].qR = cells[index].q;
				}
				else if(cells[index].cellType == CellType.GhostNoGradXL)
				{
					cells[index].q = cells[indexL].q;
					cells[index].qL = cells[indexL].qR;
					if(cells[index].corner)
					{
						fixupCorner(index, cells[indexL].qR);
					}
				}
				else if(cells[index].cellType == CellType.GhostNoGradXR) 
				{
					cells[index].q = cells[indexR].q;
					cells[index].qR = cells[indexR].qL;
					if(cells[index].corner)
					{
						fixupCorner(index, cells[indexR].qL);
					}
				}
				else if(cells[index].cellType == CellType.GhostNoGradYT)
				{
					cells[index].q = cells[indexT].q;
					cells[index].qT = cells[indexT].qB;
					if(cells[index].corner)
					{
						fixupCorner(index, cells[indexT].qB);
					}
				}
				else if(cells[index].cellType == CellType.GhostNoGradYB)
				{
					cells[index].q = cells[indexB].q;
					cells[index].qB = cells[indexB].qT;
					if(cells[index].corner)
					{
						fixupCorner(index, cells[indexB].qT);
					}
				}
				else if(cells[index].cellType == CellType.GhostConstPressureXL)
				{
					immutable auto pBoundary = getPressure(cells[index]);
					immutable auto pPlus = getPressure(cells[indexL]);
					immutable auto rhoPlus = cells[indexL].q[0];
					immutable auto sPlus = pPlus/(rhoPlus^^gamma);

					immutable auto rhoBoundary = (pBoundary/sPlus)^^(1/gamma);

					immutable auto cBoundary = getSoundSpeed(cells[index]);
					immutable auto cPlus = getSoundSpeed(cells[indexL]);

					// normal velocity
					immutable auto uPlus = getVelocity!0(cells[indexL]);
					// tangential velocity
					immutable auto vPlus = getVelocity!1(cells[indexL]);

					immutable auto jPlus = uPlus + (2*cPlus)/(gamma - 1);

					immutable auto uBoundary = jPlus - (2*cBoundary)/(gamma - 1);
					immutable auto vBoundary = vPlus;

					cells[index].q[0] = rhoBoundary;
					cells[index].q[1] = uBoundary*rhoBoundary;
					cells[index].q[2] = vBoundary*rhoBoundary;
					cells[index].q[3] = pBoundary/(gamma - 1) + 0.5*rhoBoundary*(vBoundary^^2 + uBoundary^^2);

					cells[index].qL = cells[index].q;

					//writeln(pBoundary);
				}
				else if(cells[index].cellType == CellType.GhostConstPressureXR)
				{
					immutable auto pBoundary = getPressure(cells[index]);
					immutable auto pPlus = getPressure(cells[indexR]);
					immutable auto rhoPlus = cells[indexR].q[0];
					immutable auto sPlus = pPlus/(rhoPlus^^gamma);

					immutable auto rhoBoundary = (pBoundary/sPlus)^^(1/gamma);

					immutable auto cBoundary = getSoundSpeed(cells[index]);
					immutable auto cPlus = getSoundSpeed(cells[indexR]);

					// normal velocity
					immutable auto uPlus = getVelocity!0(cells[indexR]);
					// tangential velocity
					immutable auto vPlus = getVelocity!1(cells[indexR]);

					immutable auto jPlus = uPlus + (2*cPlus)/(gamma - 1);

					immutable auto uBoundary = jPlus - (2*cBoundary)/(gamma - 1);
					immutable auto vBoundary = vPlus;

					cells[index].q[0] = rhoBoundary;
					cells[index].q[1] = uBoundary*rhoBoundary;
					cells[index].q[2] = vBoundary*rhoBoundary;
					cells[index].q[3] = pBoundary/(gamma - 1) + 0.5*rhoBoundary*(vBoundary^^2 + uBoundary^^2);

					cells[index].qR = cells[index].q;

					//writeln(pBoundary);
				}
				else if(cells[index].cellType == CellType.GhostConstPressureYT)
				{
					immutable auto pBoundary = getPressure(cells[index]);
					immutable auto pPlus = getPressure(cells[indexT]);
					immutable auto rhoPlus = cells[indexT].q[0];
					immutable auto sPlus = pPlus/(rhoPlus^^gamma);

					immutable auto rhoBoundary = (pBoundary/sPlus)^^(1/gamma);

					immutable auto cBoundary = getSoundSpeed(cells[index]);
					immutable auto cPlus = getSoundSpeed(cells[indexT]);

					// normal velocity
					immutable auto uPlus = getVelocity!1(cells[indexT]);
					// tangential velocity
					immutable auto vPlus = getVelocity!0(cells[indexT]);

					immutable auto jPlus = uPlus + (2*cPlus)/(gamma - 1);

					immutable auto uBoundary = jPlus - (2*cBoundary)/(gamma - 1);
					immutable auto vBoundary = vPlus;

					cells[index].q[0] = rhoBoundary;
					cells[index].q[2] = uBoundary*rhoBoundary;
					cells[index].q[1] = vBoundary*rhoBoundary;
					cells[index].q[3] = pBoundary/(gamma - 1) + 0.5*rhoBoundary*(vBoundary^^2 + uBoundary^^2);

					cells[index].qT = cells[index].q;

					//writeln(pBoundary);
				}
				else if(cells[index].cellType == CellType.GhostConstPressureYB)
				{
					immutable auto pBoundary = getPressure(cells[index]);
					immutable auto pPlus = getPressure(cells[indexB]);
					immutable auto rhoPlus = cells[indexB].q[0];
					immutable auto sPlus = pPlus/(rhoPlus^^gamma);

					immutable auto rhoBoundary = (pBoundary/sPlus)^^(1/gamma);

					immutable auto cBoundary = getSoundSpeed(cells[index]);
					immutable auto cPlus = getSoundSpeed(cells[indexB]);

					// normal velocity
					immutable auto uPlus = getVelocity!1(cells[indexB]);
					// tangential velocity
					immutable auto vPlus = getVelocity!0(cells[indexB]);

					immutable auto jPlus = uPlus + (2*cPlus)/(gamma - 1);

					immutable auto uBoundary = jPlus - (2*cBoundary)/(gamma - 1);
					immutable auto vBoundary = vPlus;

					cells[index].q[0] = rhoBoundary;
					cells[index].q[2] = uBoundary*rhoBoundary;
					cells[index].q[1] = vBoundary*rhoBoundary;
					cells[index].q[3] = pBoundary/(gamma - 1) + 0.5*rhoBoundary*(vBoundary^^2 + uBoundary^^2);

					cells[index].qB = cells[index].q;

					//writeln(pBoundary);
				}
				/+
				if(cells[indexL].cellType == CellType.GhostMirrorYT)
				{
					
				}
				else
				{+/
					cells[index].xSm = (cells[indexL].q - cells[index].q)/(cells[indexL].x - cells[index].x);
					cells[index].xSp = (cells[index].q - cells[indexR].q)/(cells[index].x - cells[indexR].x);
					
					cells[index].ySm = (cells[indexB].q - cells[index].q)/(cells[indexB].y - cells[index].y);
					cells[index].ySp = (cells[index].q - cells[indexT].q)/(cells[index].y - cells[indexT].y);
				//}
				
				if(cells[index].cellType != CellType.Solid)
				{
					for(int k = 0; k < 4; k++)
					{
						assert(cells[index].q[k].to!string != "nan", cells[index].cellType.to!string~" "~cells[indexL].q.to!string~" "~cells[indexR].q.to!string~" "~cells[indexT].q.to!string~" "~cells[indexB].q.to!string~" "~i.to!string~" "~j.to!string);
						assert(cells[index].qT[k].to!string != "nan", cells[index].cellType.to!string~" "~cells[indexL].qT.to!string~" "~cells[indexR].qT.to!string~" "~cells[indexT].qT.to!string~" "~cells[indexB].qT.to!string~" "~i.to!string~" "~j.to!string);
						assert(cells[index].qB[k].to!string != "nan", cells[index].cellType.to!string~" "~cells[indexL].qB.to!string~" "~cells[indexR].qB.to!string~" "~cells[indexT].qB.to!string~" "~cells[indexB].qB.to!string~" "~i.to!string~" "~j.to!string);
						assert(cells[index].qL[k].to!string != "nan", cells[index].cellType.to!string~" "~cells[indexL].qL.to!string~" "~cells[indexR].qL.to!string~" "~cells[indexT].qL.to!string~" "~cells[indexB].qL.to!string~" "~i.to!string~" "~j.to!string);
						assert(cells[index].qR[k].to!string != "nan", cells[index].cellType.to!string~" "~cells[indexL].qR.to!string~" "~cells[indexR].qR.to!string~" "~cells[indexT].qR.to!string~" "~cells[indexB].qR.to!string~" "~i.to!string~" "~j.to!string);
					}
				}
				
			}
		}
	}

	ref Cell opIndex(size_t i, size_t j)
	{
		assert(i < N);
		assert(j < M);
		// actual array index
		size_t index = M*i + j;
		
		return cells[index];
	}

	Vector!2 computeBodyForces()
	{
		auto bodyForces = Vector!2(0, 0);

		for(uint i = 0; i < N; i++)
			with(CellType)
		{
			for(uint j = 0; j < M; j++)
			{
				if(this[i,j].cellType == GhostMirrorXL)
				{
					auto p = getPressure(this[i-1,j]);
					bodyForces += p*this[i,j].dy*Vector!2(-1, 0);
				}
				else if(this[i,j].cellType == GhostMirrorXR)
				{
					auto p = getPressure(this[i+1,j]);
					bodyForces += p*this[i,j].dy*Vector!2(1, 0);
				}
				else if(this[i,j].cellType == GhostMirrorYT)
				{
					auto p = getPressure(this[i,j+1]);
					bodyForces += p*this[i,j].dx*Vector!2(0, 1);
				}
				else if(this[i,j].cellType == GhostMirrorYB)
				{
					auto p = getPressure(this[i,j-1]);
					bodyForces += p*this[i,j].dx*Vector!2(0, -1);
				}
			}
		}
		return bodyForces;
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
	// Cell everaged values
	Vector!4 q;

	uint[6] edges;
	double[6] fluxMultiplier;
	uint[6] neighborCells;
	Matrix!(2, 6) gradMat;
	Vector!2[4] gradient;
	
	Vector!4 minQ;
	Vector!4 maxQ;

	uint nNeighborCells;
	uint nEdges;
	double area = 0;
	double d = 0;
	double perim = 0;
	Vector!2 centroid;
}

double determinant(Matrix!(2, 2) mat)
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
			cells[i].q = Vector!4(0);
			cells[i].edges[] = 0;
			cells[i].fluxMultiplier[] = 1.0;
			cells[i].area = 0.0;
			cells[i].d = 0.0;
			cells[i].perim = 0;
			cells[i].centroid = Vector!2(0);
			cells[i].nNeighborCells = 0;
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
						edge.len = sqrt((nodes[ni[0]][0] - nodes[ni[1]][0])^^2 + (nodes[ni[0]][1] - nodes[ni[1]][1])^^2);
						Vector!2 normal = (1/edge.len)*Vector!2(nodes[ni[1]][1] - nodes[ni[0]][1], nodes[ni[0]][0] - nodes[ni[1]][0]);
						Vector!2 tangent = Vector!2(-normal[1], normal[0]);
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
					Vector!2 normal = (1/edge.len)*Vector!2(nodes[ni[1]][1] - nodes[ni[0]][1], nodes[ni[0]][0] - nodes[ni[1]][0]);
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

	@nogc Vector!2 computeBoundaryForces(string bTag)
	{
		auto f = Vector!2(0);

		auto bgIdx = bTags.countUntil(bTag);
		if(bgIdx > -1)
		{
			for(uint i = 0; i < bGroups[bgIdx].length; i++)
			{
				double p = getPressure(cells[edges[bGroups[bgIdx][i]].cellIdx[0]].q);
				auto len = edges[bGroups[bgIdx][i]].len;
				f += p*len*edges[bGroups[bgIdx][i]].bNormal;
			}
		}

		return f;
	}
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
		mesh.cells[i].q[0] = buffer.read!(double);
		mesh.cells[i].q[1] = buffer.read!(double);
		mesh.cells[i].q[2] = buffer.read!(double);
		mesh.cells[i].q[3] = buffer.read!(double);
	}

}

import core.stdc.stdio;
@nogc void saveMatlabSolution(ref UMesh2 mesh, char* filename)
{
	import std.experimental.allocator.mallocator : Mallocator;
	import std.bitmanip : write;

	ulong nodesSize = cast(ulong)(2*mesh.nodes.length*double.sizeof);
	ulong e2nSize = cast(ulong)(3*mesh.elements.length*double.sizeof);
	ulong slnSize = cast(ulong)(4*mesh.cells.length*double.sizeof);
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
	immutable ulong slnHeaderSize = ulong.sizeof;

	ulong totSize = nodeHeaderSize + e2nHeaderSize + ieHeaderSize + beHeaderSize + tagHeaderSize + slnHeaderSize;
	totSize += nodesSize + e2nSize + ieSize + beSize + bNameSize + slnSize;
	ubyte[] buffer = cast(ubyte[])Mallocator.instance.allocate(totSize);

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

	buffer.write!ulong(mesh.cells.length, &offset);
	for(uint i = 0; i < mesh.cells.length; i++)
	{
		buffer.write!double(mesh.cells[i].q[0], &offset);
		buffer.write!double(mesh.cells[i].q[1], &offset);
		buffer.write!double(mesh.cells[i].q[2], &offset);
		buffer.write!double(mesh.cells[i].q[3], &offset);
	}

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

char[][] readCleanLine(ref File file)
{
	return file.readln.strip.chomp.detab.split(' ').filter!(a => a != "").array;
}

UMesh2 parseXflowMesh(string meshFile)
{
	UMesh2 mesh;
	auto file = File(meshFile);

	auto headerLine = file.readCleanLine;
	immutable uint nNodes = headerLine[0].to!uint;
	immutable uint nElems = headerLine[1].to!uint;
	immutable uint nDims = headerLine[2].to!uint;

	enforce(nDims == 2, new Exception(nDims.to!string~" dimensional meshes not supported"));
/+
	writeln("Importing XFlow mesh "~meshFile);
	writeln("    nNodes = ", nNodes);
	writeln("    nElems = ", nElems);
	writeln("    nDims = ", nDims);
+/
	for(uint i = 0; i < nNodes; i++)
	{
		mesh.nodes ~= file.readCleanLine.to!(double[]);
	}

	immutable uint nBoundaryGroups = file.readCleanLine[0].to!uint;
	//writeln("    nBoundaryGroups = ", nBoundaryGroups);

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

		//writeln("        Boundary group ", i, ": faces = ", bFaces, ", nodes per face = ", nodesPerFace, ", tag = ", faceTag);

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

		//writeln("    Element group ", eGroup, ": faces = ", faces, ", q = ", q, ", subElements = ", subElements);

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
	for(uint i = 0; i < nElems; i++)
	{
		mesh.cells[i].nEdges = faces;
	}
	mesh.elements = elements;
	mesh.bNodes = bNodes;
	mesh.bGroupStart = bGroupStart;
	mesh.bTags = bTags;
	mesh.buildMesh;

	return mesh;
}

void saveMesh(ref Mesh mesh, string filename, double dt, double t)
{
	import std.bitmanip : write;
	import std.conv : to;

	size_t bufferSize = mesh.N*mesh.M*Cell.sizeof + 2*double.sizeof + 2*ulong.sizeof;
	//writeln("Requesting buffer of "~bufferSize.to!string~" bytes");
	ubyte[] buffer = new ubyte[bufferSize];
	size_t offset = 0;
	buffer.write!ulong(mesh.N, &offset);
	buffer.write!ulong(mesh.M, &offset);
	buffer.write!double(dt, &offset);
	buffer.write!double(t, &offset);
	for(int i = 0; i < mesh.N; i++)
	{
		for(int j = 0; j < mesh.M; j++)
		{
			buffer.write!double(mesh[i,j].x, &offset);
			buffer.write!double(mesh[i,j].y, &offset);
			buffer.write!double(mesh[i,j].dx, &offset);
			buffer.write!double(mesh[i,j].dy, &offset);
			buffer.write!double(mesh[i,j].q[0], &offset);
			buffer.write!double(mesh[i,j].q[1], &offset);
			buffer.write!double(mesh[i,j].q[2], &offset);
			buffer.write!double(mesh[i,j].q[3], &offset);
			
			buffer.write!double(mesh[i,j].qL[0], &offset);
			buffer.write!double(mesh[i,j].qL[1], &offset);
			buffer.write!double(mesh[i,j].qL[2], &offset);
			buffer.write!double(mesh[i,j].qL[3], &offset);
			
			buffer.write!double(mesh[i,j].qR[0], &offset);
			buffer.write!double(mesh[i,j].qR[1], &offset);
			buffer.write!double(mesh[i,j].qR[2], &offset);
			buffer.write!double(mesh[i,j].qR[3], &offset);
			
			buffer.write!double(mesh[i,j].qT[0], &offset);
			buffer.write!double(mesh[i,j].qT[1], &offset);
			buffer.write!double(mesh[i,j].qT[2], &offset);
			buffer.write!double(mesh[i,j].qT[3], &offset);
			
			buffer.write!double(mesh[i,j].qT[0], &offset);
			buffer.write!double(mesh[i,j].qT[1], &offset);
			buffer.write!double(mesh[i,j].qT[2], &offset);
			buffer.write!double(mesh[i,j].qT[3], &offset);
			
			buffer.write!double(mesh[i,j].qB[0], &offset);
			buffer.write!double(mesh[i,j].qB[1], &offset);
			buffer.write!double(mesh[i,j].qB[2], &offset);
			buffer.write!double(mesh[i,j].qB[3], &offset);
			
			buffer.write!double(mesh[i,j].xFlux[0], &offset);
			buffer.write!double(mesh[i,j].xFlux[1], &offset);
			buffer.write!double(mesh[i,j].xFlux[2], &offset);
			buffer.write!double(mesh[i,j].xFlux[3], &offset);
			
			buffer.write!double(mesh[i,j].yFlux[0], &offset);
			buffer.write!double(mesh[i,j].yFlux[1], &offset);
			buffer.write!double(mesh[i,j].yFlux[2], &offset);
			buffer.write!double(mesh[i,j].yFlux[3], &offset);
			
			buffer.write!double(mesh[i,j].xSp[0], &offset);
			buffer.write!double(mesh[i,j].xSp[1], &offset);
			buffer.write!double(mesh[i,j].xSp[2], &offset);
			buffer.write!double(mesh[i,j].xSp[3], &offset);
			
			buffer.write!double(mesh[i,j].xSm[0], &offset);
			buffer.write!double(mesh[i,j].xSm[1], &offset);
			buffer.write!double(mesh[i,j].xSm[2], &offset);
			buffer.write!double(mesh[i,j].xSm[3], &offset);
			
			buffer.write!double(mesh[i,j].ySp[0], &offset);
			buffer.write!double(mesh[i,j].ySp[1], &offset);
			buffer.write!double(mesh[i,j].ySp[2], &offset);
			buffer.write!double(mesh[i,j].ySp[3], &offset);
			
			buffer.write!double(mesh[i,j].ySm[0], &offset);
			buffer.write!double(mesh[i,j].ySm[1], &offset);
			buffer.write!double(mesh[i,j].ySm[2], &offset);
			buffer.write!double(mesh[i,j].ySm[3], &offset);
			
			buffer.write!CellType(mesh[i,j].cellType, &offset);
			
			buffer.write!bool(mesh[i,j].corner, &offset);
			buffer.write!CellType(mesh[i,j].cornerType, &offset);
		}
	}
	
	auto file = File(filename, "wb");
	ulong writeOffset = 0;
	while(writeOffset < buffer.length)
	{
		ulong chunkSize = 1024*1024*1024;
		if(buffer.length - writeOffset < chunkSize)
		{
			chunkSize = buffer.length - writeOffset;
		}
		file.rawWrite(buffer[writeOffset..writeOffset+chunkSize]);
		writeOffset += chunkSize;
	}
	file.close;
	//std.file.write(filename, cast(void[])buffer);
}

Mesh loadMesh(string file, ref double dt, ref double t)
{
	import std.bitmanip : read;
	import std.conv : to;
	ubyte[] buffer;
	size_t offset = 0;
	buffer = cast(ubyte[])std.file.read(file);
	ulong N = buffer.read!ulong;
	ulong M = buffer.read!ulong;
	dt = buffer.read!double;
	t = buffer.read!double;
	/+
	N.writeln;
	M.writeln;
	dt.writeln;
	t.writeln;
	+/
	Mesh mesh = Mesh(N, M);
	
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < M; j++)
		{
			mesh[i,j].x = buffer.read!double;
			mesh[i,j].y = buffer.read!double;
			mesh[i,j].dx = buffer.read!double;
			mesh[i,j].dy = buffer.read!double;
			mesh[i,j].q[0] = buffer.read!double;
			mesh[i,j].q[1] = buffer.read!double;
			mesh[i,j].q[2] = buffer.read!double;
			mesh[i,j].q[3] = buffer.read!double;
			mesh[i,j].qL[0] = buffer.read!double;
			mesh[i,j].qL[1] = buffer.read!double;
			mesh[i,j].qL[2] = buffer.read!double;
			mesh[i,j].qL[3] = buffer.read!double;
			mesh[i,j].qR[0] = buffer.read!double;
			mesh[i,j].qR[1] = buffer.read!double;
			mesh[i,j].qR[2] = buffer.read!double;
			mesh[i,j].qR[3] = buffer.read!double;
			mesh[i,j].qT[0] = buffer.read!double;
			mesh[i,j].qT[1] = buffer.read!double;
			mesh[i,j].qT[2] = buffer.read!double;
			mesh[i,j].qT[3] = buffer.read!double;
			mesh[i,j].qT[0] = buffer.read!double;
			mesh[i,j].qT[1] = buffer.read!double;
			mesh[i,j].qT[2] = buffer.read!double;
			mesh[i,j].qT[3] = buffer.read!double;
			mesh[i,j].qB[0] = buffer.read!double;
			mesh[i,j].qB[1] = buffer.read!double;
			mesh[i,j].qB[2] = buffer.read!double;
			mesh[i,j].qB[3] = buffer.read!double;
			mesh[i,j].xFlux[0] = buffer.read!double;
			mesh[i,j].xFlux[1] = buffer.read!double;
			mesh[i,j].xFlux[2] = buffer.read!double;
			mesh[i,j].xFlux[3] = buffer.read!double;		
			mesh[i,j].yFlux[0] = buffer.read!double;
			mesh[i,j].yFlux[1] = buffer.read!double;
			mesh[i,j].yFlux[2] = buffer.read!double;
			mesh[i,j].yFlux[3] = buffer.read!double;
			mesh[i,j].xSp[0] = buffer.read!double;
			mesh[i,j].xSp[1] = buffer.read!double;
			mesh[i,j].xSp[2] = buffer.read!double;
			mesh[i,j].xSp[3] = buffer.read!double;
			mesh[i,j].xSm[0] = buffer.read!double;
			mesh[i,j].xSm[1] = buffer.read!double;
			mesh[i,j].xSm[2] = buffer.read!double;
			mesh[i,j].xSm[3] = buffer.read!double;
			mesh[i,j].ySp[0] = buffer.read!double;
			mesh[i,j].ySp[1] = buffer.read!double;
			mesh[i,j].ySp[2] = buffer.read!double;
			mesh[i,j].ySp[3] = buffer.read!double;
			mesh[i,j].ySm[0] = buffer.read!double;
			mesh[i,j].ySm[1] = buffer.read!double;
			mesh[i,j].ySm[2] = buffer.read!double;
			mesh[i,j].ySm[3] = buffer.read!double;
			mesh[i,j].cellType = buffer.read!CellType;
			mesh[i,j].corner = buffer.read!bool;
			mesh[i,j].cornerType = buffer.read!CellType;
		}
	}
	return mesh;
}

void printMesh(ref Mesh mesh)
{
	import std.conv : to;
	for(int j = cast(int)mesh.M-1; j >= 0; j--)
	{
		write(j.to!string~"\t|");
		for(int i = 0; i < mesh.N; i++)
		{
			string cellType = "";
			switch(mesh[i,j].cellType)
				with(CellType)
			{
				case Solid:
					cellType = " S ";
					break;
				case Normal:
					cellType = " N ";
					break;
				case GhostConst:
					cellType = " GC";
					break;
				case GhostMirrorXL:
					cellType = "GML";
					break;
				case GhostMirrorXR:
					cellType = "GMR";
					break;
				case GhostMirrorYB:
					cellType = "GMB";
					break;
				case GhostMirrorYT:
					cellType = "GMT";
					break;
				case GhostNoGradXL:
					cellType = "GNL";
					break;
				case GhostNoGradXR:
					cellType = "GNR";
					break;
				case GhostNoGradYB:
					cellType = "GNB";
					break;
				case GhostNoGradYT:
					cellType = "GNT";
					break;
				case GhostConstPressureXL:
					cellType = "GPL";
					break;
				case GhostConstPressureXR:
					cellType = "GPR";
					break;
				case GhostConstPressureYB:
					cellType = "GPB";
					break;
				case GhostConstPressureYT:
					cellType = "GPT";
					break;
				default:
					cellType = "   ";
					break;
			}
			write(cellType~"|");
		}
		writeln;
	}

	write(" \t|");
	for(int i = 0; i < mesh.N; i++)
	{
		import std.format : format;
		string cellType = format("%3d", i);
		write(cellType~"|");
	}
	writeln;
}

void printMesh(ref Mesh mesh, string file)
{
	import std.conv : to;
	string output;
	for(int j = cast(int)mesh.M-1; j >= 0; j--)
	{
		//write(j.to!string~"\t|");
		output ~= j.to!string~"\t|";
		for(int i = 0; i < mesh.N; i++)
		{
			string cellType = "";
			switch(mesh[i,j].cellType)
				with(CellType)
			{
				case Solid:
					cellType = " S ";
					break;
				case Normal:
					cellType = " N ";
					break;
				case GhostConst:
					cellType = " GC";
					break;
				case GhostMirrorXL:
					cellType = "GML";
					break;
				case GhostMirrorXR:
					cellType = "GMR";
					break;
				case GhostMirrorYB:
					cellType = "GMB";
					break;
				case GhostMirrorYT:
					cellType = "GMT";
					break;
				case GhostNoGradXL:
					cellType = "GNL";
					break;
				case GhostNoGradXR:
					cellType = "GNR";
					break;
				case GhostNoGradYB:
					cellType = "GNB";
					break;
				case GhostNoGradYT:
					cellType = "GNT";
					break;
				default:
					cellType = "   ";
					break;
			}
			//write(cellType~"|");
			output ~= cellType~"|";
		}
		//writeln;
		output ~= '\n';
	}

	//write(" \t|");
	output ~= " \t|";
	for(int i = 0; i < mesh.N; i++)
	{
		import std.format : format;
		string cellType = format("%3d", i);
		//write(cellType~"|");
		output ~= cellType~"|";
	}
	//writeln;
	output ~= '\n';
	std.file.write(file, output);
}

@nogc Vec buildQ(double rho, double u, double v, double p)
{
	return Vec([rho, rho*u, rho*v, p/(gamma - 1) + 0.5*rho*(u^^2.0 + v^^2.0)]);
}

Meshgrid!double buildMeshgrid(ref Mesh mesh)
{
	Meshgrid!double meshgrid = {X: new double[][](mesh.M, mesh.N), Y: new double[][](mesh.M, mesh.N)};
	
	for(int i = 0; i < mesh.N; i++)
	{
		for(int j = 0; j < mesh.M; j++)
		{
			meshgrid.X[j][i] = mesh[i, j].x;
			meshgrid.Y[j][i] = mesh[i, j].y;
		}
	}
	return meshgrid;
}