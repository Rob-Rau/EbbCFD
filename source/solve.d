/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.solve;

import core.atomic;
import core.stdc.stdio : fopen, fwrite, fopen, printf, snprintf;
import core.sys.posix.signal;

import std.algorithm : canFind, countUntil, min, map, max, reduce, sum;
import std.conv;
import std.file;
import std.math;
import std.stdio;

import numd.utility;
import numd.linearalgebra.matrix;

import ebb.config;
import ebb.euler;
import ebb.exception;
import ebb.flux;
import ebb.integrators;
import ebb.limiters;
import ebb.mesh;
import ebb.io;
import ebb.mpid;

struct Solution(size_t dims)
{
	Vector!dims[] q;
	uint[] globalElementMap;
	double t;
}

struct SolverState
{

}

@nogc void invert(ref Matrix!(2, 2) inv, ref Matrix!(2, 2) mat)
{
	inv[0] = mat[3];
	inv[1] = -mat[1];
	inv[2] = -mat[2];
	inv[3] = mat[0];
	assert(mat.determinant != 0.0);
	inv *= (1.0/mat.determinant);
}

@nogc void LP(uint m)(ref Matrix!(m, 2) A, ref Vector!m b, ref Vector!2 c, ref Vector!2 xk)
{
	uint[2] Wk = [0, 1];
	auto AkInv = Matrix!(2, 2)(A[Wk[0],0], A[Wk[0],1], A[Wk[1],0], A[Wk[1],1]);
	auto Ak = Matrix!(2, 2)(A[Wk[0],0], A[Wk[0],1], A[Wk[1],0], A[Wk[1],1]);
	auto pk = Vector!2(0);
	auto lam = Vector!2(0);
	auto e = Vector!2(0);
	bool done = false;
	uint k = 0;
	uint[m - 2] D;
	double[m - 2] gamma;
	uint Dlen = 0;

	@nogc void invertTrans(ref Matrix!(2, 2) inv, ref Matrix!(2, 2) mat)
	{
		inv[0] = mat[3];
		inv[1] = -mat[2];
		inv[2] = -mat[1];
		inv[3] = mat[0];
		assert(mat.determinant != 0.0);
		inv *= (1.0/mat.determinant);
	}

	while(!done)
	{
		// compute lagrange multipliers
		// lam = (Ak')^-1*c
		invertTrans(AkInv, Ak);
		lam = AkInv*c;

		// If all lam >= 0, done = true
		if((lam[0] >= -10e-9) && (lam[1] >= -10e-9))
		{
			done = true;
			break;
		}

		uint q;
		if(lam[0] < lam[1])
		{
			q = Wk[0];
		}
		else
		{
			q = Wk[1];
		}

		// compute step direction
		if(lam[0] < lam[1])
		{
			e[0] = 1;
			e[1] = 0;
		}
		else
		{
			e[0] = 0;
			e[1] = 1;
		}

		// pk = (Ak)^-1*e_q
		invert(AkInv, Ak);
		pk = AkInv*e;

		// find set of decreasing constraints
		// D = (a_i'*pk < 0) note: only use i's not in Wk
		Dlen = 0;
		for(uint i = 0; i < m; i++)
		{
			if((i != Wk[0]) && (i != Wk[1]))
			{
				if(A[i, 0]*pk[0] + A[i,1]*pk[1] < -1e-11)
				{
					D[Dlen] = i;
					Dlen++;
				}
			}
		}

		// ensure D is not the null set. (shouldn't ever happen here)
		assert(abs(D[].sum) > 1.0e-11);

		// for all i in D
		// 		gamma_i = (a_i'*xk - b_i)/(-a_i'*pk)
		// alpha = min(gamma_i)
		double alpha = double.infinity;
		for(uint i = 0; i < Dlen; i++)
		{
			immutable double aixk = A[D[i],0]*xk[0] + A[D[i],1]*xk[1];
			immutable double aipk = A[D[i],0]*pk[0] + A[D[i],1]*pk[1];
			gamma[i] = (aixk - b[D[i]])/(-aipk);
			assert(!gamma[i].isNaN);
			alpha = fmin(alpha, gamma[i]);
			assert(!alpha.isNaN);
		}

		uint t = 0;
		for(uint i = 0; i < Dlen; i++)
		{
			if(gamma[i] == alpha)
			{
				t = D[i];
				break;
			}
		}

		if(q == Wk[0])
		{
			Wk[0] = t;
		}
		else if(q == Wk[1])
		{
			Wk[1] = t;
		}
		else
		{
			printf("ruh, rho, idk wtf mate\n");
		}

		Ak[0,0] = A[Wk[0],0];
		Ak[0,1] = A[Wk[0],1];
		Ak[1,0] = A[Wk[1],0];
		Ak[1,1] = A[Wk[1],1];

		xk += alpha*pk;

		assert(!pk[0].isNaN);

		k++;

		if(k > 15)
		{
			printf("k > 15, bailing out\n");
			break;
		}
	}
}

@nogc void ufvmSetup(ref UMesh2 mesh, ref Config config, double[] lastRho, double[] lastU, double[] lastV, double[] lastE, ref double t, ref double dt, string saveFile)
{
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
		printf("mu = %f\n", config.physicalConfig.mu);
	}

	printf("M = %f\tU = %f\n", M, U);

	if(saveFile == "")
	{
		// Setup initial conditions
		foreach(i; mesh.interiorCells)
		{
			mesh.q[i][0] = rho;
			mesh.q[i][1] = rho*u;
			mesh.q[i][2] = rho*v;
			mesh.q[i][3] = p/(gamma - 1.0) + 0.5*rho*(u^^2.0 + v^^2.0);
			lastRho[i] = mesh.q[i][0];
			lastU[i] = mesh.q[i][1];
			lastV[i] = mesh.q[i][2];
			lastE[i] = mesh.q[i][3];
		}
	}
	else
	{
		enforce(loadSolution(mesh, t, dt, saveFile), "Failed to load solution");
		
		foreach(i; mesh.interiorCells)
		{
			lastRho[i] = mesh.q[i][0];
			lastU[i] = mesh.q[i][1];
			lastV[i] = mesh.q[i][2];
			lastE[i] = mesh.q[i][3];
		}
	}

	// Setup bc's
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

			mesh.edges[mesh.bGroups[i][j]].boundaryType = config.boundaries[bcIdx].type;

			mesh.edges[mesh.bGroups[i][j]].bIdx = bcIdx;

			if(mesh.edges[mesh.bGroups[i][j]].boundaryType == BoundaryType.FullState)
			{
				M = config.boundaries[bcIdx].boundaryData[0];
				aoa = config.boundaries[bcIdx].boundaryData[1] * (PI/180);
				rho = config.boundaries[bcIdx].boundaryData[3];

				U = 1.0;
				a = U/M;
				p = (rho*a^^2.0)/gamma;
				u = U*cos(aoa);
				v = U*sin(aoa);

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
			else if(mesh.edges[mesh.bGroups[i][j]].boundaryType == BoundaryType.TempPresInflow)
			{
				enforce(config.boundaries[bcIdx].boundaryData.length == 3, "Constant temperature/pressure/aoa boundary data should have three elements");
				mesh.edges[mesh.bGroups[i][j]].bData = config.boundaries[bcIdx].boundaryData;
			}
		}
	}
}

Datatype vec4dataType;

@nogc updateBoundaryEdge(size_t dims)(ref UMesh2 mesh, Vector!dims[] q, uint i, Config config)
{
	if(config.order == 1)
	{
		mesh.edges[i].q[0] = q[mesh.edges[i].cellIdx[0]];
	}
	else
	{
		auto qM = q[mesh.edges[i].cellIdx[0]];
		auto grad = mesh.cells[mesh.edges[i].cellIdx[0]].gradient;
		auto centroid = mesh.cells[mesh.edges[i].cellIdx[0]].centroid;
		auto mid = mesh.edges[i].mid;

		auto dx = mid[0] - centroid[0];
		auto dy = mid[1] - centroid[1];
		
		for(uint j = 0; j < dims; j++)
		{
			mesh.edges[i].q[0][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[0]].lim[j][0]*grad[j,0]*dx + 
											mesh.cells[mesh.edges[i].cellIdx[0]].lim[j][1]*grad[j,1]*dy;
		}

		if(getPressure(mesh.edges[i].q[0]) < 0)
		{
			double[2] lim = [1.0, 1.0];
			for(uint j = 0; j < dims; j++)
			{
				lim[0] = fmin(lim[0], mesh.cells[mesh.edges[i].cellIdx[0]].lim[j][0]);
				lim[1] = fmin(lim[0], mesh.cells[mesh.edges[i].cellIdx[1]].lim[j][1]);
			}
			for(uint j = 0; j < dims; j++)
			{
				mesh.edges[i].q[0][j] = qM[j] + lim[0]*grad[j,0]*dx + lim[0]*grad[j,1]*dy;
			}
		}
	}
}

// Unstructured finite volume solver
@nogc void ufvmSolver(alias F, size_t dims)(ref Vector!dims[] R, ref Vector!dims[] q, ref UMesh2 mesh, Config config, ref double newDt, ref double Rmax, bool limit, bool dtUpdate)
{
	@nogc void buildGradients(uint[] cellList)
	{
		// Build gradients
		if(config.order > 1)
		{
			foreach(i; cellList)
			{
				Vector!MAX_EDGES[dims] du;
				for(uint j = 0; j < dims; j++)
				{
					du[j] = Vector!MAX_EDGES(0); 
				}
				mesh.cells[i].minQ = q[i];
				mesh.cells[i].maxQ = q[i];
				mesh.cells[i].lim[] = Vector!2(1.0);

				for(uint j = 0; j < mesh.cells[i].nNeighborCells; j++)
				{
					uint idx = mesh.cells[i].neighborCells[j];
					
					for(uint k = 0; k < dims; k++)
					{
						mesh.cells[i].minQ[k] = fmin(mesh.cells[i].minQ[k], q[idx][k]);
						mesh.cells[i].maxQ[k] = fmax(mesh.cells[i].maxQ[k], q[idx][k]);

						du[k][j] = q[idx][k] - q[i][k];
					}
				}

				for(uint j = 0; j < dims; j++)
				{
					auto tmpGrad = mesh.cells[i].gradMat*du[j];
					mesh.cells[i].gradient[j,0] = tmpGrad[0];
					mesh.cells[i].gradient[j,1] = tmpGrad[1];
				}

				if(config.limited && limit)
				{
					for(uint j = 0; j < mesh.cells[i].nNeighborCells; j++)
					{
						auto qM = q[i];
						auto grad = mesh.cells[i].gradient;
						auto centroid = mesh.cells[i].centroid;
						auto mid = mesh.cells[mesh.cells[i].neighborCells[j]].centroid;

						auto dx = mid[0] - centroid[0];
						auto dy = mid[1] - centroid[1];
						auto minQ = mesh.cells[i].minQ;
						auto maxQ = mesh.cells[i].maxQ;
						
						auto qE = Vector!dims(0);
						for(uint k = 0; k < dims; k++)
						{
							double s = 1.0;
							qE[k] = qM[k] + (grad[k,0]*dx + grad[k,1]*dy);

							if(qE[k] < minQ[k])
							{
								s = (minQ[k] - qM[k])/(grad[k,0]*dx + grad[k,1]*dy);
							}
							else if(qE[k] > maxQ[k])
							{
								s = (maxQ[k] - qM[k])/(grad[k,0]*dx + grad[k,1]*dy);
							}

							if(s < 0)
							{
								printf("computed negative limiter\n");
							}

							if(s > 1.0)
							{
								printf("computed limiter greater than 1.0\n");
							}

							mesh.cells[i].lim[k][0] = fmin(mesh.cells[i].lim[k][0], s);
							mesh.cells[i].lim[k][1] = mesh.cells[i].lim[k][0];
						}
					}

					if(config.lpThresh > 0)
					{
						for(uint k = 0; k < dims; k++)
						{
							if(mesh.cells[i].lim[k][0] < config.lpThresh)
							{
								if(mesh.cells[i].nNeighborCells == 3)
								{
									auto A = Matrix!(10,2)(0);
									auto b = Vector!10(0);
									auto c = Vector!2(-abs(mesh.cells[i].gradient[k,0]), -abs(mesh.cells[i].gradient[k,1]));
									auto xk = Vector!2(0);
									A[0,0] = 1;
									A[0,1] = 0;
									A[1,0] = 0;
									A[1,1] = 1;
									A[2,0] = -1;
									A[2,1] = 0;
									A[3,0] = 0;
									A[3,1] = -1;
									b[0] = 0;
									b[1] = 0;
									b[2] = -1;
									b[3] = -1;

									for(uint j = 0, cIdx = 0; j < mesh.cells[i].nNeighborCells; j++, cIdx += 2)
									{
										immutable uint cellIdx = mesh.cells[i].neighborCells[j];

										A[cIdx+4,0] = (mesh.cells[cellIdx].centroid[0] - mesh.cells[i].centroid[0])*mesh.cells[i].gradient[k,0];
										A[cIdx+4,1] = (mesh.cells[cellIdx].centroid[1] - mesh.cells[i].centroid[1])*mesh.cells[i].gradient[k,1];
										
										b[cIdx+4] = mesh.cells[i].minQ[k] - q[i][k];

										A[cIdx+5,0] = -(mesh.cells[cellIdx].centroid[0] - mesh.cells[i].centroid[0])*mesh.cells[i].gradient[k,0];
										A[cIdx+5,1] = -(mesh.cells[cellIdx].centroid[1] - mesh.cells[i].centroid[1])*mesh.cells[i].gradient[k,1];
										
										b[cIdx+5] = q[i][k] - mesh.cells[i].maxQ[k];
									}

									LP!10(A, b, c, xk);

									mesh.cells[i].lim[k][0] = xk[0];
									mesh.cells[i].lim[k][1] = xk[1];

									assert((xk[0] >= -10e-16) && (xk[0] <= (1.0 + 1e-12)));
									assert((xk[1] >= -10e-16) && (xk[1] <= (1.0 + 1e-12)));
								}
								else if(mesh.cells[i].nNeighborCells == 4)
								{
									immutable matsize = 4 + 8;
									auto A = Matrix!(matsize,2)(0);
									auto b = Vector!matsize(0);
									auto c = Vector!2(-abs(mesh.cells[i].gradient[k,0]), -abs(mesh.cells[i].gradient[k,1]));
									auto xk = Vector!2(0);
									A[0,0] = 1;
									A[0,1] = 0;
									A[1,0] = 0;
									A[1,1] = 1;
									A[2,0] = -1;
									A[2,1] = 0;
									A[3,0] = 0;
									A[3,1] = -1;
									b[0] = 0;
									b[1] = 0;
									b[2] = -1;
									b[3] = -1;

									for(uint j = 0, cIdx = 0; j < mesh.cells[i].nNeighborCells; j++, cIdx += 2)
									{
										immutable uint cellIdx = mesh.cells[i].neighborCells[j];

										A[cIdx+4,0] = (mesh.cells[cellIdx].centroid[0] - mesh.cells[i].centroid[0])*mesh.cells[i].gradient[k,0];
										A[cIdx+4,1] = (mesh.cells[cellIdx].centroid[1] - mesh.cells[i].centroid[1])*mesh.cells[i].gradient[k,1];
										
										b[cIdx+4] = mesh.cells[i].minQ[k] - q[i][k];

										A[cIdx+5,0] = -(mesh.cells[cellIdx].centroid[0] - mesh.cells[i].centroid[0])*mesh.cells[i].gradient[k,0];
										A[cIdx+5,1] = -(mesh.cells[cellIdx].centroid[1] - mesh.cells[i].centroid[1])*mesh.cells[i].gradient[k,1];
										
										b[cIdx+5] = q[i][k] - mesh.cells[i].maxQ[k];
									}

									LP!matsize(A, b, c, xk);

									mesh.cells[i].lim[k][0] = xk[0];
									mesh.cells[i].lim[k][1] = xk[1];

									assert((xk[0] >= -10e-16) && (xk[0] <= (1.0 + 1e-12)));
									assert((xk[1] >= -10e-16) && (xk[1] <= (1.0 + 1e-12)));
								}
							}
						}
					}
				
					for(uint j = 0; j < dims; j++)
					{
						mesh.cells[i].gradErr[j] = -abs(mesh.cells[i].gradient[j,0])*mesh.cells[i].lim[j][0] +
												-abs(mesh.cells[i].gradient[j,1])*mesh.cells[i].lim[j][1] + 
													abs(mesh.cells[i].gradient[j,0]) + abs(mesh.cells[i].gradient[j,1]);
					}
				}
			}
		}
	}

	foreach(commIdx, commCells; mesh.commCellSendIdx)
	{
		foreach(i, cell; commCells)
		{
			mesh.sendStateBuffers[commIdx][i] = q[cell];
		}
	}
	mesh.sendRequests.startall;
	mesh.recvRequests.startall;

	// update ghost cell states before we can compute gradients
	foreach(i; mesh.ghostCells)
	{
		auto edge = mesh.edges[mesh.cells[i].edges[0]];
		switch(edge.boundaryType)
			with(BoundaryType)
		{
			case Dirichlet:
				auto centroid = mesh.cells[i].centroid;
				q[i] = config.boundaries[edge.bIdx].dFunc(centroid[0], centroid[1]);
				break;
			case FullState:
				q[edge.cellIdx[1]] = q[edge.cellIdx[0]];
				break;
			case InviscidWall:
				auto cellIdx = edge.cellIdx[0];
				auto v = Vector!2(q[cellIdx][1]/q[cellIdx][0], q[cellIdx][2]/q[cellIdx][0]);
				immutable double x1 = mesh.nodes[edge.nodeIdx[0]][0];
				immutable double y1 = mesh.nodes[edge.nodeIdx[0]][1];
				immutable double x2 = mesh.nodes[edge.nodeIdx[1]][0];
				immutable double y2 = mesh.nodes[edge.nodeIdx[1]][1];
				auto newV = Vector!2(0.0);
				if(x1 != x2)
				{
					immutable double m = (y2 - y1)/(x2 - x1);
					auto reflection = Matrix!(2,2)(1 - m^^2.0, 2.0*m, 2.0*m, m^^2.0 - 1)*1.0/(1 + m^^2.0);
					newV = reflection*v;
				}
				else
				{
					newV[0] = -v[0];
					newV[1] = v[1];
				}

				auto cellIdx2 = edge.cellIdx[1];
				q[cellIdx2][0] = q[cellIdx][0];
				q[cellIdx2][1] = q[cellIdx][0]*newV[0];
				q[cellIdx2][2] = q[cellIdx][0]*newV[1];
				q[cellIdx2][3] = q[cellIdx][3];

			 	break;
			case ViscousWall:
				auto cellIdx = edge.cellIdx[0];

				auto cellIdx2 = edge.cellIdx[1];
				q[cellIdx2][0] = q[cellIdx][0];
				q[cellIdx2][1] = -q[cellIdx][1];
				q[cellIdx2][2] = -q[cellIdx][2];
				q[cellIdx2][3] = q[cellIdx][3];
				break;
			case ConstPressure:
				auto cellIdx = edge.cellIdx[0];
				auto cellIdx2 = edge.cellIdx[1];

				double rhop = q[cellIdx][0];
				double up = q[cellIdx][1]/q[cellIdx][0];
				double vp = q[cellIdx][2]/q[cellIdx][0];
				auto Vp = Vector!2(up, vp);
				double pp = getPressure(q[cellIdx]); // lol
				double ap = sqrt(gamma*(pp/rhop));
				double Jp = up + (2.0*ap)/(gamma - 1);
				double Sp = pp/(rhop^^gamma);
				double pb = edge.bData[0];
				double rhob = (pb/Sp)^^(1.0/gamma);
				double ab = sqrt(gamma*(pb/rhob));
				double ub = Jp - (2.0*ab)/(gamma - 1);
				double vb = vp;
				auto Vb = Vector!2(ub, vb);
				double rEb = pb/(gamma - 1) + 0.5*rhob*Vb.magnitude^^2.0;
				auto Vbg = Vp - Vp.dot(edge.normal)*edge.normal + Vb.dot(edge.normal)*edge.normal;

				q[cellIdx2][0] = rhob;
				q[cellIdx2][1] = rhob*Vbg[0];
				q[cellIdx2][2] = rhob*Vbg[1];
				q[cellIdx2][3] = rEb;

				break;

			case Symmetry:
				auto cellIdx = edge.cellIdx[0];
				auto v = Vector!2(q[cellIdx][1]/q[cellIdx][0], q[cellIdx][2]/q[cellIdx][0]);
				immutable double x1 = mesh.nodes[edge.nodeIdx[0]][0];
				immutable double y1 = mesh.nodes[edge.nodeIdx[0]][1];
				immutable double x2 = mesh.nodes[edge.nodeIdx[1]][0];
				immutable double y2 = mesh.nodes[edge.nodeIdx[1]][1];
				auto newV = Vector!2(0.0);
				if(x1 != x2)
				{
					immutable double m = (y2 - y1)/(x2 - x1);
					auto reflection = Matrix!(2,2)(1 - m^^2.0, 2.0*m, 2.0*m, m^^2.0 - 1)*1.0/(1 + m^^2.0);
					newV = reflection*v;
				}
				else
				{
					newV[0] = -v[0];
					newV[1] = v[1];
				}

				auto cellIdx2 = edge.cellIdx[1];
				q[cellIdx2][0] = q[cellIdx][0];
				q[cellIdx2][1] = q[cellIdx][0]*newV[0];
				q[cellIdx2][2] = q[cellIdx][0]*newV[1];
				q[cellIdx2][3] = q[cellIdx][3];
				break;
			
			case TempPresInflow:
				auto cellIdx = edge.cellIdx[0];
				auto cellIdx2 = edge.cellIdx[1];
				
				double Tt = edge.bData[0];
				double Pt = edge.bData[1];
				double aoa = edge.bData[2] * PI/180;

				double rhop = q[cellIdx][0];
				double up = q[cellIdx][1]/q[cellIdx][0];
				double vp = q[cellIdx][2]/q[cellIdx][0];
				auto Vp = Vector!2(up, vp);
				double pp = getPressure(q[cellIdx]); // lol
				double ap = sqrt(gamma*(pp/rhop));
				double Jp = up + (2.0*ap)/(gamma - 1);

				auto n_in = Vector!2(cos(aoa), sin(aoa));
				auto d = n_in.dot(edge.normal);
				double a = gamma*config.physicalConfig.R*Tt*d^^2.0 - 0.5*(gamma - 1)*Jp^^2.0;
				double b = (4.0*gamma*config.physicalConfig.R*Tt*d)/(gamma - 1);
				double c = (4.0*gamma*config.physicalConfig.R*Tt)/(gamma - 1)^^2.0 - Jp^^2.0;

				auto M1 = (-b + sqrt(b^^2.0 - 4*a*c))/(2.0*a);
				auto M2 = (-b - sqrt(b^^2.0 - 4*a*c))/(2.0*a);
				double Mb = 0;
				if((M1 >= 0) && (M2 >= 0))
				{
					Mb = min(M1, M2);
				}
				else if((M1 >= 0) && (M2 < 0))
				{
					Mb = M1;
				}
				else if((M1 < 0) && (M2 >= 0))
				{
					Mb = M2;
				}
				else
				{
					printf("M1 = %f, M2 = %f\n", M1, M2);
					printf("a = %f, b = %f, c = %f\n", a, b, c);
					printf("Jp = %f\n", Jp);
					printf("%f\n", (4.0*gamma*config.physicalConfig.R*Tt)/(gamma - 1)^^2.0);
					enforce(false, "Both computed mach numbers are negative");
				}

				double Tb = Tt/(1 + 0.5*(gamma - 1)*Mb^^2.0);
				double pb = Pt*(Tb/Tt)^^(gamma/(gamma - 1));
				double rhob = pb/(config.physicalConfig.R*Tb);
				double ab = sqrt(gamma*(pb/rhob));
				auto vb = Mb*ab*n_in;
				double rEb = pb/(gamma - 1) + 0.5*rhob*vb.magnitude^^2.0;

				q[cellIdx2][0] = rhob;
				q[cellIdx2][0] = rhob*vb[0];
				q[cellIdx2][0] = rhob*vb[1];
				q[cellIdx2][0] = rEb;

				break;

			default:
				enforce(false, "Unsupported boundary type");
				break;
		}
	}

	buildGradients(mesh.nonCommCells);

	mesh.recvRequests.waitall(mesh.statuses);
	foreach(commIdx, commCells; mesh.commCellRecvIdx)
	{
		foreach(i, cell; commCells)
		{
			q[cell] = mesh.recvStateBuffers[commIdx][i];
		}
	}

	buildGradients(mesh.needCommCells);

	foreach(commIdx, commCells; mesh.commCellSendIdx)
	{
		foreach(i, cell; commCells)
		{
			mesh.sendGradBuffers[commIdx][i] = mesh.cells[cell].gradient;
		}
	}
	
	mesh.sendGradRequests.startall;
	mesh.recvGradRequests.startall;

	// need to do a half update of comm edges
	foreach(commIdx, commEdges; mesh.commEdgeIdx)
	{
		foreach(i; commEdges)
		{
			if(config.order == 1)
			{
				mesh.edges[i].q[0] = q[mesh.edges[i].cellIdx[0]];
			}
			else
			{
				auto qM = q[mesh.edges[i].cellIdx[0]];
				auto grad = mesh.cells[mesh.edges[i].cellIdx[0]].gradient;
				auto centroid = mesh.cells[mesh.edges[i].cellIdx[0]].centroid;
				auto mid = mesh.edges[i].mid;
				
				auto dx = mid[0] - centroid[0];
				auto dy = mid[1] - centroid[1];
				
				for(uint j = 0; j < dims; j++)
				{
					mesh.edges[i].q[0][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[0]].lim[j][0]*grad[j,0]*dx + 
													mesh.cells[mesh.edges[i].cellIdx[0]].lim[j][1]*grad[j,1]*dy;
				}

				if(getPressure(mesh.edges[i].q[0]) < 0)
				{
					double[2] lim = [1.0, 1.0];
					for(uint j = 0; j < dims; j++)
					{
						lim[0] = fmin(lim[0], mesh.cells[mesh.edges[i].cellIdx[0]].lim[j][0]);
						lim[1] = fmin(lim[1], mesh.cells[mesh.edges[i].cellIdx[0]].lim[j][1]);
					}

					for(uint j = 0; j < dims; j++)
					{
						mesh.edges[i].q[0][j] = qM[j] + lim[0]*grad[j,0]*dx + lim[1]*grad[j,1]*dy;
					}
				}
			}
		}
	}

	foreach(commIdx, commEdges; mesh.commEdgeIdx)
	{
		foreach(i, edge; commEdges)
		{
			mesh.sendStateBuffers[commIdx][i] = mesh.edges[edge].q[0];
		}
	}

	mesh.sendRequests.startall;
	mesh.recvRequests.startall;

	// update edge values and compute fluxes on boundary edges.
	foreach(i; mesh.boundaryEdges)
	{
		mesh.updateBoundaryEdge!dims(q, i, config);

		switch(mesh.edges[i].boundaryType)
			with(BoundaryType)
		{
			case Dirichlet:
				double x = mesh.cells[mesh.edges[i].cellIdx[1]].centroid[0];
				double y = mesh.cells[mesh.edges[i].cellIdx[1]].centroid[1];

				mesh.cells[mesh.edges[i].cellIdx[1]].gradient = config.boundaries[mesh.edges[i].bIdx].dFuncDerivative(x, y);
				auto mid = mesh.edges[i].mid;
				mesh.edges[i].q[1] = config.boundaries[mesh.edges[i].bIdx].dFunc(mid[0], mid[1]);
				goto case FullState;

			case FullState:
				auto qL = mesh.edges[i].q[0];
				auto qR = mesh.edges[i].q[1];

				if(config.viscosity)
				{
					auto qAve = 0.5*(qL + qR);
					auto grad1 = mesh.cells[mesh.edges[i].cellIdx[0]].gradient;
					auto grad2 = mesh.cells[mesh.edges[i].cellIdx[1]].gradient;
					auto dqAve = 0.5*(grad1 + grad2);

					auto Fv = diffusiveFlux!dims(config.physicalConfig.Pr, config.physicalConfig.mu, qAve, dqAve, mesh.edges[i].normal);
					mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax) - Fv;
				}
				else
				{
					mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax);
				}

				immutable bool haveNan = (mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN);
				enforce!EdgeException(!haveNan, "Got NaN on fullstate edge", mesh.edges[i]);
				break;

			case InviscidWall:
				Vector!2 velP = (1/mesh.edges[i].q[0][0])*Vector!2(mesh.edges[i].q[0][1], mesh.edges[i].q[0][2]);
				auto vel = (velP - (velP.dot(mesh.edges[i].normal))*mesh.edges[i].normal).magnitude;
				double p = (gamma - 1)*(mesh.edges[i].q[0][3] - 0.5*mesh.edges[i].q[0][0]*vel^^2);
				double a = sqrt(gamma*(p/mesh.edges[i].q[0][0]));
				if(p < 0)
				{
					p = 1.0e-12;
					//printf("pressure less than 0 at wall\n");
				}
				mesh.edges[i].flux = Vector!dims(0, p*mesh.edges[i].normal[0], p*mesh.edges[i].normal[1], 0);
				mesh.edges[i].sMax = std.math.abs(a);
				
				auto qL = mesh.edges[i].q[0];
				auto qR = q[mesh.edges[i].cellIdx[1]];
				mesh.edges[i].q[1] = qR;
				immutable bool haveNan = (mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN);
				enforce!EdgeException(!haveNan, "Got NaN on inviscid wall edge", mesh.edges[i]);

				break;

			case ViscousWall:
				auto vel = 0.0;
				double p = (gamma - 1)*(mesh.edges[i].q[0][3] - 0.5*mesh.edges[i].q[0][0]*vel^^2);
				double a = sqrt(gamma*(p/mesh.edges[i].q[0][0]));
				if(p < 0)
				{
					p = 1.0e-12;
					//printf("pressure less than 0 at wall\n");
				}

				auto qL = mesh.edges[i].q[0];
				auto qR = qL;//q[mesh.edges[i].cellIdx[1]];
				qR[1] *= -1;
				qR[2] *= -1;
				
				mesh.edges[i].q[1] = qR;

				auto qAve = 0.5*(qL + qR);

				qAve[1] = 0.0;
				qAve[2] = 0.0;


				auto dqAve = mesh.cells[mesh.edges[i].cellIdx[0]].gradient;
				// modify energy term to be an adiabatic wall.
				// TODO: Generalize this. Equations in EbbCFD book 1 page 5
				dqAve[3,0] = (qAve[3]/qAve[0])*dqAve[0,0];
				dqAve[3,1] = (qAve[3]/qAve[0])*dqAve[0,1];
				auto Fv = diffusiveFlux!dims(config.physicalConfig.Pr, config.physicalConfig.mu, qAve, dqAve, mesh.edges[i].normal);

				mesh.edges[i].flux = Vector!dims(0, p*mesh.edges[i].normal[0], p*mesh.edges[i].normal[1], 0) - Fv;
				mesh.edges[i].sMax = std.math.abs(a);
				
				immutable bool haveNan = (mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN);
				enforce!EdgeException(!haveNan, "Got NaN on viscous wall edge", mesh.edges[i]);
				break;
			case ConstPressure:
				auto cellIdx2 = mesh.edges[i].cellIdx[1];

				auto p = getPressure(q[cellIdx2]);
				auto u = q[cellIdx2][1]/q[cellIdx2][0];
				auto v = q[cellIdx2][2]/q[cellIdx2][0];
				double a = sqrt(gamma*(p/q[cellIdx2][0]));
				mesh.edges[i].q[1] = q[cellIdx2];
				mesh.edges[i].flux = convectiveFlux!4(p, u, v, q[cellIdx2][0], q[cellIdx2][3], mesh.edges[i].normal);
				//mesh.edges[i].flux = Vector!dims(0, p*mesh.edges[i].normal[0], p*mesh.edges[i].normal[1], 0) - Fv;
				mesh.edges[i].sMax = std.math.abs(a);
				immutable bool haveNan = (mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN);
				enforce!EdgeException(!haveNan, "Got NaN on viscous wall edge", mesh.edges[i]);

				break;

			case Symmetry:
				auto edge = mesh.edges[i];
				auto v = Vector!2(mesh.edges[i].q[0][1]/mesh.edges[i].q[0][0], mesh.edges[i].q[0][2]/mesh.edges[i].q[0][0]);
				immutable double x1 = mesh.nodes[edge.nodeIdx[0]][0];
				immutable double y1 = mesh.nodes[edge.nodeIdx[0]][1];
				immutable double x2 = mesh.nodes[edge.nodeIdx[1]][0];
				immutable double y2 = mesh.nodes[edge.nodeIdx[1]][1];
				auto newV = Vector!2(0.0);
				if(x1 != x2)
				{
					immutable double m = (y2 - y1)/(x2 - x1);
					auto reflection = Matrix!(2,2)(1 - m^^2.0, 2.0*m, 2.0*m, m^^2.0 - 1)*1.0/(1 + m^^2.0);
					newV = reflection*v;
				}
				else
				{
					newV[0] = -v[0];
					newV[1] = v[1];
				}

				mesh.edges[i].q[1][0] = mesh.edges[i].q[0][0];
				mesh.edges[i].q[1][1] = mesh.edges[i].q[0][0]*newV[0];
				mesh.edges[i].q[1][2] = mesh.edges[i].q[0][0]*newV[1];
				mesh.edges[i].q[1][3] = mesh.edges[i].q[0][3];

				auto qL = mesh.edges[i].q[0];
				auto qR = mesh.edges[i].q[1];

				if(config.viscosity)
				{
					auto qAve = 0.5*(qL + qR);
					auto cellIdx = mesh.edges[i].cellIdx[0];
					auto gradp = mesh.cells[cellIdx].gradient;
					auto nOuter = edge.normal*edge.normal.transpose;
					auto A = Matrix!(dims,dims).Identity;
					auto V = Matrix!(dims-2,dims-2).Identity - nOuter;
					A[1,1] = V[0,0];
					A[1,2] = V[0,1];
					A[2,1] = V[1,0];
					A[2,2] = V[1,1];
					auto tmp1 = A*gradp;
					auto tmp2 = Matrix!(dims-2,dims-2).Identity - 2*nOuter;
					auto tmp3 = tmp1*tmp2;
					auto bGrad = gradp*nOuter + tmp3;

					auto Fv = diffusiveFlux!dims(config.physicalConfig.Pr, config.physicalConfig.mu, qAve, bGrad, mesh.edges[i].normal);
					mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax) - Fv;
				}
				else
				{
					mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax);
				}

				immutable bool haveNan = (mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN);
				enforce!EdgeException(!haveNan, "Got NaN on interior edge", mesh.edges[i]);

				break;
			case TempPresInflow:

				mesh.edges[i].q[1] = q[mesh.edges[i].cellIdx[1]];
				auto qL = mesh.edges[i].q[0];
				auto qR = mesh.edges[i].q[1];

				if(config.viscosity)
				{
					auto qAve = 0.5*(qL + qR);
					auto grad1 = mesh.cells[mesh.edges[i].cellIdx[0]].gradient;
					//auto grad2 = mesh.cells[mesh.edges[i].cellIdx[1]].gradient;
					//auto dqAve = 0.5*(grad1 + grad2);

					auto Fv = diffusiveFlux!dims(config.physicalConfig.Pr, config.physicalConfig.mu, qAve, grad1, mesh.edges[i].normal);
					mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax) - Fv;
				}
				else
				{
					mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax);
				}

				immutable bool haveNan = (mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN);
				enforce!EdgeException(!haveNan, "Got NaN on stagnation temp/press inflow edge", mesh.edges[i]);

				break;

			default:
				enforce(false, "Unsupported boundary type");
		}
	}

	// update edge values and compute fluxes on interior edges.
	foreach(i; mesh.interiorEdges)
	{
		if(config.order == 1)
		{
			mesh.edges[i].q[0] = q[mesh.edges[i].cellIdx[0]];
			mesh.edges[i].q[1] = q[mesh.edges[i].cellIdx[1]];
		}
		else
		{
			for(uint k = 0; k < 2; k++)
			{
				auto qM = q[mesh.edges[i].cellIdx[k]];
				auto grad = mesh.cells[mesh.edges[i].cellIdx[k]].gradient;
				auto centroid = mesh.cells[mesh.edges[i].cellIdx[k]].centroid;
				auto mid = mesh.edges[i].mid;

				auto dx = mid[0] - centroid[0];
				auto dy = mid[1] - centroid[1];

				for(uint j = 0; j < dims; j++)
				{
					mesh.edges[i].q[k][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[k]].lim[j][0]*grad[j,0]*dx + 
													mesh.cells[mesh.edges[i].cellIdx[k]].lim[j][1]*grad[j,1]*dy;
				}

				if(getPressure(mesh.edges[i].q[k]) < 0)
				{
					double[2] lim = [1.0, 1.0];
					for(uint j = 0; j < dims; j++)
					{
						lim[0] = fmin(lim[0], mesh.cells[mesh.edges[i].cellIdx[k]].lim[j][0]);
						lim[1] = fmin(lim[1], mesh.cells[mesh.edges[i].cellIdx[k]].lim[j][1]);
					}

					for(uint j = 0; j < dims; j++)
					{
						mesh.edges[i].q[k][j] = qM[j] + lim[0]*grad[j,0]*dx + lim[1]*grad[j,1]*dy;
					}
				}
			}
		}

		auto qL = mesh.edges[i].q[0];
		auto qR = mesh.edges[i].q[1];

		if(config.viscosity)
		{
			auto qAve = 0.5*(qL + qR);
			auto grad1 = mesh.cells[mesh.edges[i].cellIdx[0]].gradient;
			auto grad2 = mesh.cells[mesh.edges[i].cellIdx[1]].gradient;
			auto dqAve = 0.5*(grad1 + grad2);

			auto Fv = diffusiveFlux!dims(config.physicalConfig.Pr, config.physicalConfig.mu, qAve, dqAve, mesh.edges[i].normal);
			mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax) - Fv;
		}
		else
		{
			mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax);
		}

		immutable bool haveNan = (mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN);
		enforce!EdgeException(!haveNan, "Got NaN on interior edge", mesh.edges[i]);
	}

	mesh.recvRequests.waitall(mesh.statuses);
	foreach(commIdx, commEdges; mesh.commEdgeIdx)
	{
		foreach(i, edge; commEdges)
		{
			mesh.edges[edge].q[1] = mesh.recvStateBuffers[commIdx][i];
		}
	}

	mesh.recvGradRequests.waitall(mesh.statuses);
	foreach(commIdx, commCells; mesh.commCellRecvIdx)
	{
		foreach(i, cell; commCells)
		{
			mesh.cells[cell].gradient = mesh.recvGradBuffers[commIdx][i];
		}
	}

	foreach(commIdx, commEdges; mesh.commEdgeIdx)
	{
		foreach(i; commEdges)
		{
			auto qL = mesh.edges[i].q[0];
			auto qR = mesh.edges[i].q[1];

			if(config.viscosity)
			{
				auto qAve = 0.5*(qL + qR);
				auto grad1 = mesh.cells[mesh.edges[i].cellIdx[0]].gradient;
				auto grad2 = mesh.cells[mesh.edges[i].cellIdx[1]].gradient;

				auto dqAve = 0.5*(grad1 + grad2);
				auto Fv = diffusiveFlux!dims(config.physicalConfig.Pr, config.physicalConfig.mu, qAve, dqAve, mesh.edges[i].normal);
				mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax) - Fv;
			}
			else
			{
				mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax);
			}

			immutable bool haveNan = (mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN);
			enforce!EdgeException(!haveNan, "Got NaN on comm edge", mesh.edges[i]);
		}
	}

	newDt = double.infinity;

	foreach(i; mesh.interiorCells)
	{
		R[i] = Vector!dims(0);
		double sAve = 0;
		// integrate fluxes over cell edges
		for(uint j = 0; j < mesh.cells[i].nEdges; j++)
		{
			R[i] += mesh.cells[i].fluxMultiplier[j]*mesh.edges[mesh.cells[i].edges[j]].len*mesh.edges[mesh.cells[i].edges[j]].flux;
			sAve += mesh.edges[mesh.cells[i].edges[j]].len*mesh.edges[mesh.cells[i].edges[j]].sMax;
		}
		sAve /= mesh.cells[i].perim;

		newDt = fmin(newDt, (config.CFL*mesh.cells[i].d)/sAve);

		if(config.localTimestep)
		{
			if(dtUpdate)
			{
				mesh.cells[i].dt = config.CFL*mesh.cells[i].d/sAve;
				if(mesh.cells[i].dt.isNaN)
				{
					//printf("dt is nan\n");
					mesh.cells[i].dt = 0;
				}
			}
		}

		for(uint j = 0; j < dims; j++)
		{
			if(std.math.abs(R[i][j]) > Rmax)
			{
				Rmax = std.math.abs(R[i][j]);
			}
		}
		
		R[i] *= -(1.0/mesh.cells[i].area);
	}

	mesh.comm.allreduce(newDt, newDt, MPI_MIN);
	mesh.comm.allreduce(Rmax, Rmax, MPI_MAX);
}
