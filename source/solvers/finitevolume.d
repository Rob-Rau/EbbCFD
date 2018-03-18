/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.solvers.finitevolume;

import core.atomic;
import core.stdc.stdio : fopen, fwrite, fopen, printf, snprintf;
import core.sys.posix.signal;

import std.algorithm;
import std.conv;
import std.exception;
import std.file;
import std.math;
import std.range;
import std.stdio;
import std.typecons;

import numd.utility;
import numd.linearalgebra.matrix;

import ebb.config;
import ebb.exception;
import ebb.flux;
import ebb.mesh;
import ebb.mpid;
import ebb.solve;

struct FiniteVolume(SolverParams params)
{
	immutable size = params.size;
	immutable dims = params.dims;

	Config config;

	Communication!params comm;
	IPhysics!params physics;
	immutable Mesh mesh;

	this(immutable Mesh mesh, IPhysics!params physics, Communication!params comm)
	{
		this.mesh = mesh;
		this.physics = physics;
		this.comm = comm;

		qL = new Vector!size[mesh.edges.length];
		qR = new Vector!size[mesh.edges.length];
		fluxes = new Vector!size[mesh.edges.length];
		waveSpeeds = new double[mesh.edges.length];
		if(config.order > 1)
		{
			limits = new Matrix!(size,dims)[mesh.edges.length];
			minQ = new Vector!size[mesh.edges.length];
			maxQ = new Vector!size[mesh.edges.length];
			dQ = new Matrix!(size,dims)[mesh.cells.length];
			if(physics.needEdgeGradients)
			{
				dQl = new Matrix!(size,dims)[mesh.edges.length];
				dQr = new Matrix!(size,dims)[mesh.edges.length];
			}
		}
	}

	@property size_t integrationSize() { return mesh.interiorCells.length; }

	@nogc void updateTimestep(double[] dts)
	{
		foreach(i; mesh.interiorCells)
		{
			double sAve = 0;
			for(uint j = 0; j < mesh.cells[i].nEdges; j++)
			{
				auto edge = mesh.edges[mesh.cells[i].edges[j]];
				sAve += edge.len*waveSpeeds[j];
			}
			sAve /= mesh.cells[i].perim;
			dts[i] = (config.CFL*mesh.cells[i].d)/sAve;
		}
	}

	@nogc double updateTimestep()
	{
		double dt = double.infinity;
		foreach(i; mesh.interiorCells)
		{
			double sAve = 0;
			for(uint j = 0; j < mesh.cells[i].nEdges; j++)
			{
				auto edge = mesh.edges[mesh.cells[i].edges[j]];
				sAve += edge.len*waveSpeeds[j];
			}
			sAve /= mesh.cells[i].perim;
			dt = fmin(dt, (config.CFL*mesh.cells[i].d)/sAve);
		}
		return dt;
	}

	@nogc private void invert(ref Matrix!(dims,dims) inv, ref Matrix!(dims,dims) mat)
	{
		inv[0] = mat[3];
		inv[1] = -mat[1];
		inv[2] = -mat[2];
		inv[3] = mat[0];
		assert(mat.determinant != 0.0);
		inv *= (1.0/mat.determinant);
	}

	@nogc private void LP(uint m)(ref Matrix!(m, dims) A, ref Vector!m b, ref Vector!dims c, ref Vector!dims xk)
	{
		uint[2] Wk = [0, 1];
		auto AkInv = Matrix!(dims,dims)(A[Wk[0],0], A[Wk[0],1], A[Wk[1],0], A[Wk[1],1]);
		auto Ak = Matrix!(dims,dims)(A[Wk[0],0], A[Wk[0],1], A[Wk[1],0], A[Wk[1],1]);
		auto pk = Vector!dims(0);
		auto lam = Vector!dims(0);
		auto e = Vector!dims(0);
		bool done = false;
		uint k = 0;
		uint[m - 2] D;
		double[m - 2] gamma;
		uint Dlen = 0;

		@nogc void invertTrans(ref Matrix!(dims,dims) inv, ref Matrix!(dims,dims) mat)
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
	
	// For boundary edges, this is the edge value of the actual cell
	private Vector!size[] qL; // qL -> edge.q[0]
	// For boundary edges, this is the edge value of the ghost cell
	private Vector!size[] qR; // qR -> edge.q[1]
	private Vector!size[] fluxes;
	private double[] waveSpeeds;
	private Matrix!(size,dims)[] dQ; // cell centered
	private Matrix!(size,dims)[] dQl; // edge gradients
	private Matrix!(size,dims)[] dQr; // edge gradients
	private Matrix!(size,dims)[] limits;
	private Vector!size[] minQ;
	private Vector!size[] maxQ;
	
	@nogc void buildGradients(bool limit)(Vector!size[] Q, immutable(size_t[]) cellList)
	{
		foreach(i; cellList)
		{
			Vector!MAX_EDGES[dims] du;
			for(uint j = 0; j < dims; j++)
			{
				du[j] = Vector!MAX_EDGES(0); 
			}
			minQ[i] = Q[i];
			maxQ[i] = Q[i];

			limits[i] = Matrix!(size,dims)(1.0);

			for(uint j = 0; j < mesh.cells[i].nNeighborCells; j++)
			{
				auto idx = mesh.cells[i].neighborCells[j];
				
				for(uint k = 0; k < dims; k++)
				{
					minQ[i][k] = fmin(minQ[i][k], Q[idx][k]);
					maxQ[i][k] = fmax(maxQ[i][k], Q[idx][k]);

					du[k][j] = Q[idx][k] - Q[i][k];
				}
			}

			for(uint j = 0; j < dims; j++)
			{
				auto tmpGrad = mesh.cells[i].gradMat*du[j];
				dQ[i][j,0] = tmpGrad[0];
				dQ[i][j,1] = tmpGrad[1];
			}

			if(config.limited && limit)
			{
				for(uint j = 0; j < mesh.cells[i].nNeighborCells; j++)
				{
					auto qM = Q[i];
					auto grad = dQ[i];
					auto centroid = mesh.cells[i].centroid;
					auto mid = mesh.cells[mesh.cells[i].neighborCells[j]].centroid;

					auto dx = mid[0] - centroid[0];
					auto dy = mid[1] - centroid[1];
					auto minQ = minQ[i];
					auto maxQ = maxQ[i];
					
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

						limits[i][k, 0] = fmin(limits[i][k, 0], s);
						limits[i][k, 1] = limits[i][k, 0];
					}
				}

				if(config.lpThresh > 0)
				{
					for(uint k = 0; k < dims; k++)
					{
						if(limits[i][k,0] < config.lpThresh)
						{
							if(mesh.cells[i].nNeighborCells == 3)
							{
								auto A = Matrix!(10,dims)(0);
								auto b = Vector!10(0);
								auto c = Vector!dims(-abs(dQ[i][k,0]), -abs(dQ[i][k,1]));
								auto xk = Vector!dims(0);
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
									auto cellIdx = mesh.cells[i].neighborCells[j];

									A[cIdx+4,0] = (mesh.cells[cellIdx].centroid[0] - mesh.cells[i].centroid[0])*dQ[i][k,0];
									A[cIdx+4,1] = (mesh.cells[cellIdx].centroid[1] - mesh.cells[i].centroid[1])*dQ[i][k,1];
									
									b[cIdx+4] = minQ[i][k] - Q[i][k];

									A[cIdx+5,0] = -(mesh.cells[cellIdx].centroid[0] - mesh.cells[i].centroid[0])*dQ[i][k,0];
									A[cIdx+5,1] = -(mesh.cells[cellIdx].centroid[1] - mesh.cells[i].centroid[1])*dQ[i][k,1];
									
									b[cIdx+5] = Q[i][k] - maxQ[i][k];
								}

								LP!10(A, b, c, xk);

								limits[i][k,0] = xk[0];
								limits[i][k,1] = xk[1];

								assert((xk[0] >= -10e-16) && (xk[0] <= (1.0 + 1e-12)));
								assert((xk[1] >= -10e-16) && (xk[1] <= (1.0 + 1e-12)));
							}
							else if(mesh.cells[i].nNeighborCells == 4)
							{
								immutable matsize = 4 + 8;
								auto A = Matrix!(matsize,dims)(0);
								auto b = Vector!matsize(0);
								auto c = Vector!dims(-abs(dQ[i][k,0]), -abs(dQ[i][k,1]));
								auto xk = Vector!dims(0);
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
									auto cellIdx = mesh.cells[i].neighborCells[j];

									A[cIdx+4,0] = (mesh.cells[cellIdx].centroid[0] - mesh.cells[i].centroid[0])*dQ[i][k,0];
									A[cIdx+4,1] = (mesh.cells[cellIdx].centroid[1] - mesh.cells[i].centroid[1])*dQ[i][k,1];
									
									b[cIdx+4] = minQ[i][k] - Q[i][k];

									A[cIdx+5,0] = -(mesh.cells[cellIdx].centroid[0] - mesh.cells[i].centroid[0])*dQ[i][k,0];
									A[cIdx+5,1] = -(mesh.cells[cellIdx].centroid[1] - mesh.cells[i].centroid[1])*dQ[i][k,1];
									
									b[cIdx+5] = Q[i][k] - maxQ[i][k];
								}

								LP!matsize(A, b, c, xk);

								limits[i][k,0] = xk[0];
								limits[i][k,1] = xk[1];

								assert((xk[0] >= -10e-16) && (xk[0] <= (1.0 + 1e-12)));
								assert((xk[1] >= -10e-16) && (xk[1] <= (1.0 + 1e-12)));
							}
						}
					}
				}
			}
		}
	}

	enum UpdateMode
	{
		Right,
		Left,
		Both
	}

	@nogc void updateEdgeState(UpdateMode mode)(immutable(size_t[]) edgeList, Vector!size[] Q)
	{
		@nogc void computeEdgeState(ref Vector!size q, size_t edgeIdx, size_t cellIdx)
		{
			auto qM = Q[cellIdx];
			auto grad = dQ[cellIdx];
			auto centroid = mesh.cells[cellIdx].centroid;
			auto mid = mesh.edges[edgeIdx].mid;

			auto dx = mid[0] - centroid[0];
			auto dy = mid[1] - centroid[1];

			for(uint j = 0; j < size; j++)
			{
				q[j] = qM[j] + limits[cellIdx][j, 0]*grad[j,0]*dx + 
							   limits[cellIdx][j, 1]*grad[j,1]*dy;
			}

			if(physics.checkPhysics(q) < 0)
			{
				double[2] lim = [1.0, 1.0];
				for(uint j = 0; j < dims; j++)
				{
					lim[0] = fmin(lim[0], limits[cellIdx][j, 0]);
					lim[1] = fmin(lim[1], limits[cellIdx][j, 1]);
				}

				for(uint j = 0; j < dims; j++)
				{
					q[j] = qM[j] + lim[0]*grad[j,0]*dx + lim[1]*grad[j,1]*dy;
				}
			}
		}

		if(config.order == 1)
		{
			foreach(i; edgeList)
			{
				static if((mode == UpdateMode.Right) || (mode == UpdateMode.Both))
				{
					qR[i] = Q[mesh.edges[i].cellIdxR];
				}

				static if((mode == UpdateMode.Left) || (mode == UpdateMode.Both))
				{
					qL[i] = Q[mesh.edges[i].cellIdxL];
				}
			}
		}
		else
		{
			foreach(i; edgeList)
			{
				static if((mode == UpdateMode.Right) || (mode == UpdateMode.Both))
				{
					computeEdgeState(qR[i], i, mesh.edges[i].cellIdxR);
				}

				static if((mode == UpdateMode.Left) || (mode == UpdateMode.Both))
				{
					computeEdgeState(qL[i], i, mesh.edges[i].cellIdxL);
				}
			}

			if(physics.needEdgeGradients)
			{
				foreach(i; edgeList)
				{
					static if((mode == UpdateMode.Right) || (mode == UpdateMode.Both))
					{
						dQr[i] = dQ[mesh.edges[i].cellIdxR];
					}

					static if((mode == UpdateMode.Left) || (mode == UpdateMode.Both))
					{
						dQl[i] = dQ[mesh.edges[i].cellIdxL];
					}
				}
			}
		}
	}

	@nogc void solve(bool limit)(ref Vector!size[] dQdt, Vector!size[] Q)
	{
		comm.startSendRecvState(Q, mesh.commCellSendIdx);

		/// update ghost cells
		physics.updateGhostCells(Q, mesh);

		/// build gradients regular cells
		if(config.order > 1)
		{
			buildGradients!limit(Q, mesh.nonCommCells);
		}

		comm.finishRecvState(Q, mesh.commCellRecvIdx);

		/// build gradients comm cells
		if(config.order > 1)
		{
			buildGradients!limit(Q, mesh.needCommCells);
		}

		if(physics.needEdgeGradients)
		{
			comm.startSendRecvStateGradient(dQ, mesh.commCellSendIdx);
		}

		/// update qR for comm cells
		foreach(commIdx, commEdges; mesh.commEdgeIdx)
		{
			updateEdgeState!(UpdateMode.Right)(commEdges, Q);
		}

		comm.startSendRecvState(qR, mesh.commEdgeIdx);

		updateEdgeState!(UpdateMode.Left)(mesh.boundaryEdges, Q);
		/// Update qR for boundaries
		physics.updateBoundaryEdgeFluxes(fluxes, waveSpeeds, qR, qL.assumeUnique, dQl.assumeUnique, dQr.assumeUnique, Q.assumeUnique, mesh);

		/// updated qL and qR
		updateEdgeState!(UpdateMode.Both)(mesh.interiorEdges, Q);

		comm.finishRecvState(qL, mesh.commEdgeIdx);
		if(physics.needEdgeGradients)
		{
			comm.finishRecvStateGradient(dQ, mesh.commCellRecvIdx);
		}

		// Compute the rest of the edge fluxes on all non-boundary edges.
		auto sliceEnd = mesh.interiorEdges.length + mesh.commEdgeIdx.fold!((res, x) => res += x.length)(0);
		physics.computeFluxes(fluxes[0..sliceEnd], waveSpeeds[0..sliceEnd], qL[0..sliceEnd], qR[0..sliceEnd], dQl[0..sliceEnd], dQr[0..sliceEnd], mesh.normals[0..sliceEnd]);

		// integrate fluxes and update dQdt

		foreach(i; mesh.interiorCells)
		{
			auto initLen = mesh.cells[i].fluxMultiplier[0]*mesh.edges[mesh.cells[i].edges[0]].len;
			auto initFluxIdx = mesh.cells[i].edges[0];
			dQdt[i] = initLen*fluxes[initFluxIdx];

			// integrate fluxes over cell edges
			for(uint j = 1; j < mesh.cells[i].nEdges; j++)
			{
				auto len = mesh.cells[i].fluxMultiplier[j]*mesh.edges[mesh.cells[i].edges[j]].len;
				auto fluxIdx = mesh.cells[i].edges[j];
				dQdt[i] += len*fluxes[fluxIdx];
			}

			dQdt[i] *= -(1.0/mesh.cells[i].area);
		}
	}
}
