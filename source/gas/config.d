/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.gas.config;

import ebb.euler;
import ebb.exception;
import ebb.gas.flux;
import ebb.mesh;
import ebb.solve;

import numd.linearalgebra.matrix;

import std.algorithm;
import std.conv;
import std.json;
import std.math;
import std.stdio;

struct GasPhysicalConfig
{
	immutable double gamma;
	immutable double mu;
	immutable double Pr;
	immutable double R;

	private Vector!2[uint] forces;
	/+this()
	{
		force = Vector!2(0);
	}

	@nogc Vector!2 getBoundaryValues()
	{
		auto tmp = force;
		force = Vector!2(0);
		return tmp;
	}+/

	@nogc bool needEdgeGradients()
	{
		return mu != 0.0;
	}

	@nogc bool checkPhysics(Vec)(ref Vec q)
	{
		return getPressure(q, gamma) < 0;
	}

	@nogc void updateGhostCells(Vec)(ref Vec[] Q, ref immutable Mesh mesh)
	{
		foreach(i; mesh.ghostCells)
		{
			auto edge = mesh.edges[mesh.cells[i].edges[0]];
			switch(edge.boundaryType)
				with(BoundaryType)
			{
				case Dirichlet:
					auto centroid = mesh.cells[i].centroid;
					//Q[i] = config.boundaries[edge.bIdx].dFunc(centroid[0], centroid[1]);
					break;
				case FullState:
					// We shouldn't be updating on a fullstate boundary.
					//Q[edge.cellIdxR] = Q[edge.cellIdxL];
					break;
				case InviscidWall:
					auto cellIdx = edge.cellIdxL;
					auto v = Vector!2(Q[cellIdx][1]/Q[cellIdx][0], Q[cellIdx][2]/Q[cellIdx][0]);
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

					auto cellIdx2 = edge.cellIdxR;
					Q[cellIdx2][0] = Q[cellIdx][0];
					Q[cellIdx2][1] = Q[cellIdx][0]*newV[0];
					Q[cellIdx2][2] = Q[cellIdx][0]*newV[1];
					Q[cellIdx2][3] = Q[cellIdx][3];

					break;
				case ViscousWall:
					auto cellIdx = edge.cellIdxL;

					auto cellIdx2 = edge.cellIdxR;
					Q[cellIdx2][0] = Q[cellIdx][0];
					Q[cellIdx2][1] = -Q[cellIdx][1];
					Q[cellIdx2][2] = -Q[cellIdx][2];
					Q[cellIdx2][3] = Q[cellIdx][3];
					break;
				case ConstPressure:
					auto cellIdx = edge.cellIdxL;
					auto cellIdx2 = edge.cellIdxR;

					double rhop = Q[cellIdx][0];
					double up = Q[cellIdx][1]/Q[cellIdx][0];
					double vp = Q[cellIdx][2]/Q[cellIdx][0];
					auto Vp = Vector!2(up, vp);
					double pp = getPressure(Q[cellIdx], gamma); // lol
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

					Q[cellIdx2][0] = rhob;
					Q[cellIdx2][1] = rhob*Vbg[0];
					Q[cellIdx2][2] = rhob*Vbg[1];
					Q[cellIdx2][3] = rEb;

					break;

				case Symmetry:
					auto cellIdx = edge.cellIdxL;
					auto v = Vector!2(Q[cellIdx][1]/Q[cellIdx][0], Q[cellIdx][2]/Q[cellIdx][0]);
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

					auto cellIdx2 = edge.cellIdxR;
					Q[cellIdx2][0] = Q[cellIdx][0];
					Q[cellIdx2][1] = Q[cellIdx][0]*newV[0];
					Q[cellIdx2][2] = Q[cellIdx][0]*newV[1];
					Q[cellIdx2][3] = Q[cellIdx][3];
					break;
				
				case TempPresInflow:
					auto cellIdx = edge.cellIdxL;
					auto cellIdx2 = edge.cellIdxR;
					
					double Tt = edge.bData[0];
					double Pt = edge.bData[1];
					double aoa = edge.bData[2] * PI/180;

					double rhop = Q[cellIdx][0];
					double up = Q[cellIdx][1]/Q[cellIdx][0];
					double vp = Q[cellIdx][2]/Q[cellIdx][0];
					auto Vp = Vector!2(up, vp);
					double pp = getPressure(Q[cellIdx], gamma); // lol
					double ap = sqrt(gamma*(pp/rhop));
					double Jp = up + (2.0*ap)/(gamma - 1);

					auto n_in = Vector!2(cos(aoa), sin(aoa));
					auto d = n_in.dot(edge.normal);
					double a = gamma*R*Tt*d^^2.0 - 0.5*(gamma - 1)*Jp^^2.0;
					double b = (4.0*gamma*R*Tt*d)/(gamma - 1);
					double c = (4.0*gamma*R*Tt)/(gamma - 1)^^2.0 - Jp^^2.0;

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
						printf("%f\n", (4.0*gamma*R*Tt)/(gamma - 1)^^2.0);
						enforce(false, "Both computed mach numbers are negative");
					}

					double Tb = Tt/(1 + 0.5*(gamma - 1)*Mb^^2.0);
					double pb = Pt*(Tb/Tt)^^(gamma/(gamma - 1));
					double rhob = pb/(R*Tb);
					double ab = sqrt(gamma*(pb/rhob));
					auto vb = Mb*ab*n_in;
					double rEb = pb/(gamma - 1) + 0.5*rhob*vb.magnitude^^2.0;

					Q[cellIdx2][0] = rhob;
					Q[cellIdx2][0] = rhob*vb[0];
					Q[cellIdx2][0] = rhob*vb[1];
					Q[cellIdx2][0] = rEb;

					break;

				default:
					enforce(false, "Unsupported boundary type");
					break;
			}
		}
	}

	@nogc FluxResult!size updateBoundaryEdgeFluxes(size_t size, FluxFunction)(
			ref Vector!size qR,
			immutable Vector!size qL,
			immutable Matrix!(size,size-2) dQl,
			immutable Matrix!(size,size-2) dQr,
			immutable(Vector!size[]) Q,
			ref immutable Edge edge,
			ref immutable Mesh mesh,
			FluxFunction fluxFunction
		)
	{
		switch(edge.boundaryType)
			with(BoundaryType)
		{
			case Dirichlet:
				/*
				double x = mesh.cells[mesh.edges[i].cellIdxR].centroid[0];
				double y = mesh.cells[mesh.edges[i].cellIdxR].centroid[1];

				mesh.cells[mesh.edges[i].cellIdxR].gradient = config.boundaries[mesh.edges[i].bIdx].dFuncDerivative(x, y);
				auto mid = mesh.edges[i].mid;
				qR[i] = config.boundaries[mesh.edges[i].bIdx].dFunc(mid[0], mid[1]);
				*/
				enforce(false, "Dirichlet boudaries are not supported yet");
				goto case FullState;

			case FullState:
				return fluxFunction(qL, qR, dQl, dQr, edge.normal);

			case InviscidWall:
				Vector!2 velP = (1/qL[0])*Vector!2(qL[1], qL[2]);
				auto vel = (velP - (velP.dot(edge.normal))*edge.normal).magnitude;
				double p = (gamma - 1)*(qL[3] - 0.5*qL[0]*vel^^2);
				double a = sqrt(gamma*(p/qL[0]));
				if(p < 0)
				{
					p = 1.0e-12;
					//printf("pressure less than 0 at wall\n");
				}
				qR = Q[edge.cellIdxR];
				auto flux = Vector!size(0, p*edge.normal[0], p*edge.normal[1], 0);
				auto force = edge.bTag in forces;
				if(force)
				{
					(*force)[0] += edge.len*flux[1];
					(*force)[1] += edge.len*flux[2];
				}
				return fluxResult(flux, abs(a));

			case ViscousWall:
				auto vel = 0.0;
				double p = (gamma - 1)*(qL[3] - 0.5*qL[0]*vel^^2);
				double a = sqrt(gamma*(p/qL[0]));
				if(p < 0)
				{
					p = 1.0e-12;
					//printf("pressure less than 0 at wall\n");
				}

				qR = qL;
				qR[1] *= -1;
				qR[2] *= -1;
				
				auto qAve = 0.5*(qL + qR);

				qAve[1] = 0.0;
				qAve[2] = 0.0;

				
				Matrix!(size,size-2) dqAve = dQl;
				// modify energy term to be an adiabatic wall.
				// TODO: Generalize this. Equations in EbbCFD book 1 page 5
				dqAve[3,0] = (qAve[3]/qAve[0])*dqAve[0,0];
				dqAve[3,1] = (qAve[3]/qAve[0])*dqAve[0,1];

				auto Fv = diffusiveFlux!size(qAve, dqAve, edge.normal, this);
				auto flux = Vector!size(0, p*edge.normal[0], p*edge.normal[1], 0) - Fv.flux;
				// TODO: Check to see if this is correct. There should be some force
				// contribution from the viscous flux... I think.
				auto force = edge.bTag in forces;
				if(force)
				{
					(*force)[0] += edge.len*flux[1];
					(*force)[1] += edge.len*flux[2];
				}
				return fluxResult(flux, abs(a));

			case ConstPressure:
				auto cellIdx2 = edge.cellIdxR;

				auto p = getPressure(Q[cellIdx2], gamma);
				auto u = Q[cellIdx2][1]/Q[cellIdx2][0];
				auto v = Q[cellIdx2][2]/Q[cellIdx2][0];
				double a = sqrt(gamma*(p/Q[cellIdx2][0]));
				qR = Q[cellIdx2];
				return fluxResult(convectiveFlux!4(p, u, v, Q[cellIdx2][0], Q[cellIdx2][3], edge.normal), abs(a));
				
			case Symmetry:
				auto v = Vector!2(qL[1]/qL[0], qL[2]/qL[0]);
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

				qR[0] = qL[0];
				qR[1] = qL[0]*newV[0];
				qR[2] = qL[0]*newV[1];
				qR[3] = qL[3];

				if(mu != 0)
				{
					auto cellIdx = edge.cellIdxL;
					auto gradp = dQl;
					auto nOuter = edge.normal*edge.normal.transpose;
					auto A = Matrix!(size,size).identity;
					auto V = Matrix!(size-2,size-2).identity - nOuter;
					A[1,1] = V[0,0];
					A[1,2] = V[0,1];
					A[2,1] = V[1,0];
					A[2,2] = V[1,1];
					auto tmp1 = A*gradp;
					auto tmp2 = Matrix!(size-2,size-2).identity - 2*nOuter;
					auto tmp3 = tmp1*tmp2;
					auto bGrad = gradp*nOuter + tmp3;

					return fluxFunction(qL, qR, bGrad, bGrad, edge.normal);
				}
				else
				{
					return fluxFunction(qL, qR, dQl, dQr, edge.normal);
				}

			case TempPresInflow:
				qR = Q[edge.cellIdxR];
				return fluxFunction(qL, qR, dQl, dQl, edge.normal);

			default:
				enforce(false, "Unsupported boundary type");
				return fluxResult(Vector!size(0), -double.infinity);
		}
	}
}
