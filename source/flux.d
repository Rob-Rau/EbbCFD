/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.flux;

import std.algorithm;
import std.math : abs, fmax, sqrt;
import std.meta : aliasSeqOf;

import numd.linearalgebra.matrix;

import ebb.euler;

alias fluxList = aliasSeqOf!(["rusanovFlux", "roeFlux"]);

/++
	Examples:
	auto fL = convectiveFlux!dims(pL, uL, vL, rhoL, qL[3], n);
	auto fR = convectiveFlux!dims(pR, uR, vR, rhoR, qR[3], n);
+/
@nogc Vector!dims convectiveFlux(size_t dims)(double p, double u, double v, double rho, double e, Vector!2 n)
{
	Vector!dims f = Vector!dims([rho*u, p + rho*u^^2.0, rho*u*v, u*(p + e)]);
	Vector!dims h = Vector!dims([rho*v, rho*u*v, p + rho*v^^2.0, v*(p + e)]);
	return f*n[0] + h*n[1];
}

/++
+/
@nogc Vector!dims diffusiveFlux(size_t dims)(double Pr, double mu, Vector!dims q, Matrix!(dims,2) dq, Vector!2 n)
{
	double rho = q[0];
	double u = q[1]/q[0];
	double v = q[2]/q[0];
	double E = q[3]/q[0];

	auto dqdx = Vector!dims(dq[0,0], dq[1,0], dq[2,0], dq[3,0]);
	auto dqdy = Vector!dims(dq[0,1], dq[1,1], dq[2,1], dq[3,1]);

	double tmp = (-4.0/3.0 + gamma/Pr)*u^^2.0 + (-1.0 + gamma/Pr)*v^^2.0 - (gamma/Pr)*E;
	auto A00 = (mu/rho)*Matrix!(dims, dims)(0, 0, 0, 0, 
											-4.0/3.0*u, 4.0/3.0, 0, 0,
											-v, 0, 1, 0,
											tmp, (4.0/3.0 - gamma/Pr)*u, (1.0 - gamma/Pr)*v, gamma/Pr);

	auto A01 = (mu/rho)*Matrix!(dims, dims)(0, 0, 0, 0,
											2.0/3.0*v, 0, -2.0/3.0, 0,
											-u, 1, 0, 0,
											(2.0/3.0 - 1.0)*u*v, v, -2.0/3.0*u, 0);

	auto G0 = A00*dqdx + A01*dqdy;

	auto A10 = (mu/rho)*Matrix!(dims, dims)(0, 0, 0, 0,
											-v, 0, 1, 0,
											2.0/3.0*u, 0, -2.0/3.0, 0,
											(2.0/3.0 - 1.0)*u*v, -2.0/3.0*v, u, 0);

	tmp = (-1.0 + gamma/Pr)*u^^2.0 + (-4.0/3.0 + gamma/Pr)*v^^2.0 - (gamma/Pr)*E;
	auto A11 = (mu/rho)*Matrix!(dims, dims)(0, 0, 0, 0, 
											-u, 1, 0, 0,
											-4.0/3.0*v, 0, 4.0/3.0, 0,
											tmp, (1.0 - gamma/Pr)*u, (4.0/3.0 - gamma/Pr)*v, gamma/Pr);

	auto G1 = A10*dqdx + A11*dqdy;
	return G0*n[0] + G1*n[1];
}

@nogc Vector!dims rusanovFlux(size_t dims)(Vector!dims qL, Vector!dims qR, Vector!(dims - 2) n, ref double sMax)
{
	// Left state variables
	immutable double rhoL = (qL[0]);
	immutable double uL = qL[1]/rhoL;
	immutable double vL = qL[2]/rhoL;
	auto rhoVelL = rhoL*Vector!2(uL, vL);

	immutable double pL = ((gamma - 1)*(qL[3] - 0.5*rhoL*(uL^^2.0 + vL^^2.0)));
	immutable double aL = sqrt(gamma*(pL/rhoL));
	immutable double hL = (qL[3] + pL)/rhoL;
	
	// Right state variables
	immutable double rhoR = (qR[0]);
	immutable double uR = qR[1]/rhoR;
	immutable double vR = qR[2]/rhoR;
	auto rhoVelR = rhoR*Vector!2(uR, vR);
	immutable double pR = ((gamma - 1)*(qR[3] - 0.5*rhoR*(uR^^2.0 + vR^^2.0)));
	immutable double aR = sqrt(gamma*(pR/rhoR));
	immutable double hR = (qR[3] + pR)/rhoR;

	// average values
	immutable double rho = sqrt(rhoL*rhoR);
	auto vel = Vector!2(0);
	vel[0] = (sqrt((rhoL))*uL + sqrt((rhoR))*uR)/(sqrt((rhoL)) + sqrt((rhoR)));
	vel[1] = (sqrt((rhoL))*vL + sqrt((rhoR))*vR)/(sqrt((rhoL)) + sqrt((rhoR)));
	immutable double u = vel.dot(n);
	immutable double p = (sqrt(rhoL)*pL + sqrt(rhoR)*pR)/(sqrt(rhoL) + sqrt(rhoR));
	immutable double h = (sqrt((rhoL))*hL + sqrt((rhoR))*hR)/(sqrt((rhoL)) + sqrt((rhoR)));
	immutable double a = sqrt(gamma*(p/rho));

	auto dRhoV = rhoVelR - rhoVelL;
	auto dRho = rhoR - rhoL;
	auto dRhoE = qR[3] - qL[3];

	double lam1 = abs(u + a);
	double lam2 = abs(u - a);
	double lam3 = abs(u);
	double lam4 = abs(u);

	immutable double eps = 0.01*a;
	if(lam1 < eps)
	{
		lam1 = (eps^^2.0 + lam1^^2.0)/(2.0*eps);
	}
	if(lam2 < eps)
	{
		lam2 = (eps^^2.0 + lam2^^2.0)/(2.0*eps);
	}
	if(lam3 < eps)
	{
		lam3 = (eps^^2.0 + lam3^^2.0)/(2.0*eps);
	}
	if(lam4 < eps)
	{
		lam4 = (eps^^2.0 + lam4^^2.0)/(2.0*eps);
	}

	sMax = max(abs(lam1), abs(lam2), abs(lam3), abs(lam4));

	immutable double s1 = 0.5*(sMax + sMax);
	immutable double s2 = 0.5*(sMax - sMax);

	immutable double q = sqrt(vel[0]^^2 + vel[1]^^2);
	immutable double G1 = (gamma - 1)*((q^^2/2)*dRho - vel.dot(dRhoV) + dRhoE);
	immutable double G2 = -u*dRho + dRhoV.dot(n);

	immutable double C1 = (G1/a^^2)*(s1 - sMax) + G2/a*s2;
	immutable double C2 = (G1/a)*s2 + (s1 - sMax)*G2;

	auto fL = convectiveFlux!dims(pL, uL, vL, rhoL, qL[3], n);
	auto fR = convectiveFlux!dims(pR, uR, vR, rhoR, qR[3], n);

	auto flux = 0.5*(fL + fR) - 0.5*Vector!4(sMax*dRho + C1, sMax*(dRhoV[0]) + C1*vel[0] + C2*n[0], sMax*(dRhoV[1]) + C1*vel[1] + C2*n[1], sMax*dRhoE + C1*h + C2*u);
	return flux;
}

unittest 
{
	import ebb.manufacturedsolution;
	import std.stdio : writeln;
	immutable double eps = 1.0e-5;

	double x1 = 1;
	double y1 = 1;
	double x2 = 2;
	double y2 = 1;

	auto n = Vector!2(1, 0);
	double sMax;

	auto qL = solution(x1, y1);
	auto qR = solution(x2, y2);
	auto F = roeFlux!4(qL, qR, n, sMax);

	writeln(abs(F[0] - 0.087503926254201));
	writeln(abs(F[1] - 1.039732984117400));
	writeln(abs(F[2] - 0.008215342318582));
	writeln(abs(F[3] - 0.367160441647966));
	writeln;
	//writeln(F.ToString);

	x1 = 1;
	y1 = 1;
	x2 = 1;
	y2 = 2;

	n = Vector!2(0, 1);
	qL = solution(x1, y1);
	qR = solution(x2, y2);
	F = roeFlux!4(qL, qR, n, sMax);

	writeln(abs(F[0] - 0.117173193168260));
	writeln(abs(F[1] - 0.010523092695679));
	writeln(abs(F[2] - 1.005792391773638));
	writeln(abs(F[3] - 0.482469106067425));
	writeln;
	//writeln(F.ToString);

	x1 = 1;
	y1 = 1;
	x2 = 2;
	y2 = 2;

	n = Vector!2(sqrt(2.0)/2.0, sqrt(2.0)/2.0);
	qL = solution(x1, y1);
	qR = solution(x2, y2);
	F = roeFlux!4(qL, qR, n, sMax);

	writeln(abs(F[0] - 0.116317039552757));
	writeln(abs(F[1] - 0.743804395816748));
	writeln(abs(F[2] - 0.743854568709977));
	writeln(abs(F[3] - 0.490010270963920));
	writeln;
	//writeln(F.ToString);

	x1 = 1;
	y1 = 1;
	x2 = 1.1;
	y2 = 1.1;

	n = Vector!2(sqrt(2.0)/2.0, sqrt(2.0)/2.0);
	qL = solution(x1, y1);
	qR = solution(x2, y2);
	F = roeFlux!4(qL, qR, n, sMax);

	writeln(abs(F[0] - 0.115207589649517));
	writeln(abs(F[1] - 0.744728044864156));
	writeln(abs(F[2] - 0.744751077184082));
	writeln(abs(F[3] - 0.485506505793277));
	writeln;
	//writeln(F.ToString);

}

/++
	Uses my matrix operations
+/
@nogc Vector!dims roeFlux(size_t dims)(Vector!dims qL, Vector!dims qR, Vector!(dims - 2) n, ref double sMax)
{
	// Left state variables
	immutable double rhoL = (qL[0]);
	immutable double uL = qL[1]/rhoL;
	immutable double vL = qL[2]/rhoL;
	auto rhoVelL = rhoL*Vector!2(uL, vL);

	immutable double pL = ((gamma - 1)*(qL[3] - 0.5*rhoL*(uL^^2.0 + vL^^2.0)));
	immutable double aL = sqrt(gamma*(pL/rhoL));
	immutable double hL = (qL[3] + pL)/rhoL;
	
	// Right state variables
	immutable double rhoR = (qR[0]);
	immutable double uR = qR[1]/rhoR;
	immutable double vR = qR[2]/rhoR;
	auto rhoVelR = rhoR*Vector!2(uR, vR);
	immutable double pR = ((gamma - 1)*(qR[3] - 0.5*rhoR*(uR^^2.0 + vR^^2.0)));
	immutable double aR = sqrt(gamma*(pR/rhoR));
	immutable double hR = (qR[3] + pR)/rhoR;

	immutable double di = sqrt(rhoR/rhoL);
	immutable double d1 = 1.0/(1.0 + di);
	// average values
	auto vel = Vector!2(0);
	vel[0] = (di*uR + uL)*d1;
	vel[1] = (di*vR + vL)*d1;
	immutable double u = vel.dot(n);
	immutable double h = (di*hR + hL)*d1;
	immutable double a = sqrt((gamma - 1.0)*(h - 0.5*(vel[0]*vel[0] + vel[1]*vel[1])));

	auto dRhoV = rhoVelR - rhoVelL;
	auto dRho = rhoR - rhoL;
	auto dRhoE = qR[3] - qL[3];

	double lam1 = abs(u + a);
	double lam2 = abs(u - a);
	double lam3 = abs(u);
	double lam4 = abs(u);

	immutable double eps = 0.01*a;
	if(lam1 < eps)
	{
		lam1 = (eps^^2.0 + lam1^^2.0)/(2.0*eps);
	}
	if(lam2 < eps)
	{
		lam2 = (eps^^2.0 + lam2^^2.0)/(2.0*eps);
	}
	if(lam3 < eps)
	{
		lam3 = (eps^^2.0 + lam3^^2.0)/(2.0*eps);
	}
	if(lam4 < eps)
	{
		lam4 = (eps^^2.0 + lam4^^2.0)/(2.0*eps);
	}

	sMax = max(abs(lam1), abs(lam2), abs(lam3), abs(lam4));

	immutable double s1 = 0.5*(lam1 + lam2);
	immutable double s2 = 0.5*(lam1 - lam2);

	immutable double q = sqrt(vel[0]^^2.0 + vel[1]^^2.0);
	immutable double G1 = (gamma - 1)*(((q^^2.0)/2.0)*dRho - vel.dot(dRhoV) + dRhoE);
	immutable double G2 = -u*dRho + dRhoV.dot(n);

	immutable double C1 = (G1/a^^2)*(s1 - lam3) + G2/a*s2;
	immutable double C2 = (G1/a)*s2 + (s1 - lam3)*G2;

	auto fL = convectiveFlux!dims(pL, uL, vL, rhoL, qL[3], n);
	auto fR = convectiveFlux!dims(pR, uR, vR, rhoR, qR[3], n);

	auto flux = 0.5*(fL + fR) - 0.5*Vector!4(lam3*dRho + C1,
											 lam3*dRhoV[0] + C1*vel[0] + C2*n[0],
											 lam3*dRhoV[1] + C1*vel[1] + C2*n[1],
											 lam3*dRhoE + C1*h + C2*u);

	return flux;
}
