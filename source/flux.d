/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.flux;

import std.algorithm;
import std.math : abs, fmax, sqrt;
import std.meta : aliasSeqOf;

import numd.linearalgebra.matrix;

import ebb.euler;

//alias fluxList = aliasSeqOf!(["rusanovFlux", "roeFlux"]);
alias fluxList = aliasSeqOf!(["roeFlux"]);

enum fluxDir
{
	xDir,
	yDir
}

@nogc Vector!dims physicalFlux(size_t dims)(double p, double u, double v, double rho, double e, Vector!2 n)
{
	Vector!dims f = Vector!dims([rho*u, p + rho*u^^2.0, rho*u*v, u*(p + e)]);
	Vector!dims h = Vector!dims([rho*v, rho*u*v, p + rho*v^^2.0, v*(p + e)]);
	return f*n[0] + h*n[1];
	/+
	static if(dir == fluxDir.xDir)
	{
		return Vector!dims([rho*u, p + rho*u^^2.0, rho*u*v, u*(p + e)]);
	}
	else static if(dir == fluxDir.yDir)
	{
		return Vector!dims([rho*v, rho*u*v, p + rho*v^^2.0, v*(p + e)]);
	}
	+/
}
/+
Vector!dims physicalFlux(size_t dims)(Vector!dims u)
{
	double rho = u[0];
	double v = u[1]/rho;
	double e = u[2];
	double p = (gamma - 1)*(e - 0.5*rho*v^^2);
	return Vector!dims([rho*v, p + rho*v^^2.0, v*(p + e)]);
}
+/
//massFlux(rhoL, aL, pL, pm1);
double massFlux(double rho, double a, double pQ, double pM)
{
	immutable double gam1 = 0.5*(gamma + 1.0)/gamma;
	immutable double gam2 = 0.5*(gamma - 1.0)/gamma;
	
	if(pM/pQ >= (1.0 - 1.0e-15))
	{
		return (rho*a)*sqrt(1.0 + gam1*(pM/pQ - 1.0));
	}
	else
	{
		return (rho*a)*gam2*(1.0 - pM/pQ)/(1.0 - (pM/pQ)^^gam2);
	}
}

    //call sonic(vR,aR,PR,vm,amR,vR+aR,SmR)
Vector!3 sonic(double uR, double aR, double pR, double uM, double amR, double waveSpeed1, double waveSpeed2)
{
	Vector!3 u;
	
	double r1 = waveSpeed2/(waveSpeed2 - waveSpeed1);
	double r2 = -waveSpeed1/(waveSpeed2 - waveSpeed1);
	double uS = r1*uR + r2*uM;
	double aS = r1*aR + r2*amR;
	double pS = (aS/aR)^^(2.0*gamma/(gamma - 1.0));
	double rhoS = gamma*pS/(aS^^2.0);
	
	u[0] = rhoS;
	u[1] = uS*rhoS;
	u[2] = pS/(gamma - 1.0) + 0.5*rhoS*uS^^2.0;
	
	return u;
}

Vector!dims godunovFlux(size_t dims)(Vector!dims uL, Vector!dims uR, double dt, double dx)
{
	alias Vec = Vector!dims;
	alias Mat = Matrix!(dims, dims);
	
	Vec flux;
	
	// Left state variables
	double rhoL = uL[0];
	double vL = uL[1]/rhoL;
	double pL = (gamma - 1)*(uL[2] - 0.5*rhoL*vL^^2);
	double aL = sqrt(gamma*(pL/rhoL));
	double hL = (uL[2] + pL)/rhoL;
	
	// Right state variables
	double rhoR = uR[0];
	double vR = uR[1]/rhoR;
	double pR = (gamma - 1)*(uR[2] - 0.5*rhoR*vR^^2);
	double aR = sqrt(gamma*(pR/rhoR));
	double hR = (uR[2] + pR)/rhoR;
	
	// check for sonic/super sonic fluxes
	if(vL/aL >= 1.0)
	{
		flux = physicalFlux!dims(pL, vL, rhoL, uL[2]);
		return flux;
	}
	else if(vR/aR <= -1.0)
	{
		flux = physicalFlux!dims(pR, vR, rhoR, uR[2]);
		return flux;
	}
	
	double pm1 = ((0.5*(vL-vR)*(gamma - 1) + aL + aR)/( aL*pL^^(0.5*(gamma - 1)/gamma) + aR*pR^^(0.5*(gamma - 1)/gamma)))^^(2.0*gamma/(gamma - 1));
	double pm2;
	double mR;
	double mL;
	int k = 0;
	import std.math : approxEqual;
	do
	{
		mL = massFlux(rhoL, aL, pL, pm1);
		mR = massFlux(rhoR, aR, pR, pm1);
		
		pm2 = (mL*pR + mR*pL - mL*mR*(vR - vL))/(mL + mR);
		if(k > 100)
		{
			writeln("shit aint right: pm1 = ", pm1, " pm2 = ", pm2);
		}
		k++;
		pm1 = pm2;
	} while(!approxEqual(pm1, pm2));
	
	//writeln("done");
	mL = massFlux(rhoL, aL, pL, pm1);
	mR = massFlux(rhoR, aR, pR, pm1);
	double vM = (mL*vL + mR*vR - (pR - pL))/(mL + mR);
	immutable double gam2 = 0.5*(gamma - 1.0)/gamma;
	immutable double gam = (gamma + 1.0)/(gamma - 1.0);
	double rhomL;
	double rhomR;
	import std.meta : aliasSeqOf;
	foreach(l; aliasSeqOf!(['L', 'R']))
	{
		immutable pStr = "p"~l;
		immutable rhomStr = "rhom"~l;
		immutable rhoStr = "rho"~l;
		if(pm2/mixin(pStr) >= 1.0)
		{
			//writeln("rho1");
			mixin(rhomStr) = mixin(rhoStr)*(1.0 + gam*pm2/mixin(pStr))/(gam + pm2/mixin(pStr));
		}
		else
		{
			//writeln("rho2");
			mixin(rhomStr) = mixin(rhoStr)*(pm2/mixin(pStr))^^(1.0/gamma);
		}
	}
	
	double rhomI;
	if(vM >= 0.0)
	{
		//writeln("L");
		rhomI = rhomL;
	}
	else
	{
		//writeln("R");
		rhomI = rhomR;
	}
	
	double amL = sqrt(gamma*pm2/rhomL);
	double amR = sqrt(gamma*pm2/rhomR);
	double waveSpeedL = vM - amL;
	double waveSpeedR = vM + amR;
	
	Vec u;
	
	u[0] = rhomI;
	if((waveSpeedL <= 0.0) && (waveSpeedR >= 0.0))
	{
		//writeln("no sonic");
		u[1] = rhomI*vM;
		u[2] = pm2/(gamma - 1.0) + 0.5*rhomI*vM^^2.0;
	}
	else if((waveSpeedL > 0.0) && ((vL - aL) < 0.0))
	{
		//writeln("sonic1");
		u = sonic(vL, aL, pL, vM, amL, vL - aL, waveSpeedL);
	}
	else if((waveSpeedR < 0.0) && ((vR + aR) > 0.0))
	{
		//writeln("sonic2");
		u = sonic(vR, aR, pR, vM, amR, vR + aR, waveSpeedR);
	}
	
	flux = physicalFlux!dims(u);
	return flux;
}

// Vector!dims rusanovFlux(size_t dims)(Vector!dims qL, Vector!dims qR, Vector!(dims - 2) n)
Vector!dims rusanovFlux(size_t dims, int kx, int ky)(Vector!dims qL, Vector!dims qR, double dt, double dx)
{
	alias Vec = Vector!dims;
	alias Mat = Matrix!(dims, dims);
	
	Vec flux;

	// Left state variables
	double rhoL = (qL[0]);
	double uL = qL[1]/rhoL;
	double vL = qL[2]/rhoL;
	double pL = ((gamma - 1)*(qL[3] - 0.5*rhoL*(uL^^2.0 + vL^^2.0)));
	//double aL = sqrt(gamma*(pL/rhoL));
	double hL = (qL[3] + pL)/rhoL;
	
	// Right state variables
	double rhoR = (qR[0]);
	double uR = qR[1]/rhoR;
	double vR = qR[2]/rhoR;
	double pR = ((gamma - 1)*(qR[3] - 0.5*rhoR*(uR^^2.0 + vR^^2.0)));
	//double aR = sqrt(gamma*(pR/rhoR));
	double hR = (qR[3] + pR)/rhoR;

	// average values
	double rho = sqrt(rhoL*rhoR);
	double u = (sqrt((rhoL))*uL + sqrt((rhoR))*uR)/(sqrt((rhoL)) + sqrt((rhoR)));
	double v = (sqrt((rhoL))*vL + sqrt((rhoR))*vR)/(sqrt((rhoL)) + sqrt((rhoR)));
	double h = (sqrt((rhoL))*hL + sqrt((rhoR))*hR)/(sqrt((rhoL)) + sqrt((rhoR)));
	double a = sqrt( (gamma - 1.0)*(h - 0.5*(v^^2.0 + u^^2.0)));
	static if(kx == 1)
	{
		Vec Fl = physicalFlux!(fluxDir.xDir, dims)(pL, uL, vL, rhoL, qL[3]);
		Vec Fr = physicalFlux!(fluxDir.xDir, dims)(pR, uR, vR, rhoR, qR[3]);
		flux = 0.5*((Fl + Fr) - (abs(u)+a)*(qR - qL));
	}
	else static if(ky == 1)
	{
		Vec Fl = physicalFlux!(fluxDir.yDir, dims)(pL, uL, vL, rhoL, qL[3]);
		Vec Fr = physicalFlux!(fluxDir.yDir, dims)(pR, uR, vR, rhoR, qR[3]);
		flux = 0.5*((Fl + Fr) - (abs(v)+a)*(qR - qL));
	}
	
	return flux;
}

/++
	Uses my matrix operations
+/
@nogc Vector!dims roeFlux(size_t dims)(Vector!dims qL, Vector!dims qR, Vector!(dims - 2) n, ref double sMax)
{
	alias Vec = Vector!dims;
	alias Mat = Matrix!(dims, dims);
	
	Vec flux;

	// Left state variables
	double rhoL = (qL[0]);
	double uL = qL[1]/rhoL;
	double vL = qL[2]/rhoL;
	auto rhoVelL = rhoL*Vector!2(uL, vL);

	double pL = ((gamma - 1)*(qL[3] - 0.5*rhoL*(uL^^2.0 + vL^^2.0)));
	double aL = sqrt(gamma*(pL/rhoL));
	double hL = (qL[3] + pL)/rhoL;
	
	// Right state variables
	double rhoR = (qR[0]);
	double uR = qR[1]/rhoR;
	double vR = qR[2]/rhoR;
	auto rhoVelR = rhoR*Vector!2(uR, vR);
	double pR = ((gamma - 1)*(qR[3] - 0.5*rhoR*(uR^^2.0 + vR^^2.0)));
	double aR = sqrt(gamma*(pR/rhoR));
	double hR = (qR[3] + pR)/rhoR;

	// average values
	double rho = sqrt(rhoL*rhoR);
	auto vel = Vector!2(0);
	vel[0] = (sqrt((rhoL))*uL + sqrt((rhoR))*uR)/(sqrt((rhoL)) + sqrt((rhoR)));
	vel[1] = (sqrt((rhoL))*vL + sqrt((rhoR))*vR)/(sqrt((rhoL)) + sqrt((rhoR)));
	double u = vel.dot(n);

	double h = (sqrt((rhoL))*hL + sqrt((rhoR))*hR)/(sqrt((rhoL)) + sqrt((rhoR)));
	double a = sqrt( (gamma - 1.0)*(h - 0.5*(vel.magnitude^^2.0)));

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
		lam1 = (eps^^2 + lam1^^2)/(2.0*eps);
	}
	if(lam2 < eps)
	{
		lam2 = (eps^^2 + lam2^^2)/(2.0*eps);
	}
	if(lam3 < eps)
	{
		lam3 = (eps^^2 + lam3^^2)/(2.0*eps);
	}
	if(lam4 < eps)
	{
		lam4 = (eps^^2 + lam4^^2)/(2.0*eps);
	}

	sMax = max(abs(lam1), abs(lam2), abs(lam3), abs(lam4));

	double s1 = 0.5*(lam1 + lam2);
	double s2 = 0.5*(lam1 - lam2);

	double q = sqrt(vel[0]^^2 + vel[1]^^2);
	double G1 = (gamma - 1)*((q^^2/2)*dRho - vel.dot(dRhoV) + dRhoE);
	double G2 = -u*dRho + dRhoV.dot(n);

	double C1 = (G1/a^^2)*(s1 - lam3) + G2/a*s2;
	double C2 = (G1/a)*s2 + (s1 - lam3)*G2;

	auto fL = physicalFlux!dims(pL, uL, vL, rhoL, qL[3], n);
	auto fR = physicalFlux!dims(pR, uR, vR, rhoR, qR[3], n);

	flux = 0.5*(fL + fR) - 0.5*Vector!4(lam3*dRho + C1, lam3*(dRhoV[0]) + C1*vel[0] + C2*n[0], lam3*(dRhoV[1]) + C1*vel[1] + C2*n[1], lam3*dRhoE + C1*h + C2*u);
	return flux;
}
