/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.solver;

import core.atomic;
import core.stdc.stdio : fopen, fwrite, fopen, printf, snprintf;
import core.sys.posix.signal;

import std.algorithm : canFind, countUntil, min, max, reduce, sum;
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

static shared bool interrupted = false;

@nogc @system nothrow extern(C) void handle(int sig)
{
	printf("Signal received\n");
	atomicStore(interrupted, true);
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

@nogc void runIntegrator(alias setup, alias solver, alias integrator)(ref UMesh2 mesh, Config config, string saveFile, SolverException ex)
{
	import std.experimental.allocator.mallocator : Mallocator;
	import std.bitmanip : write;

	double residRhoLast = double.infinity;

	uint residRhoIncIters = 0;

	uint iterations = 0;
	uint saveItr = 0;
	double t = 0;
	double dt = config.dt;

	auto rotMat = Matrix!(2, 2)(cos(config.ic[1] * (PI/180)), -sin(config.ic[1] * (PI/180)), sin(config.ic[1] * (PI/180)), cos(config.ic[1] * (PI/180)));
	auto ld = Vector!2(0);

	immutable uint buffSize = 3*1024*1024*double.sizeof;
	size_t buffPos = 0;

	FILE* forceFile;
	if(mesh.mpiRank == 0)
	{
		forceFile = fopen("boundaryForces.frc", "wb");
	}

	double residRho = 0.0;
	double residU = 0.0;
	double residV = 0.0;
	double residE = 0.0;
	double residMax = 0.0;

	double[] lastRho = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] thisRho = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] lastU = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] thisU = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] lastV = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] thisV = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] lastE = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] thisE = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	double[] tmp = cast(double[])Mallocator.instance.allocate(mesh.cells.length*double.sizeof);
	ubyte[] forceBuffer = cast(ubyte[])Mallocator.instance.allocate(buffSize);
	Vector!4[] R = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);

	lastRho[] = 0.0;
	thisRho[] = 0.0;
	lastU[] = 0.0;
	thisU[] = 0.0;
	lastV[] = 0.0;
	thisV[] = 0.0;
	lastE[] = 0.0;
	thisE[] = 0.0;
	tmp[] = 0.0;

	Vector!4[] q0;
	Vector!4[] q1;
	Vector!4[] q2;
	Vector!4[] dem;
	if(config.aitkenTol > 0)
	{
		q0 = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
		q1 = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
		q2 = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);
		dem = cast(Vector!4[])Mallocator.instance.allocate(mesh.cells.length*Vector!4.sizeof);

		for(uint i = 0; i < mesh.cells.length; i++)
		{
			q0[i] = Vector!4(0);
			q1[i] = Vector!4(0);
			q2[i] = Vector!4(0);
			dem[i] = Vector!4(0);
		}
		scope(exit) Mallocator.instance.deallocate(q0);
		scope(exit) Mallocator.instance.deallocate(q1);
		scope(exit) Mallocator.instance.deallocate(q2);
		scope(exit) Mallocator.instance.deallocate(dem);
	}

	scope(exit) Mallocator.instance.deallocate(forceBuffer);
	scope(exit) Mallocator.instance.deallocate(R);
	scope(exit) Mallocator.instance.deallocate(lastRho);
	scope(exit) Mallocator.instance.deallocate(thisRho);
	scope(exit) Mallocator.instance.deallocate(lastU);
	scope(exit) Mallocator.instance.deallocate(thisU);
	scope(exit) Mallocator.instance.deallocate(lastV);
	scope(exit) Mallocator.instance.deallocate(thisV);
	scope(exit) Mallocator.instance.deallocate(lastE);
	scope(exit) Mallocator.instance.deallocate(thisE);
	scope(exit) Mallocator.instance.deallocate(tmp);

	// let the integrator do any neccessary initialization
	integrator.init(mesh);

	double lastRmax = 0.0;

	// Setup IC's and BC's
	setup(mesh, config, lastRho, lastU, lastV, lastE, t, dt, saveFile, ex);
	MPI_Barrier(mesh.comm);
	while((t < config.tEnd) && !atomicLoad(interrupted))
	{
		double Rmax = 0;
		double newDt = dt;

		if((config.aitkenTol < 0) || ((config.aitkenTol > 0) && (iterations < 500)))
		{
			integrator.step!solver(R, mesh.q, mesh, config, newDt, Rmax, ex);
		}
		else
		{
			q0[] = mesh.q[];

			integrator.step!solver(R, q1, mesh, config, newDt, Rmax, ex);
			newDt = dt;
			mesh.q[] = q1[];

			integrator.step!solver(R, q2, mesh, config, newDt, Rmax, ex);

			for(uint i = 0; i < mesh.cells.length; i++)
			{
				for(uint k = 0; k < 4; k++)
				{
					printf("%f, %f, %f\n", q0[i][k], q1[i][k], q2[i][k]);
				}
				printf("\n");
			}
		}

		if(config.dynamicDt)
		{
			dt = newDt;
		}

		foreach(i; mesh.interiorCells)
		{
			thisRho[i] = mesh.q[i][0];
			thisU[i] = mesh.q[i][1];
			thisV[i] = mesh.q[i][2];
			thisE[i] = mesh.q[i][3];
			
			if(mesh.q[i][0].isNaN || mesh.q[i][1].isNaN || mesh.q[i][2].isNaN || mesh.q[i][3].isNaN)
			{
				printf("p = %f\n", getPressure(mesh.q[i]));
				printf("q = [%f, %f, %f, %f]\n", mesh.q[i][0], mesh.q[i][1], mesh.q[i][2], mesh.q[i][3]);
				printf("cell = %d\n", i);
				printf("iteration = %d\n", iterations);
				printf("dt = %f\n", dt);
				printf("t = %f\n", t);
				ex.msg = "Got nan on cell average value";
				ex.file = __FILE__;
				ex.line = __LINE__;
				throw ex;
			}
		}

		import std.algorithm : sum;

		@nogc double computeResidual(double[] now, double[] last)
		{
			tmp[] = (now[] - last[])^^2;
			double tmpSum = tmp.sum;
			MPI_Allreduce(&tmpSum, &tmpSum, 1, MPI_DOUBLE, MPI_SUM, mesh.comm);
			uint totLen = cast(uint)mesh.interiorCells.length;
			MPI_Allreduce(&totLen, &totLen, 1, MPI_UINT32_T, MPI_SUM, mesh.comm);
			double valSum = now.sum;
			MPI_Allreduce(&valSum, &valSum, 1, MPI_DOUBLE, MPI_SUM, mesh.comm);
			return sqrt(cast(double)totLen*tmpSum)/valSum;
		}

		residRho = computeResidual(thisRho, lastRho);
		residU = computeResidual(thisU, lastU);
		residV = computeResidual(thisV, lastV);
		residE = computeResidual(thisE, lastE);
		/*
		tmp[] = (thisRho[] - lastRho[])^^2;
		residRho = sqrt(mesh.cells.length*tmp.sum)/thisRho.sum;
		tmp[] = (thisU[] - lastU[])^^2;
		residU = sqrt(mesh.cells.length*tmp.sum)/thisU.sum;
		tmp[] = (thisV[] - lastV[])^^2;
		residV = sqrt(mesh.cells.length*tmp.sum)/thisV.sum;
		tmp[] = (thisE[] - lastE[])^^2;
		residE = sqrt(mesh.cells.length*tmp.sum)/thisE.sum;
		*/
		//MPI_Allreduce(&residRho, &residRho, 1, MPI_DOUBLE, MPI_)
		lastRho[] = thisRho[];
		lastU[] = thisU[];
		lastV[] = thisV[];
		lastE[] = thisE[];

		auto f = mesh.computeBoundaryForces(config.forceBoundary);
		if(mesh.mpiRank == 0)
		{
			ld = rotMat*f;

			forceBuffer.write!double(t, &buffPos);
			forceBuffer.write!double(ld[0], &buffPos);
			forceBuffer.write!double(ld[1], &buffPos);
		}

		if(buffPos == buffSize)
		{
			fwrite(forceBuffer.ptr, ubyte.sizeof, buffSize, forceFile);
			buffPos = 0;
		}

		residMax = max(residRho, residU, residV, residE);
		if(residRhoLast < residMax)
		{
			residRhoIncIters++;
		}
		else
		{
			if(residRhoIncIters >= 1)
			{
				residRhoIncIters--;
			}
		}

		if(residRhoIncIters >= 2000)
		{
			config.CFL *= 0.99;
			printf("Max RMS residual hasn't decreased in 2000 iterations, decreasing CFL to %f\n", config.CFL);
			residRhoIncIters = 0;
		}

		residRhoLast = residMax;

		if(iterations % config.plotIter == 0)
		{
			if(mesh.mpiRank == 0)
			{
				printf("lift force = %f\t drag force = %f\t t = %f\n", ld[1], ld[0], t);
				printf("rho_RMS = %.10e\tu_RMS = %.10e\tv_RMS = %.10e\tE_RMS = %.10e\tFlux_R = %.10e\t dt = %10.10f\n", residRho, residU, residV, residE, Rmax, dt);
			}
		}
		
		if(config.saveIter != -1)
		{
			if(iterations % config.saveIter == 0)
			{
				char[512] filename;
				filename[] = 0;
				snprintf(filename.ptr, 512, "save_%d_%d.esln", saveItr, mesh.mpiRank);
				saveSolution(mesh, filename.ptr, t, dt);
				filename[] = 0;
				snprintf(filename.ptr, 512, "save_%d_%d.lsln", saveItr, mesh.mpiRank);
				saveLimits(mesh, filename.ptr, t, dt);
				saveItr++;
			}
		}
		t += dt;
		iterations++;
		lastRmax = Rmax;
	}

	if(mesh.mpiRank == 0)
	{
		if(buffPos != 0)
		{
			fwrite(forceBuffer.ptr, ubyte.sizeof, buffPos, forceFile);
			buffPos = 0;
		}
		fclose(forceFile);
	}

	char[512] filename;
	filename[] = 0;
	snprintf(filename.ptr, 512, "final_%d.esln", mesh.mpiRank);
	saveSolution(mesh, filename.ptr, t, dt);
	snprintf(filename.ptr, 512, "final_%d.lsln", mesh.mpiRank);
	saveLimits(mesh, filename.ptr, t, dt);
	if(mesh.mpiRank == 0)
	{
		printf("lift force = %f\t drag force = %f\t t = %f\n", ld[1], ld[0], t);
		printf("rho_RMS = %.10e\tu_RMS = %.10e\tv_RMS = %.10e\tE_RMS = %.10e\tFlux_R = %.10e\t dt = %10.10f\n", residRho, residU, residV, residE, lastRmax, dt);
	}
}

@nogc void ufvmSetup(ref UMesh2 mesh, Config config, double[] lastRho, double[] lastU, double[] lastV, double[] lastE, ref double t, ref double dt, string saveFile, SolverException ex)
{
	double M = config.ic[0];
	double aoa = config.ic[1] * (PI/180);
	double p = config.ic[2];
	double rho = config.ic[3];
	double a = sqrt(gamma*(p/rho));
	double U = M*a;
	double u = U*cos(aoa);
	double v = U*sin(aoa);

	printf("M = %f\tU = %f\n", M, U);

	if(saveFile == "")
	{
		// Setup initial conditions
		foreach(i; mesh.interiorCells)
		{
			//mesh.q[i] = buildQ(rho, u, v, p);
			mesh.q[i][0] = rho;
			mesh.q[i][1] = M*cos(aoa);
			mesh.q[i][2] = M*sin(aoa);
			mesh.q[i][3] = 1.0/((gamma - 1.0)*gamma) + M^^2.0/2.0;

			lastRho[i] = mesh.q[i][0];
			lastU[i] = mesh.q[i][1];
			lastV[i] = mesh.q[i][2];
			lastE[i] = mesh.q[i][3];
		}
	}
	else
	{
		if(loadSolution(mesh, t, dt, saveFile))
		{
			foreach(i; mesh.interiorCells)
			{
				lastRho[i] = mesh.q[i][0];
				lastU[i] = mesh.q[i][1];
				lastV[i] = mesh.q[i][2];
				lastE[i] = mesh.q[i][3];
			}
		}
		else
		{
			ex.msg = "Failed to load solution";
			ex.line = __LINE__;
			ex.file = __FILE__;
			throw ex;
		}
	}

	// Setup bc's
	for(uint i = 0; i < mesh.bGroups.length; i++)
	{
		@nogc uint findBcIndex(string tag)
		{
			for(uint j = 0; j < config.bTags.length; j++)
			{
				if(config.bTags[j] == tag)
				{
					return j;
				}
			}
			ex.msg = "Could not find match boundary condition tag";
			ex.file = __FILE__;
			ex.line = __LINE__;
			throw ex;
		}

		uint bcIdx = findBcIndex(mesh.bTags[i]);

		for(uint j = 0; j < mesh.bGroups[i].length; j++)
		{
			if(!mesh.edges[mesh.bGroups[i][j]].isBoundary)
			{
				ex.msg = "Edge not boundary edge but should be";
				ex.file = __FILE__;
				ex.line = __LINE__;
				throw ex;
			}

			if(mesh.edges[mesh.bGroups[i][j]].boundaryTag != config.bTags[bcIdx])
			{
				ex.msg = "Incorrect boundary tag";
				ex.file = __FILE__;
				ex.line = __LINE__;
				throw ex;
			}

			M = config.bc[bcIdx][0];
			aoa = config.bc[bcIdx][1] * (PI/180);
			p = config.bc[bcIdx][2];
			rho = config.bc[bcIdx][3];
			a = sqrt(gamma*(p/rho));
			U = M*a;
			u = U*cos(aoa);
			v = U*sin(aoa);
			mesh.edges[mesh.bGroups[i][j]].boundaryType = config.bTypes[bcIdx];

			if(mesh.edges[mesh.bGroups[i][j]].boundaryType == BoundaryType.FullState)
			{
				//mesh.edges[mesh.bGroups[i][j]].q[1] = buildQ(rho, u, v, p);
				//printf("bGroup = %d\n", mesh.bGroups[i][j]);
				mesh.edges[mesh.bGroups[i][j]].q[1][0] = rho;
				mesh.edges[mesh.bGroups[i][j]].q[1][1] = M*cos(aoa);
				mesh.edges[mesh.bGroups[i][j]].q[1][2] = M*sin(aoa);
				mesh.edges[mesh.bGroups[i][j]].q[1][3] = 1.0/((gamma - 1.0)*gamma) + M^^2.0/2.0;
				//printf("i = %d; q0 = %f, %f, %f, %f\n", mesh.bGroups[i][j], mesh.edges[mesh.bGroups[i][j]].q[1][0], mesh.edges[mesh.bGroups[i][j]].q[1][1], mesh.edges[mesh.bGroups[i][j]].q[1][2], mesh.edges[mesh.bGroups[i][j]].q[1][3]);
			}
		}
	}
}

MPI_Datatype vec4dataType;

// Unstructured finite volume solver
@nogc void ufvmSolver(alias S, alias F, size_t dims)(ref Vector!4[] R, ref Vector!4[] q, ref UMesh2 mesh, Config config, ref double newDt, ref double Rmax, bool limit, bool dtUpdate, SolverException ex)
{
	//MPI_Barrier(mesh.comm);

	foreach(commIdx, commEdges; mesh.commEdgeIdx)
	{
		foreach(i, edge; commEdges)
		{
			mesh.stateBuffers[commIdx][i] = q[mesh.edges[edge].cellIdx[0]]; 
		}
		MPI_Send(mesh.stateBuffers[commIdx].ptr, cast(uint)mesh.stateBuffers[commIdx].length, vec4dataType, mesh.commProc[commIdx], mesh.meshTag, mesh.comm);
	}

	for(uint c = 0; c < mesh.commProc.length; c++)
	{
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, mesh.meshTag, mesh.comm, &status);

		// determine which proc this is comming from, will be different order
		// than commProc
		if(mesh.commProc.canFind(status.MPI_SOURCE))
		{
			auto commIdx = mesh.commProc.countUntil(status.MPI_SOURCE);
			auto commEdges = mesh.commEdgeIdx[commIdx];

			MPI_Recv(mesh.stateBuffers[commIdx].ptr, cast(int)mesh.stateBuffers[commIdx].length, vec4dataType, mesh.commProc[commIdx], mesh.meshTag, mesh.comm, &status);

			foreach(i, edge; commEdges)
			{
				q[mesh.edges[edge].cellIdx[1]] = mesh.stateBuffers[commIdx][i];
			}
		}
		else
		{
			printf("Unexpected source message from %d\n", status.MPI_SOURCE);
		}
	}

	// update ghost cell states before we can compute gradients
	foreach(i; mesh.ghostCells)
	{
		auto edge = mesh.edges[mesh.cells[i].edges[0]];
		switch(edge.boundaryType)
			with(BoundaryType)
		{
			case FullState:
				//printf("Cell idx 0 = %d\n", edge.cellIdx[0]);
				//printf("Cell idx 1 = %d\n", edge.cellIdx[1]);
				//printf("i = %d; q0 = %f, %f, %f, %f\n", i, q[edge.cellIdx[0]][0], q[edge.cellIdx[0]][1], q[edge.cellIdx[0]][2], q[edge.cellIdx[0]][3]);
				q[edge.cellIdx[1]] = q[edge.cellIdx[0]];
				break;
			case InviscidWall:
				/+
				auto cellIdx = edge.cellIdx[0];
				auto v = Vector!2(q[cellIdx][1]/q[cellIdx][0], q[cellIdx][2]/q[cellIdx][0]);
				// rotate velocity into edge frame.
				auto localV = edge.rotMat*v;
				// flip normal velocity component;
				localV[1] = -localV[1];
				auto inv = Matrix!(2,2)(0);
				invert(inv, edge.rotMat);
				// rotate back into global frame
				auto newV = inv*localV;
				auto cellIdx2 = edge.cellIdx[1];
				// update ghost cell
				q[cellIdx2] = q[cellIdx];
				// TODO: Fix this for wall boundaries
				+/
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
				
				// reflect
			 	break;

			default:
				printf("HIT DEFAULT!!\n");
				break;
		}
	}

	// Build gradients
	if(config.order > 1)
	{
		foreach(i; mesh.interiorCells)
		{
			Vector!6[4] du;
			for(uint j = 0; j < 4; j++)
			{
				du[j] = Vector!6(0); 
			}
			mesh.cells[i].minQ = q[i];
			mesh.cells[i].maxQ = q[i];
			mesh.cells[i].lim[] = Vector!2(1.0);

			for(uint j = 0; j < mesh.cells[i].nNeighborCells; j++)
			{
				uint idx = mesh.cells[i].neighborCells[j];
				
				for(uint k = 0; k < 4; k++)
				{
					mesh.cells[i].minQ[k] = fmin(mesh.cells[i].minQ[k], q[idx][k]);
					mesh.cells[i].maxQ[k] = fmax(mesh.cells[i].maxQ[k], q[idx][k]);

					du[k][j] = q[idx][k] - q[i][k];
				}
			}

			for(uint j = 0; j < 4; j++)
			{
				mesh.cells[i].gradient[j] = mesh.cells[i].gradMat*du[j];
			}

			if(config.limited && limit)
			{
				//for(uint j = 0; j < mesh.cells[i].nEdges; j++)
				for(uint j = 0; j < mesh.cells[i].nNeighborCells; j++)
				{
					//uint eIdx = mesh.cells[i].edges[j];
					auto qM = q[i];
					auto grad = mesh.cells[i].gradient;
					auto centroid = mesh.cells[i].centroid;
					//auto mid = mesh.edges[eIdx].mid;
					auto mid = mesh.cells[mesh.cells[i].neighborCells[j]].centroid;

					auto dx = mid[0] - centroid[0];
					auto dy = mid[1] - centroid[1];
					auto minQ = mesh.cells[i].minQ;
					auto maxQ = mesh.cells[i].maxQ;
					
					auto qE = Vector!4(0);
					for(uint k = 0; k < dims; k++)
					{
						double s = 1.0;
						qE[k] = qM[k] + (grad[k][0]*dx + grad[k][1]*dy);

						if(qE[k] < minQ[k])
						{
							s = (minQ[k] - qM[k])/(grad[k][0]*dx + grad[k][1]*dy);
						}
						else if(qE[k] > maxQ[k])
						{
							s = (maxQ[k] - qM[k])/(grad[k][0]*dx + grad[k][1]*dy);
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
							//printf("Using LP limiter\n");
							auto A = Matrix!(10,2)(0);
							auto b = Vector!10(0);
							auto c = Vector!2(-abs(mesh.cells[i].gradient[k][0]), -abs(mesh.cells[i].gradient[k][1]));
							//auto xk = Vector!2(mesh.cells[i].lim[k][0], mesh.cells[i].lim[k][1]);
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

								A[cIdx+4,0] = (mesh.cells[cellIdx].centroid[0] - mesh.cells[i].centroid[0])*mesh.cells[i].gradient[k][0];
								A[cIdx+4,1] = (mesh.cells[cellIdx].centroid[1] - mesh.cells[i].centroid[1])*mesh.cells[i].gradient[k][1];
								b[cIdx+4] = mesh.cells[i].minQ[k] - q[i][k];

								A[cIdx+5,0] = -(mesh.cells[cellIdx].centroid[0] - mesh.cells[i].centroid[0])*mesh.cells[i].gradient[k][0];
								A[cIdx+5,1] = -(mesh.cells[cellIdx].centroid[1] - mesh.cells[i].centroid[1])*mesh.cells[i].gradient[k][1];
								b[cIdx+5] = q[i][k] - mesh.cells[i].maxQ[k];
							}

							LP!10(A, b, c, xk);

							//printf("xk = [%.10e, %.10e], lim = %.10e\n\n", xk[0], xk[1], mesh.cells[i].lim[k][0]);
							mesh.cells[i].lim[k][0] = xk[0];
							mesh.cells[i].lim[k][1] = xk[1];

							assert((xk[0] >= -10e-16) && (xk[0] <= (1.0 + 1e-12)));
							assert((xk[1] >= -10e-16) && (xk[1] <= (1.0 + 1e-12)));
						}
					}
				}
			
				for(uint j = 0; j < 4; j++)
				{
					mesh.cells[i].gradErr[j] = -abs(mesh.cells[i].gradient[j][0])*mesh.cells[i].lim[j][0] +
											   -abs(mesh.cells[i].gradient[j][1])*mesh.cells[i].lim[j][1] + 
											    abs(mesh.cells[i].gradient[j][0]) + abs(mesh.cells[i].gradient[j][1]);
				}
			}
		}
	}

	// update edge values and compute fluxes on boundary edges.
	foreach(i; mesh.boundaryEdges)
	{
		switch(mesh.edges[i].boundaryType)
			with(BoundaryType)
		{
			case FullState:
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
					//auto lim = mesh.cells[mesh.edges[i].cellIdx[0]].lim;
					
					for(uint j = 0; j < dims; j++)
					{
						mesh.edges[i].q[0][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[0]].lim[j][0]*grad[j][0]*dx + 
														mesh.cells[mesh.edges[i].cellIdx[0]].lim[j][1]*grad[j][1]*dy;
						//mesh.edges[i].q[0][j] = qM[j] + lim*(grad[j][0]*dx + grad[j][1]*dy);
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
							//mesh.edges[i].q[0][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[0]].lim[j]*(grad[j][0]*dx + grad[j][1]*dy);
							mesh.edges[i].q[0][j] = qM[j] + lim[0]*grad[j][0]*dx + lim[0]*grad[j][1]*dy;
						}
					}
				}

				auto qL = mesh.edges[i].q[0];
				auto qR = mesh.edges[i].q[1];

				//printf("edge = %d\n", i);
				mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax);
				if(mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN)
				{
					ex.SetException(SolverException.SExceptionType.EdgeException,
									"Got NaN on fullstate edge",
									SolverException.EdgeException(getPressure(mesh.edges[i].q[0]), 
																	getPressure(mesh.edges[i].q[1]),
																	mesh.edges[i].flux,
																	mesh.edges[i].q[0],
																	mesh.edges[i].q[1],
																	mesh.edges[i].normal,
																	mesh.edges[i].cellIdx[0],
																	mesh.edges[i].cellIdx[1]));
					throw ex;
				}
				break;
			case InviscidWall:
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
					//auto lim = mesh.cells[mesh.edges[i].cellIdx[0]].lim;
					
					for(uint j = 0; j < dims; j++)
					{
						mesh.edges[i].q[0][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[0]].lim[j][0]*grad[j][0]*dx + 
														mesh.cells[mesh.edges[i].cellIdx[0]].lim[j][1]*grad[j][1]*dy;
						//mesh.edges[i].q[0][j] = qM[j] + lim*(grad[j][0]*dx + grad[j][1]*dy);
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
							//mesh.edges[i].q[0][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[0]].lim[j]*(grad[j][0]*dx + grad[j][1]*dy);
							mesh.edges[i].q[0][j] = qM[j] + lim[0]*grad[j][0]*dx + lim[1]*grad[j][1]*dy;
						}
					}
				}

				Vector!2 velP = (1/mesh.edges[i].q[0][0])*Vector!2(mesh.edges[i].q[0][1], mesh.edges[i].q[0][2]);
				auto vel = (velP - (velP.dot(mesh.edges[i].normal))*mesh.edges[i].normal).magnitude;
				double p = (gamma - 1)*(mesh.edges[i].q[0][3] - 0.5*mesh.edges[i].q[0][0]*vel^^2);
				double a = sqrt(gamma*(p/mesh.edges[i].q[0][0]));
				if(p < 0)
				{
					p = 1.0e-12;
					//printf("pressure less than 0 at wall\n");
				}
				mesh.edges[i].flux = Vector!4(0, p*mesh.edges[i].normal[0], p*mesh.edges[i].normal[1], 0);
				mesh.edges[i].sMax = std.math.abs(a);
				
				auto qL = mesh.edges[i].q[0];
				auto qR = q[mesh.edges[i].cellIdx[1]];
				//mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax);
				
				if(mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN)
				{
					ex.SetException(SolverException.SExceptionType.EdgeException,
									"Got NaN on wall edge",
									SolverException.EdgeException(getPressure(mesh.edges[i].q[0]), 
																	getPressure(mesh.edges[i].q[1]),
																	mesh.edges[i].flux,
																	qL,
																	qR,
																	mesh.edges[i].normal,
																	mesh.edges[i].cellIdx[0],
																	mesh.edges[i].cellIdx[1]));
					throw ex;
				}
				break;
			default:
				ex.msg = "Unsupported boundary type";
				ex.file = __FILE__;
				ex.line = __LINE__;
				throw ex;
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
				//auto mid = mesh.edges[eIdx].mid;
				//Vector!2 mid = mesh.cells[mesh.edges[i].cellIdx[(k+1)%2]].centroid;

				auto dx = mid[0] - centroid[0];
				auto dy = mid[1] - centroid[1];
				//auto lim = mesh.cells[mesh.edges[i].cellIdx[k]].lim;
				
				for(uint j = 0; j < dims; j++)
				{
					mesh.edges[i].q[k][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[k]].lim[j][0]*grad[j][0]*dx + 
													mesh.cells[mesh.edges[i].cellIdx[k]].lim[j][1]*grad[j][1]*dy;
					//mesh.edges[i].q[k][j] = qM[j] + lim*(grad[j][0]*dx + grad[j][1]*dy);
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
						//mesh.edges[i].q[k][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[k]].lim[j]*(grad[j][0]*dx + grad[j][1]*dy);
						mesh.edges[i].q[k][j] = qM[j] + lim[0]*grad[j][0]*dx + lim[1]*grad[j][1]*dy;
					}
				}
			}
		}

		auto qL = mesh.edges[i].q[0];
		auto qR = mesh.edges[i].q[1];

		mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax);

		if(mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN)
		{
			ex.SetException(SolverException.SExceptionType.EdgeException,
							"Got NaN on interior edge",
							SolverException.EdgeException(getPressure(mesh.edges[i].q[0]), 
															getPressure(mesh.edges[i].q[1]),
															mesh.edges[i].flux,
															mesh.edges[i].q[0],
															mesh.edges[i].q[1],
															mesh.edges[i].normal,
															mesh.edges[i].cellIdx[0],
															mesh.edges[i].cellIdx[1]));
			throw ex;
		}
	}

	// need to do a half update of comm edges
	foreach(commIdx, commEdges; mesh.commEdgeIdx)
	{
		foreach(i; commEdges)
		{
			if(config.order == 1)
			{
				mesh.edges[i].q[0] = q[mesh.edges[i].cellIdx[0]];
				//mesh.edges[i].q[1] = q[mesh.edges[i].cellIdx[1]];
			}
			else
			{
				auto qM = q[mesh.edges[i].cellIdx[0]];
				auto grad = mesh.cells[mesh.edges[i].cellIdx[0]].gradient;
				auto centroid = mesh.cells[mesh.edges[i].cellIdx[0]].centroid;
				auto mid = mesh.edges[i].mid;
				//auto mid = mesh.edges[eIdx].mid;
				//Vector!2 mid = mesh.cells[mesh.edges[i].cellIdx[(0+1)%2]].centroid;

				auto dx = mid[0] - centroid[0];
				auto dy = mid[1] - centroid[1];
				//auto lim = mesh.cells[mesh.edges[i].cellIdx[0]].lim;
				
				for(uint j = 0; j < dims; j++)
				{
					mesh.edges[i].q[0][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[0]].lim[j][0]*grad[j][0]*dx + 
													mesh.cells[mesh.edges[i].cellIdx[0]].lim[j][1]*grad[j][1]*dy;
					//mesh.edges[i].q[0][j] = qM[j] + lim*(grad[j][0]*dx + grad[j][1]*dy);
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
						//mesh.edges[i].q[0][j] = qM[j] + mesh.cells[mesh.edges[i].cellIdx[0]].lim[j]*(grad[j][0]*dx + grad[j][1]*dy);
						mesh.edges[i].q[0][j] = qM[j] + lim[0]*grad[j][0]*dx + lim[1]*grad[j][1]*dy;
					}
				}
			}
		}
	}

	MPI_Barrier(mesh.comm);
	foreach(commIdx, commEdges; mesh.commEdgeIdx)
	{
		foreach(i, edge; commEdges)
		{
			mesh.stateBuffers[commIdx][i] = mesh.edges[edge].q[0];
		}
		MPI_Send(mesh.stateBuffers[commIdx].ptr, cast(uint)mesh.stateBuffers[commIdx].length, vec4dataType, mesh.commProc[commIdx], mesh.meshTag, mesh.comm);
	}

	for(uint c = 0; c < mesh.commProc.length; c++)
	{
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, mesh.meshTag, mesh.comm, &status);

		// determine which proc this is comming from, will be different order
		// than commProc
		if(mesh.commProc.canFind(status.MPI_SOURCE))
		{
			auto commIdx = mesh.commProc.countUntil(status.MPI_SOURCE);
			auto commEdges = mesh.commEdgeIdx[commIdx];

			MPI_Recv(mesh.stateBuffers[commIdx].ptr, cast(int)mesh.stateBuffers[commIdx].length, vec4dataType, mesh.commProc[commIdx], mesh.meshTag, mesh.comm, &status);

			foreach(i, edge; commEdges)
			{
				mesh.edges[edge].q[1] = mesh.stateBuffers[commIdx][i];
			}
		}
		else
		{
			printf("Unexpected source message from %d\n", status.MPI_SOURCE);
		}
	}

	foreach(commIdx, commEdges; mesh.commEdgeIdx)
	{
		foreach(i; commEdges)
		{
			auto qL = mesh.edges[i].q[0];
			auto qR = mesh.edges[i].q[1];

			mesh.edges[i].flux = F!dims(qL, qR, mesh.edges[i].normal, mesh.edges[i].sMax);

			if(mesh.edges[i].flux[0].isNaN || mesh.edges[i].flux[1].isNaN || mesh.edges[i].flux[2].isNaN || mesh.edges[i].flux[3].isNaN)
			{
				ex.SetException(SolverException.SExceptionType.EdgeException,
								"Got NaN on interior edge",
								SolverException.EdgeException(getPressure(mesh.edges[i].q[0]), 
																getPressure(mesh.edges[i].q[1]),
																mesh.edges[i].flux,
																mesh.edges[i].q[0],
																mesh.edges[i].q[1],
																mesh.edges[i].normal,
																mesh.edges[i].cellIdx[0],
																mesh.edges[i].cellIdx[1]));
				throw ex;
			}
		}
	}

	newDt = double.infinity;

	foreach(i; mesh.interiorCells)
	{
		R[i] = Vector!4(0);
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
					mesh.cells[i].dt = 0;
					//printf("R = [%.10e, %.10e, %.10e, %.10e]\n", R[0], R[1], R[2], R[3]);
					/+
					printf("%.10e\n", mesh.cells[i].perim);
					//printf("%.10e\n", mesh.cells[i].perim);
					printf("%.10e\n", sAve);
					ex.msg = "dt is NaN";
					ex.line = __LINE__;
					ex.file = __FILE__;
					throw ex;
					+/
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

	MPI_Allreduce(&newDt, &newDt, 1, MPI_DOUBLE, MPI_MIN, mesh.comm);
	MPI_Allreduce(&Rmax, &Rmax, 1, MPI_DOUBLE, MPI_MAX, mesh.comm);
}

void startComputation(Config config, string saveFile, uint p, uint id)
{
	try
	{
		auto umesh = UMesh2(MPI_COMM_SELF, id);

		double dt = config.dt;
		double t = 0;

		auto ex = new SolverException("No error");

		if(id == 0)
		{
			if(config.meshFile.canFind(".gri"))
			{
				umesh = parseXflowMesh(config.meshFile);
				umesh.comm = MPI_COMM_SELF;
				umesh.mpiRank = id;
			}
			else
			{
				writeln("Unsupported mesh format, exiting");
				return;
			}
		}

		if(p > 1)
		{
			umesh = partitionMesh(umesh, p, id, MPI_COMM_WORLD);
			umesh.comm = MPI_COMM_WORLD;
			umesh.mpiRank = id;
			//umesh.buildMesh;
			
			//core.stdc.stdlib.abort;
			//return;
		}

		umesh.buildMesh;

		import std.array : split;
		saveMatlabMesh(umesh, config.meshFile.split('.')[0]~"_"~id.to!string~".mmsh");

		final switch(config.limiter)
		{
			foreach(lim; limiterList)
			{
				case lim:
				{
					final switch(config.flux)
					{
						foreach(fl; fluxList)
						{
							case fl:
							{
								final switch(config.integrator)
								{
									foreach(inte; integratorList)
									{
										case inte:
											writeln("Running 2D finite volume solver");
											writeln("-limited: ", config.limited);
											writeln("-order: ", config.order);
											writeln("-limiter: "~lim);
											writeln("-flux: "~fl);
											writeln("-integrator: "~inte);
											runIntegrator!(ufvmSetup, ufvmSolver!(mixin(lim), mixin(fl), 4), mixin(inte))(umesh, config, saveFile, ex);
											break;
									}
								}
								break;
							}
						}
					}
					break;
				}
			}
		}
	}
	catch(SolverException ex)
	{
		writeln("Solver encountered an error: ", ex.msg);
		if(ex.exceptionType == SolverException.SExceptionType.EdgeException)
		{
			writeln("pL = ", ex.eExcept.pL);
			writeln("pR = ", ex.eExcept.pR);
			writeln("Flux = ", ex.eExcept.flux);
			writeln("qL = ", ex.eExcept.qL);
			writeln("qR = ", ex.eExcept.qR);
			writeln("cell L = ", ex.eExcept.cellL);
			writeln("cell R = ", ex.eExcept.cellR);
			writeln("normal = ", ex.eExcept.normal);
		}
		writeln("exiting");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
}

import mpi;
import mpi.util;

int main(string[] args)
{
	import std.getopt;

	string configFile;
	string saveFile = "";

	int argc = cast(int)args.length;
    auto argv = args.toArgv;

	int p;
	int id;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	vec4dataType = toMPIType!(Vector!4);

	double startTime = MPI_Wtime();

	signal(SIGINT, &handle);
	signal(SIGUSR1, &handle);
	auto res = getopt(args, "c|config", "config file to read", &configFile, 
							"s|save", "Save file to start from", &saveFile);

	auto configStr = readText(configFile);
	auto config = loadConfig(configStr);

	startComputation(config, saveFile, p, id);

	MPI_Barrier(MPI_COMM_WORLD);
	double elapsed = MPI_Wtime() - startTime;
	if(id == 0)
	{
		writeln("total time: ", elapsed);
	}
	
	return MPI_Shutdown;
}
