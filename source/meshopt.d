module ebb.meshopt;

import core.atomic;
import core.stdc.stdio : fopen, fwrite, fopen, printf, snprintf;
import core.sys.posix.signal;

import std.algorithm : canFind, countUntil, min, max, reduce, sum;
import std.complex;
import std.conv;
import std.file;
import std.math;
import std.stdio;

import numd.utility;
import numd.linearalgebra.matrix;
import numd.optimization.ObjectiveFunction;
import numd.optimization.Gradient;
import numd.optimization.Derivative;
import numd.optimization.ComplexStep;
import numd.optimization.FiniteDifference;

import ebb.config;
import ebb.euler;
import ebb.exception;
import ebb.flux;
import ebb.integrators;
import ebb.limiters;
import ebb.mesh;
import ebb.io;
import ebb.mpid;

class MeshOpt : ObjectiveFunction
{
	this()
	{
		Constraints = 1;
	}

	final override Complex!double Compute(Complex!double[] designVar)
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

			@nogc double computeRMSResidual(double[] now, double[] last)
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

			residRho = computeRMSResidual(thisRho, lastRho);
			residU = computeRMSResidual(thisU, lastU);
			residV = computeRMSResidual(thisV, lastV);
			residE = computeRMSResidual(thisE, lastE);

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

		return complex(0.0, 0.0);
	}

	final override Complex!double[] Constraint(Complex!double[] designVar)
	{
		Complex!double[] c = new Complex!double[1];

		return c;
	}
}