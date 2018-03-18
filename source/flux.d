/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.flux;

import ebb.gas.flux;
import ebb.mesh;
import ebb.solve;

import std.algorithm;
import std.array;
import std.exception;
import std.math : abs, fmax, sqrt;
import std.meta : aliasSeqOf;
import std.range;
import std.typecons;

import numd.linearalgebra.matrix;

enum FluxMode
{
	Init,
	Add
}

interface IPhysics(SolverParams params)
{
	alias Vec = Vector!(cast(ulong)params.size);
	alias Normal = Vector!(cast(ulong)params.dims);
	alias Mat = Matrix!(cast(ulong)params.size, cast(ulong)params.dims);

	@nogc bool needEdgeGradients();
	@nogc void computeFluxes(Vec[] residuals, double[] waveSpeeds, Vec[] Ql, Vec[] Qr, Mat[] dQl, Mat[] dQr, immutable(Normal[]) normal);
	@nogc bool checkPhysics(ref Vec Q);
	@nogc void updateGhostCells(ref Vec[] Q, ref immutable Mesh mesh);
	@nogc void updateBoundaryEdgeFluxes(ref Vec[] residuals, ref double[] waveSpeeds, ref Vec[] qR, immutable(Vec[]) qL, immutable(Mat[]) dQl, immutable(Mat[]) dQr, immutable(Vec[]) Q, ref immutable Mesh mesh);
}

class Physics(SolverParams params, alias PhysicalConfig) : IPhysics!params
{
	immutable size = params.size;

	private FluxFunctionOverRange!(size, PhysicalConfig)[] rangeFluxFunctions;
	private FluxFunction!(size, PhysicalConfig)[] fluxFunctions;
	private PhysicalConfig config;

	alias ComputeFluxFunction = @nogc FluxResult!size delegate(Vec Ql, Vec Qr, Mat dQl, Mat dQr, immutable Normal normal);
	private ComputeFluxFunction computeFlux;

	this(string[] functions, PhysicalConfig config)
	{
		this.config = config;
		rangeFluxFunctions ~= getRangeFluxFunction!(size, FluxMode.Init, PhysicalConfig)(functions[0]);
		rangeFluxFunctions ~= functions[1..$].map!(a => a.getRangeFluxFunction!(size, FluxMode.Add, PhysicalConfig)).array;

		fluxFunctions = functions[0..$].map!(a => a.getFluxFunction!(size, PhysicalConfig)).array;

		computeFlux = (Vec Ql, Vec Qr, Mat dQl, Mat dQr, immutable Normal normal) {
			return fluxFunctions.fold!((res, func) {
					auto result = func(Ql, Qr, dQl, dQr, normal, config);
					res.flux += result.flux;
					res.waveSpeed = fmax(res.waveSpeed, result.waveSpeed);
					return res;
				})(fluxResult(Vector!size(0), -double.infinity));
		};
	}

	@nogc bool needEdgeGradients()
	{
		return config.needEdgeGradients;
	}

	@nogc void computeFluxes(Vec[] residuals, double[] waveSpeeds, Vec[] Ql, Vec[] Qr, Mat[] dQl, Mat[] dQr, immutable(Normal[]) normal)
	{
		rangeFluxFunctions.each!(func => func(residuals, waveSpeeds, Ql, Qr, dQl, dQr, normal, config));
	}

	@nogc bool checkPhysics(ref Vec Q)
	{
		return config.checkPhysics(Q);
	}

	@nogc void updateGhostCells(ref Vec[] Q, ref immutable Mesh mesh)
	{
		config.updateGhostCells(Q, mesh);
	}

	@nogc void updateBoundaryEdgeFluxes(ref Vec[] residuals, ref double[] waveSpeeds, ref Vec[] qR, immutable(Vec[]) qL, immutable(Mat[]) dQl, immutable(Mat[]) dQr, immutable(Vec[]) Q, ref immutable Mesh mesh)
	{
		foreach(i; mesh.boundaryEdges)
		{
			auto result = config.updateBoundaryEdgeFluxes!(size, ComputeFluxFunction)(qR[i], qL[i], dQl[i], dQr[i], Q, mesh.edges[i], mesh, computeFlux);
			residuals[i] = result.flux;
			waveSpeeds[i] = result.waveSpeed;
		}
	}
}

alias FluxFunction(uint size, alias PhysicalConfig) = @nogc FluxResult!size function(Vector!size Ql, Vector!size Qr, Matrix!(size,2) dQl, Matrix!(size,2) dQr, immutable Vector!(size - 2) normal, PhysicalConfig config);
alias FluxFunctionOverRange(uint size, alias PhysicalConfig) = @nogc void function(ref Vector!size[] R, ref double[] waveSpeeds, Vector!size[] Ql, Vector!size[] Qr, Matrix!(size,2)[] dQl, Matrix!(size,2)[] dQr, immutable(Vector!(size - 2)[]) normal, PhysicalConfig config);

@nogc void ComputeFluxes(uint size, FluxMode mode, alias PhysicalConfig, alias fluxFunction)(ref Vector!size[] R, ref double[] waveSpeeds, Vector!size[] Ql, Vector!size[] Qr, Matrix!(size,2)[] dQl, Matrix!(size,2)[] dQr, immutable(Vector!(size - 2)[]) normal, PhysicalConfig config)
{
	static if(mode == FluxMode.Init)
	{
		for(size_t i = 0; i < R.length; i++)
		{
			auto fluxResult = fluxFunction!size(Ql[i], Qr[i], dQl[i], dQr[i], normal[i], config);
			R[i] = fluxResult.flux;
			waveSpeeds[i] = fluxResult.waveSpeed;
		}
	}
	else
	{
		for(size_t i = 0; i < R.length; i++)
		{
			auto fluxResult = fluxFunction!size(Ql[i], Qr[i], dQl[i], dQr[i], normal[i], config);
			R[i] += fluxResult.flux;
			waveSpeeds[i] = fmax(fluxResult.waveSpeed, waveSpeeds[i]);
		}
	}
}

@nogc FluxResult!size ComputeFlux(uint size, alias PhysicalConfig, alias fluxFunction)(Vector!size Ql, Vector!size Qr, Matrix!(size,size-2) dQl, Matrix!(size,size-2) dQr, immutable Vector!(size - 2) normal, PhysicalConfig config)
{
	return fluxFunction!size(Ql, Qr, dQl, dQr, normal, config);
}

FluxFunctionOverRange!(size, PhysicalConfig) getRangeFluxFunction(uint size, FluxMode mode, alias PhysicalConfig)(string name)
{
	if(name == "roeFlux")
	{
		return &ComputeFluxes!(size, mode, PhysicalConfig, roeFlux);
	}
	else if(name == "rusanovFlux")
	{
		return &ComputeFluxes!(size, mode, PhysicalConfig, rusanovFlux);
	}
	else if(name == "averagedDiffusiveFlux")
	{
		return &ComputeFluxes!(size, mode, PhysicalConfig, averagedDiffusiveFlux);
	}
	throw new Exception("whoops");
}

FluxFunction!(size, PhysicalConfig) getFluxFunction(uint size, alias PhysicalConfig)(string name)
{
	if(name == "roeFlux")
	{
		return &ComputeFlux!(size, PhysicalConfig, roeFlux);
	}
	else if(name == "rusanovFlux")
	{
		return &ComputeFlux!(size, PhysicalConfig, rusanovFlux);
	}
	else if(name == "averagedDiffusiveFlux")
	{
		return &ComputeFlux!(size, PhysicalConfig, averagedDiffusiveFlux);
	}
	throw new Exception("whoops");
}
