/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.config;

import std.conv;
import std.json;
import std.stdio;

import numd.linearalgebra.matrix;

import ebb.gas.config;
import ebb.mesh;

/+
{
	"mesh": "mesh.gmsh",
	"plotIter": -1,
	"saveIter": -1,
	"integrator": "Euler",
	"integratorConfig": {
		"timestep": 0.001,
		"localTimestep": false,
		"dynamicDt": false,
		"steadyConvergenceConfig": {
			"residuals": []
		}
	},
	"spacialSolver": "finiteVolume",
	"finiteVolumeConfig": {
		"order": 2,
		"limit": true,
		"lpThresh": -1,
		"forceBoundaries": ["Wall1", "Wall2"],
		"CFL": 0.3,
		"adjustCFL": false,
		"flux": ["roeFlux", "averagedDiffusiveFlux"],
		"gasPhysicalConfig": {
			"Pr": 0.7,
			"gamma": 1.4,
			"mu": 0.0001
		}
	},
	"initialConditions": [
		{
			"tag": "tag1",
			"state": [0.1, 0, 1.1, 1] //[rho, rhoU, rhoV, E]
		},
		{
			"tag": "tag2",
			"state": [0.1, 0, 1.1, 1] //[rho, rhoU, rhoV, E]
		}
	],
	"boudaryConditions": [
		{
			"tag": "Freestream",
			"type": "fullState",
			"state": [0.1, 0, 1.1, 1] // [rho, rhoU, rhoV, E]
		},
		{
			"tag": "Wall",
			"type": "viscidWall",
			"state": [0.1, 0, 1.1, 1] // [rho, rhoU, rhoV, E]
		},
		{
			"tag": "Outflow",
			"type": "constP",
			"state": [71.429] // p
		},
		{
			"tag": "Symmetric",
			"type": "symmetry",
			"state": [71.429] // p
		}
	]
}
+/

alias @nogc Vector!4 function(double x, double y) DirichletFunc;
alias @nogc Matrix!(4, 2) function(double x, double y) DirichletFuncDerivative;

struct BoundaryData
{
	double[] boundaryData;
	string tag;
	BoundaryType type;
	DirichletFunc dFunc;
	DirichletFuncDerivative dFuncDerivative;
}


struct Config
{
	string meshFile;
	double dt;
	bool dynamicDt;
	double tEnd;
	string flux;
	long saveIter;
	long plotIter;
	string[] forceBoundary;
	string integrator;
	bool limited;
	double lpThresh;
	long order;
	double CFL;
	double[4] ic;
	BoundaryData[] boundaries;
	bool localTimestep = false;
	bool multistageLimiting = false;
	bool cflAdjust = false;
	bool viscosity;
	//PhysicalConfig physicalConfig;
}
/+
static double getDouble(T)(T val)
{
	if(val.type == JSON_TYPE.INTEGER)
	{
		return val.integer.to!double;
	}
	else if(val.type == JSON_TYPE.FLOAT)
	{
		return val.floating;
	}
	else
	{
		assert(false, "invalid type");
	}
}

Config loadConfig(string conf)
{
	JSONValue jConfig = parseJSON(conf);
	Config config;

	config.meshFile = jConfig["mesh"].str;
	config.flux = jConfig["flux"].str;
	config.dt = jConfig["dt"].floating;
	config.tEnd = jConfig["tEnd"].getDouble;
	config.saveIter = jConfig["saveIter"].integer;
	config.plotIter = jConfig["plotIter"].integer;

	try
	{
		auto forceBoundaries = jConfig["forceBoundary"].array;
		for(uint i = 0; i < forceBoundaries.length; i++)
		{
			config.forceBoundary ~= forceBoundaries[i].str;
		}
	}
	catch(Exception ex)
	{
		config.forceBoundary ~= jConfig["forceBoundary"].str;
	}

	try
	{
		config.CFL = jConfig["CFL"].getDouble;
	}
	catch(Exception ex)
	{
		writeln("CFL not provided, setting to 0.5");
		config.CFL = 0.5;
	}

	try
	{
		config.dynamicDt = (jConfig["dynamicDt"].type == JSON_TYPE.TRUE);
	}
	catch(Exception ex)
	{
		writeln("dynamicDt not provided, enabling");
		config.dynamicDt = true;
	}

	try
	{
		config.integrator = jConfig["integrator"].str;
	}
	catch(Exception ex)
	{
		writeln("Integrator not provided, using forward euler");
		config.integrator = "Euler";
	}

	try
	{
		config.limited = (jConfig["limited"].type == JSON_TYPE.TRUE);
	}
	catch(Exception ex)
	{
		writeln("Limited option not provided, setting to true");
		config.limited = true;
	}

	try
	{
		config.localTimestep = (jConfig["localTimestep"].type == JSON_TYPE.TRUE);
	}
	catch(Exception ex)
	{
		writeln("Local timestepping option not provided, setting to false");
		config.localTimestep = false;
	}

	try
	{
		config.lpThresh = jConfig["lpThresh"].getDouble;
	}
	catch(Exception ex)
	{
		writeln("LP threshold not provided, turning off (-1)");
		config.lpThresh = -1;
	}

	try
	{
		config.order = jConfig["order"].integer;

		if((config.order != 1) && (config.order != 2))
		{
			writeln("Invalid order supplied, setting order 2");
			config.order = 2;
		}
	}
	catch(Exception ex)
	{
		writeln("Order option not provided, setting order 2");
		config.order = 2;
	}

	try
	{
		config.multistageLimiting = (jConfig["multistageLimiting"].type == JSON_TYPE.TRUE);
	}
	catch(Exception ex)
	{
		writeln("Multistage Limiting option not provided, enabling");
		config.multistageLimiting = true;
	}

	try
	{
		config.cflAdjust = (jConfig["cflAdjust"].type == JSON_TYPE.TRUE);
	}
	catch(Exception ex)
	{
		writeln("CFL adjustment option not provided, enabling");
		config.cflAdjust = true;
	}

	try
	{
		config.viscosity = (jConfig["viscosity"].type == JSON_TYPE.TRUE);
	}
	catch(Exception ex)
	{
		writeln("viscosity not provided, disabling");
		config.viscosity = false;
	}

	try
	{
		auto physConf = jConfig["physicalConfig"];
		config.physicalConfig.Pr = physConf["Pr"].getDouble;
		config.physicalConfig.R = physConf["R"].getDouble;
		config.physicalConfig.gamma = physConf["gamma"].getDouble;
		config.physicalConfig.Re = physConf["Re"].getDouble;
		config.physicalConfig.L = physConf["L"].getDouble;
	}
	catch(Exception ex)
	{
		writeln("physicalConfig not provided, setting to defaults");
		config.physicalConfig.Pr = 0.7;
		config.physicalConfig.R = 8.3145;
		config.physicalConfig.gamma = 1.4;
		config.physicalConfig.Re = 1000;
		config.physicalConfig.L = 1;
	}

	auto ics = jConfig["initialConditions"].array;
	config.ic[0] = ics[0].getDouble;
	config.ic[1] = ics[1].getDouble;
	config.ic[2] = ics[2].getDouble;
	config.ic[3] = ics[3].getDouble;

	auto bcs = jConfig["boudaryConditions"].array;
	config.boundaries.length = bcs.length;
	for(uint i = 0; i < bcs.length; i++)
	{
		config.boundaries[i].bTag = bcs[i]["tag"].str;
		immutable string bType = bcs[i]["type"].str;
		if(bType == "fullState")
		{
			config.boundaries[i].type = BoundaryType.FullState;
		}
		else if(bType == "inviscidWall")
		{
			config.boundaries[i].type = BoundaryType.InviscidWall;
		}
		else if(bType == "viscidWall")
		{
			config.boundaries[i].type = BoundaryType.ViscousWall;
		}
		else if(bType == "constP")
		{
			config.boundaries[i].type = BoundaryType.ConstPressure;
		}
		else if(bType == "symmetry")
		{
			config.boundaries[i].type = BoundaryType.Symmetry;
		}
		else if(bType == "tpInflow")
		{
			config.boundaries[i].type = BoundaryType.TempPresInflow;
		}
		else if(bType == "dirichlet")
		{
			config.boundaries[i].type = BoundaryType.Dirichlet;
		}

		auto state = bcs[i]["q"].array;
		//config.bc ~= [state[0].getDouble, state[1].getDouble, state[2].getDouble, state[3].getDouble];
		//config.boundaries[i].boundaryData = state;
		foreach(bc; bcs[i]["q"].array)
		{
			config.boundaries[i].boundaryData ~= bc.getDouble;
		}
	}

	return config;
}
+/
