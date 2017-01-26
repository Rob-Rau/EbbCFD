module ebb.config;

import std.conv;
import std.json;
import std.stdio;

import ebb.mesh;
/+
{
	"mesh": "../box.mesh",
	"dt": 0.01,
	"tEnd": 750.0,
	"limiter": "minmodS",
	"flux": "rusanovFlux",
	"plotIter": 50,
	"saveIter": 50,
	"initialConditions": [0.85, 0, 1, 1.4], (M, aoa, p, rho)
	"boudaryConditions": [
		{
			"tag": "Bottom",
			"type": "const",
			"q": [0.85, 0, 1, 1.4]
		},
		{
			"tag": "Top",
			"type": "const",
			"q": [0.85, 0, 1, 1.4]
		},
		{
			"tag": "Left",
			"type": "const",
			"q": [0.85, 0, 1, 1.4]
		},
		{
			"tag": "Right",
			"type": "const",
			"q": [0.85, 0, 1, 1.4]
		},
		{
			"tag": "Airfoil",
			"type": "wall",
			"q": [0.85, 0, 1, 1.4]
		}
	]
}
+/

struct PhysicalConfig
{
	double Pr;
	double R;
	double gamma;
	double mu;
	double Re;
	double L;
}

struct Config
{
	string meshFile;
	double dt;
	bool dynamicDt;
	double tEnd;
	string limiter;
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
	string[] bTags;
	BoundaryType[] bTypes;
	double[][] bc;
	double aitkenTol = -1;
	bool localTimestep = false;
	bool multistageLimiting = false;
	bool cflAdjust = false;
	bool viscosity;
	PhysicalConfig physicalConfig;
}

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
	config.limiter = jConfig["limiter"].str;
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
		config.aitkenTol = jConfig["aitkenTol"].getDouble;
	}
	catch(Exception ex)
	{
		writeln("Aitken accelerator tolerance not supplied, disabling (-1)");
		config.aitkenTol = -1;
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
	for(uint i = 0; i < bcs.length; i++)
	{
		config.bTags ~= bcs[i]["tag"].str;
		immutable string bType = bcs[i]["type"].str;
		if(bType == "fullState")
		{
			config.bTypes ~= BoundaryType.FullState;
		}
		else if(bType == "inviscidWall")
		{
			config.bTypes ~= BoundaryType.InviscidWall;
		}
		else if(bType == "viscidWall")
		{
			config.bTypes ~= BoundaryType.ViscousWall;
		}
		else if(bType == "constP")
		{
			config.bTypes ~= BoundaryType.ConstPressure;
		}

		auto state = bcs[i]["q"].array;
		config.bc ~= [state[0].getDouble, state[1].getDouble, state[2].getDouble, state[3].getDouble];
	}

	return config;
}