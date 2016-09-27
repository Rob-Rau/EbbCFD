/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.solver;

import std.algorithm : min, max, reduce;
import std.conv;
import std.json;
import std.file;
import std.math;
import std.stdio;

import rpp.client.rpc;

import numd.utility;
import numd.linearalgebra.matrix;

import ebb.limiters;
import ebb.euler;
import ebb.flux;
import ebb.mesh;

double[] abs(double[] data)
{
	import std.math : abs;
	for(int i = 0; i < data.length; i++)
	{
		data[i] = abs(data[i]); 
	}
	return data;
}

alias solverList = aliasSeqOf!(["ufvmSolver", "sfvmSolver"]);

// Unstructured finite volume solver
void ufvmSolver(alias S, alias F, size_t dims)(ref UMesh2 mesh, Config config, double t)
{
	double dt = config.dt;
	uint iterations = 0;

	// Setup initial conditions
	for(uint i = 0; i < mesh.cells.length; i++)
	{
		double M = config.ic[0];
		double aoa = config.ic[1] * (PI/180);
		double p = config.ic[2];
		double rho = config.ic[3];
		double a = sqrt(gamma*(p/rho));
		double U = M*a;
		double u = U*cos(aoa);
		double v = U*sin(aoa);

		mesh.cells[i].q = buildQ(rho, u, v, p);
	}

	for(uint i = 0; i < mesh.bGroups.length; i++)
	{
		uint findBcIndex(string tag)
		{
			for(uint j = 0; j < config.bTags.length; j++)
			{
				if(config.bTags[j] == tag)
				{
					return j;
				}
			}
			assert(false, "Could not find match boundary condition tag");
		}

		uint bcIdx = findBcIndex(mesh.bTags[i]);

		for(uint j = 0; j < mesh.bGroups[i].length; j++)
		{
			assert(mesh.edges[mesh.bGroups[i][j]].isBoundary, "Edge not boundary edge but should be");
			assert(mesh.edges[mesh.bGroups[i][j]].boundaryTag == config.bTags[bcIdx], "Incorrect boundary tag");

			double M = config.bc[bcIdx][0];
			double aoa = config.bc[bcIdx][1] * (PI/180);
			double p = config.bc[bcIdx][2];
			double rho = config.bc[bcIdx][3];
			double a = sqrt(gamma*(p/rho));
			double U = M*a;
			double u = U*cos(aoa);
			double v = U*sin(aoa);
			mesh.edges[mesh.bGroups[i][j]].boundaryType = config.bTypes[bcIdx];

			if(mesh.edges[mesh.bGroups[i][j]].boundaryType == BoundaryType.FullState)
			{
				mesh.edges[mesh.bGroups[i][j]].q[1] = buildQ(rho, u, v, p);
			}
		}
	}

	while(!approxEqual(t, config.tEnd))
	{
		//writeln;
		//writeln;
		for(uint i = 0; i < mesh.edges.length; i++)
		{
			if(mesh.edges[i].isBoundary)
			{
				switch(mesh.edges[i].boundaryType)
					with(BoundaryType)
				{
					case FullState:
						mesh.edges[i].q[0] = mesh.cells[mesh.edges[i].cellIdx[0]].q;

						auto qL = mesh.edges[i].q[0];
						auto qR = mesh.edges[i].q[1];

						Vector!2 velL = mesh.edges[i].rotMat*Vector!2(qL[1], qL[2]);
						Vector!2 velR = mesh.edges[i].rotMat*Vector!2(qR[1], qR[2]);

						qL[1] = velL[0];
						qL[2] = velL[1];

						qR[1] = velR[0];
						qR[2] = velR[1];

						mesh.edges[i].flux = F!(dims, 1, 0)(qR, qL, dt, 0);
						mesh.edges[i].flux.writeln;
						break;
					case InviscidWall:
						mesh.edges[i].q[0] = mesh.cells[mesh.edges[i].cellIdx[0]].q;
						Vector!2 velL = (1/mesh.edges[i].q[0][0])*mesh.edges[i].rotMat*Vector!2(mesh.edges[i].q[0][1], mesh.edges[i].q[0][2]);
						double p = (gamma - 1)*(mesh.edges[i].q[0][3] - 0.5*mesh.edges[i].q[0][0]*(velL[1]^^2.0));
						mesh.edges[i].flux = Vector!4(0, p, 0, 0);
						//writeln("p = ", p);
						break;
					default:
						assert(false, "Unsupported boundary type");
				}
			}
			else
			{
				mesh.edges[i].q[0] = mesh.cells[mesh.edges[i].cellIdx[0]].q;
				mesh.edges[i].q[1] = mesh.cells[mesh.edges[i].cellIdx[1]].q;

				auto qL = mesh.edges[i].q[0];
				auto qR = mesh.edges[i].q[1];

				Vector!2 velL = mesh.edges[i].rotMat*Vector!2(qL[1], qL[2]);
				Vector!2 velR = mesh.edges[i].rotMat*Vector!2(qR[1], qR[2]);

				qL[1] = velL[0];
				qL[2] = velL[1];

				qR[1] = velR[0];
				qR[2] = velR[1];

				mesh.edges[i].flux = F!(dims, 1, 0)(qR, qL, dt, 0);
				//mesh.edges[i].flux.writeln;
			}
		}

		for(uint i = 0; i < mesh.cells.length; i++)
		{
			auto R = Vector!4(0);
			// integrate fluxes over cell edges
			for(uint j = 0; j < mesh.cells[i].nEdges; j++)
			{
				R += mesh.cells[i].fluxMultiplier[j]*mesh.edges[mesh.cells[i].edges[j]].len*mesh.edges[mesh.cells[i].edges[j]].flux;
			}
			//R.writeln;
			//mesh.cells[i].area.writeln;
			mesh.cells[i].q = mesh.cells[i].q - (dt/mesh.cells[i].area)*R;
			//mesh.cells[i].q.writeln;
		}

		auto rotMat = Matrix!(2, 2)(cos(config.ic[1] * (PI/180)), -sin(config.ic[1] * (PI/180)), sin(config.ic[1] * (PI/180)), cos(config.ic[1] * (PI/180)));
		auto f = mesh.computeBoundaryForces("Airfoil");
		auto ld = rotMat*f;

		if(iterations % config.plotIter == 0)
		{
			import core.stdc.stdio;
			printf("lift force = %f\t t = %f\n", ld[0], t);
		}
		t += dt;
		iterations++;
	}
}

// Structured finite volume solver
void sfvmSolver(alias S, alias F, size_t dims)(ref Mesh mesh, Config config, double t)
{
	alias Vec = Vector!dims;
	alias Mat = Matrix!(dims, dims);
	
	double dt = config.dt;
	immutable double tEnd = config.tEnd;
	
	int iterations = 0;//round(tEnd/dt).to!int;
	int saveIterations = 0;
	
	double newdt = ((mesh[0,0].dx*mesh[0,0].dy)/(mesh.getSpeed2.reduce!max + mesh.getSoundSpeed.reduce!max));
	writeln(newdt);

	while(!approxEqual(t, tEnd))// && (t < tEnd))
	{
		mesh.updateGhosts();
		mesh.updateGhosts();
		//dt = 0.5*((mesh[0,0].dx*mesh[0,0].dy)/(mesh.getSpeed2.reduce!max + mesh.getSoundSpeed.reduce!max));
		//writeln(newdt);
		for(int i = 0; i < mesh.N; i++)
		{
			for(int j = 0; j < mesh.M; j++)
			{
				if(mesh[i, j].cellType == CellType.Normal)
				//if(mesh[i, j].cellType != CellType.Solid)
				{
					Vec xSlopes;
					Vec ySlopes;
					for(int n = 0; n < dims; n++)
					{
						xSlopes[n] = S(mesh[i, j].xSp[n], mesh[i, j].xSm[n]);
						ySlopes[n] = S(mesh[i, j].ySp[n], mesh[i, j].ySm[n]);
					}

					// half timestep update in x dir
					auto Lmat = L!(1, 0)(mesh[i, j].q);
					auto Rmat = R!(1, 0)(mesh[i, j].q);
					auto LamMat = Lam!(1, 0)(mesh[i, j].q);
					Vec qHalfx = mesh[i, j].q - 0.5*(dt/mesh[i, j].dx)*(Rmat*LamMat*Lmat)*xSlopes;

					//qHalfx.ToString.writeln;
					mesh[i, j].qL = qHalfx + (-0.5*mesh[i, j].dx)*xSlopes;
					mesh[i, j].qR = qHalfx + (0.5*mesh[i, j].dx)*xSlopes;

					if(mesh[i, j].qL[0] < 0)
					{
						writeln("whoops, we got a negative density on left edge "~i.to!string~" "~j.to!string);
						mesh[i, j].qL[0] = 0.0001;
					}
					if(mesh[i, j].qR[0] < 0)
					{
						writeln("whoops, we got a negative density on right edge "~i.to!string~" "~j.to!string);
						mesh[i, j].qR[0] = 0.0001;
					}
					// half timestep update in y dir
					Lmat = L!(0, 1)(mesh[i, j].q);
					Rmat = R!(0, 1)(mesh[i, j].q);
					LamMat = Lam!(0, 1)(mesh[i, j].q);
					Vec qHalfy = mesh[i, j].q - 0.5*(dt/mesh[i, j].dy)*(Rmat*LamMat*Lmat)*ySlopes;

					mesh[i, j].qB = qHalfy + (-0.5*mesh[i, j].dy)*ySlopes;
					mesh[i, j].qT = qHalfy + (0.5*mesh[i, j].dy)*ySlopes;
					if(mesh[i, j].qB[0] < 0)
					{
						writeln("whoops, we got a negative density on bottom edge "~i.to!string~" "~j.to!string);
						mesh[i, j].qB[0] = 0.0001;
					}
					if(mesh[i, j].qT[0] < 0)
					{
						writeln("whoops, we got a negative density on top edge "~i.to!string~" "~j.to!string);
						mesh[i, j].qT[0] = 0.0001;
					}
					
					if(mesh[i-1, j].cellType == CellType.GhostMirrorYB)
					{
						mesh[i, j].xFlux = F!(dims, 1, 0)(mesh[i-1,j].q, mesh[i, j].q, dt, mesh[i, j].dx);
						mesh[i, j].yFlux = F!(dims, 0, 1)(mesh[i,j-1].qT, mesh[i, j].qB, dt, mesh[i, j].dy);
					}
					else if((mesh[i-1, j].cellType == CellType.GhostMirrorYT) && !mesh[i-1, j].corner)
					{
						mesh[i, j].xFlux = F!(dims, 1, 0)(mesh[i-1,j].q, mesh[i, j].qL, dt, mesh[i, j].dx);
						//mesh[i, j].qB = mesh[i, j].q;
						mesh[i, j-1].qT = mesh[i, j].qB;
						mesh[i, j-1].qT[2] = -mesh[i, j].qB[2];
						mesh[i, j].yFlux = F!(dims, 0, 1)(mesh[i,j-1].qT, mesh[i, j].qB, dt, mesh[i, j].dy);
					}
					else if((mesh[i, j-1].cellType == CellType.GhostMirrorYT) && mesh[i, j-1].corner) 
					{
						mesh[i, j].xFlux = F!(dims, 1, 0)(mesh[i-1,j].qR, mesh[i, j].qL, dt, mesh[i, j].dx);
						mesh[i, j].yFlux = F!(dims, 0, 1)(mesh[i,j-1].q, mesh[i, j].qB, dt, mesh[i, j].dy);
					}
					else
					{
						if(mesh[i,j-1].cellType == CellType.GhostNoGradYT)
						{
							mesh[i,j-1].qT = mesh[i, j].qB;
						}
						// update fluxes for this cell.
						mesh[i, j].xFlux = F!(dims, 1, 0)(mesh[i-1,j].qR, mesh[i, j].qL, dt, mesh[i, j].dx);
						mesh[i, j].yFlux = F!(dims, 0, 1)(mesh[i,j-1].qT, mesh[i, j].qB, dt, mesh[i, j].dy);
					}
				}
				else if((mesh[i, j].cellType == CellType.GhostMirrorXL) || (mesh[i, j].cellType == CellType.GhostNoGradXL) || (mesh[i, j].cellType == CellType.GhostConstPressureXL))
				{
					mesh[i, j].xFlux = F!(dims, 1, 0)(mesh[i-1,j].qR, mesh[i, j].qL, dt, mesh[i, j].dx);
				}
				else if((mesh[i, j].cellType == CellType.GhostMirrorYB) && (mesh[i, j].corner))
				{
					mesh[i, j].yFlux = F!(dims, 0, 1)(mesh[i,j-1].qT, mesh[i, j].qB, dt, mesh[i, j].dy);
					if(mesh[i,j].cornerType == CellType.GhostMirrorXL)
					{
						mesh[i, j].xFlux = F!(dims, 1, 0)(mesh[i-1,j].qR, mesh[i, j].qL, dt, mesh[i, j].dx);
					}
				}
				else if(((mesh[i, j].cellType == CellType.GhostMirrorYB) || (mesh[i, j].cellType == CellType.GhostNoGradYB) ||
						 (mesh[i, j].cellType == CellType.GhostConstPressureYB)) && !mesh[i, j].corner)
				{
					if(mesh[i, j].cellType == CellType.GhostMirrorYB)
					{
						mesh[i, j].qB = mesh[i,j-1].qT;
						mesh[i, j].qB[2] = -mesh[i,j-1].qT[2];
					}
					else if(mesh[i, j].cellType == CellType.GhostNoGradYB)
					{
						mesh[i, j].qB = mesh[i,j-1].qT;
					}
					
					mesh[i, j].yFlux = F!(dims, 0, 1)(mesh[i,j-1].qT, mesh[i, j].qB, dt, mesh[i, j].dy);
				}
				else if(((mesh[i, j].cellType == CellType.GhostMirrorYT) || (mesh[i, j].cellType == CellType.GhostNoGradYT)) && (mesh[i, j].cornerType == CellType.GhostMirrorXL))
				{
					mesh[i, j].xFlux = F!(dims, 1, 0)(mesh[i-1,j].qR, mesh[i, j].qL, dt, mesh[i, j].dx);
				}
				else if((mesh[i, j].cellType == CellType.GhostConst) && (mesh[i,j-1].cellType == CellType.Normal))
				{
					mesh[i, j].yFlux = F!(dims, 0, 1)(mesh[i,j-1].qT, mesh[i, j].qB, dt, mesh[i, j].dy);
				}
			}
		}

		for(int i = 0; i < mesh.N; i++)
		{
			for(int j = 0; j < mesh.M; j++)
			{
				if(mesh[i, j].cellType == CellType.Normal)
				{
					mesh[i, j].q = mesh[i, j].q - (dt/mesh[i, j].dx)*(mesh[i+1,j].xFlux - mesh[i, j].xFlux) - (dt/mesh[i, j].dy)*(mesh[i,j+1].yFlux - mesh[i, j].yFlux);

					if(mesh[i, j].q[0] < 0)
					{
						writeln("whoops, we got a negative density");
						mesh[i, j].q[0] = 0;
					}
				}
			}
		}

		t += dt;
		//iterations.writeln;
		iterations++;
		
		import std.format : format;
		import std.string : rightJustify;
		if(config.plotIter > 0)
		{
			if(iterations % config.plotIter == 0)
			{
				auto meshgrid = buildMeshgrid(mesh);
				
				//auto p = getVelocity!0(mesh);
				//auto p = getVelocity!1(mesh);
				//auto p = getDensity(mesh);
				auto p = getPressure(mesh);
				//auto p = getMach(mesh);
				contourf(meshgrid.X, meshgrid.Y, p, 50, `LineStyle`, `none`);

				//caxis([0.85, 1.3]);
				//caxis([0.0, 3.0]);
				//caxis([0.9, 1.15]);
				//caxis([1.0, 1.25]);
				colorbar;
				axis("equal");
				hold!"on";
				title(format("t = %4.8f", t));
			}
		}

		if(config.saveIter > 0)
		{
			if((iterations % config.saveIter) == 0)
			{
				saveMesh(mesh, format("save_%s.mesh", saveIterations.to!string.rightJustify(7, '0')), dt, t);
				saveIterations++;
			}
		}
	}
}

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
struct Config
{
	string meshFile;
	double dt;
	double tEnd;
	string limiter;
	string flux;
	long saveIter;
	long plotIter;
	string solver;
	double[4] ic;
	string[] bTags;
	BoundaryType[] bTypes;
	double[][] bc;
}

Config loadConfig(string conf, string saveFile)
{
	JSONValue jConfig = parseJSON(conf);
	Config config;

	if(saveFile == "")
	{
		config.meshFile = jConfig["mesh"].str;
	}
	else
	{
		config.meshFile = saveFile;
	}

	config.limiter = jConfig["limiter"].str;
	config.flux = jConfig["flux"].str;
	config.dt = jConfig["dt"].floating;
	config.tEnd = jConfig["tEnd"].floating;
	config.saveIter = jConfig["saveIter"].integer;
	config.plotIter = jConfig["plotIter"].integer;
	config.solver = jConfig["solver"].str;

	if(config.solver == "ufvmSolver")
	{
		double getDouble(T)(T val)
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

		auto ics = jConfig["initialConditions"].array;
		config.ic[0] = getDouble(ics[0]);
		config.ic[1] = getDouble(ics[1]);
		config.ic[2] = getDouble(ics[2]);
		config.ic[3] = getDouble(ics[3]);

		auto bcs = jConfig["boudaryConditions"].array;
		for(uint i = 0; i < bcs.length; i++)
		{
			config.bTags ~= bcs[i]["tag"].str;
			string bType = bcs[i]["type"].str;
			if(bType == "fullState")
			{
				config.bTypes ~= BoundaryType.FullState;
			}
			else if(bType == "inviscidWall")
			{
				config.bTypes ~= BoundaryType.InviscidWall;
			}
			else if(bType == "constP")
			{
				config.bTypes ~= BoundaryType.ConstPressure;
			}

			auto state = bcs[i]["q"].array;
			config.bc ~= [getDouble(state[0]), getDouble(state[1]), getDouble(state[2]), getDouble(state[3])];
		}
	}
	return config;
}

void startComputation(Config config)
{
	Mesh mesh;
	UMesh2 umesh;

	double dt = config.dt;
	double tEnd = config.tEnd;
	double t = 0;

	if(config.solver == "sfvmSolver")
	{
		mesh = loadMesh(config.meshFile, dt, t);
		mesh.updateGhosts();
	}
	else if(config.solver == "ufvmSolver")
	{
		umesh = parseXflowMesh(config.meshFile);
	}

	switch(config.limiter)
	{
		foreach(lim; limiterList)
		{
			case lim:
				switch(config.flux)
				{
					foreach(fl; fluxList)
					{
						case fl:
							writeln("Running 2D finite volume solver");
							writeln("-limiter: "~lim);
							writeln("-flux: "~fl);
							if(config.solver == "sfvmSolver")
							{
								writeln("-solver: sfvmSolver");
								sfvmSolver!(mixin(lim), mixin(fl), 4)(mesh, config, t);
							}
							else if(config.solver == "ufvmSolver")
							{
								writeln("-solver: ufvmSolver");
								ufvmSolver!(mixin(lim), mixin(fl), 4)(umesh, config, t);
							}
							break;
					}
					default:
						writeln("Invalid flux function");
						break;
				}
				break;
		}
		default:
			writeln("Invalid limiter function");
			break;
	}
}

void main(string[] args)
{
	import std.getopt;
	string configFile;
	string plotAddr = "127.0.0.1";
	string saveFile = "";
	ushort plotPort = 54000;

	auto res = getopt(args, "c|config", "config file to read", &configFile, "pa|plotAddr", "IP address to plot to", &plotAddr, 
							"pp|plotPort", "Port to plot to", &plotPort, "s|save", "Save file to start from", &saveFile);

	initRPP(plotAddr, plotPort);

	auto configStr = readText(configFile);
	auto config = loadConfig(configStr, saveFile);

	startComputation(config);
}
