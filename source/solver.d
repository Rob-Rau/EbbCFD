/+ Copyright (c) 2016 Robert F. Rau II +/
import std.json;
import std.file;
import std.math : fmax, fmin, sgn, sqrt, exp, log;
import std.stdio;

import rpp.client.rpc;

import numd.utility;
import numd.linearalgebra.matrix;

import limiters;
import euler;
import flux;
import mesh;

double[] abs(double[] data)
{
	import std.math : abs;
	for(int i = 0; i < data.length; i++)
	{
		data[i] = abs(data[i]); 
	}
	return data;
}

void finiteVolumeSolver(alias S, alias F, size_t dims)(ref Mesh mesh, Config config, double t)
{
	import std.algorithm : min, max, reduce;
	import std.conv : to;
	import std.math : approxEqual, round;
		
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
				else if((mesh[i, j].cellType == CellType.GhostMirrorXL) || (mesh[i, j].cellType == CellType.GhostNoGradXL))
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
				else if(((mesh[i, j].cellType == CellType.GhostMirrorYB) || (mesh[i, j].cellType == CellType.GhostNoGradYB)) && !mesh[i, j].corner)
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
	"flux": "rusanovFlux"
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
	
	return config;
}

void startComputation(Config config)
{
	Mesh mesh;
	double dt = config.dt;
	double tEnd = config.tEnd;
	double t = 0;

	mesh = loadMesh(config.meshFile, dt, t);

	mesh.updateGhosts();

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
							finiteVolumeSolver!(mixin(lim), mixin(fl), 4)(mesh, config, t);
							break;
					}
					default:
						break;
				}
				break;
		}
		default:
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
