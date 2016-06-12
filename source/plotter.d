module plotter;

import std.getopt;
import std.stdio;

import rpp.client.rpc;

import numd.utility;

import mesh;
import euler;

static this()
{
	initRPP("127.0.0.1", 54000);
}

void main(string[] args)
{
	string file, outputFile;
	bool plotU = false, plotV = false, plotP = false, plotRho = false, plotE = false, plotM = false, plotSpeed = false, movie = false;
	
	auto res = getopt(args, std.getopt.config.caseSensitive, "file|f", "file to plot", &file, "u|u", "plot u velocity component", &plotU, "v|v", "plot v velocity component", &plotV,
							"p|p", "plot pressure", &plotP, "d|d", "plot density", &plotRho, "e|e", "plot energy", &plotE,
							"M|M", "plot Mach number", &plotM, "s|s", "plot flow speed", &plotSpeed,
							"m|m", "generate output images for movie", &movie, "o|o", "output file name", &outputFile);
	if(res.helpWanted)
	{
		writeln("plotter options:");
		foreach(opt; res.options)
		{
			writeln(opt.optShort, " | ", opt.optLong, "\t\t", opt.help);
		}
		return;
	}

	import std.file : dirEntries, mkdir, SpanMode;
	import std.array : array;
	import std.algorithm : sort;
	
	if(file == "")
	{
		writeln("no input file, exiting");
		return;
	}
	
	movie.writeln;
	plotM.writeln;
	Mesh mesh;
	double dt, t;
	mesh = loadMesh(file, dt, t);

	Meshgrid!double meshgrid = buildMeshgrid(mesh);
	
	import core.thread : Thread;
	import std.datetime : dur;
	import std.format : format;
	
	if(plotU)
	{
		auto u = getVelocity!0(mesh);
		
		figure;
		
		contourf(meshgrid.X, meshgrid.Y, u, 350, `LineStyle`, `none`);
		colorbar;
		title(format("U velocity, t = %4.4f", t));
		xlabel("x");
		ylabel("y");
		if(outputFile != "")
		{
			print!"-dpdf"(outputFile~"_uvel");
		}
		Thread.sleep(100.dur!"msecs");
		
	}
	
	if(plotV)
	{
		auto v = getVelocity!1(mesh);
		
		figure;
		
		contourf(meshgrid.X, meshgrid.Y, v, 350, `LineStyle`, `none`);
		colorbar;
		title(format("V velocity, t = %4.4f", t));
		xlabel("x");
		ylabel("y");
		if(outputFile != "")
		{
			print!"-dpdf"(outputFile~"_vvel");
		}
		Thread.sleep(100.dur!"msecs");
	}
	
	if(plotP && !movie)
	{
		auto p = getPressure(mesh);
		
		figure;
		
		contourf(meshgrid.X, meshgrid.Y, p, 350, `LineStyle`, `none`);
		colorbar;
		caxis([0.9, 1.15]);
		title(format("Pressure, t = %4.4f", t));
		xlabel("x");
		ylabel("y");
		if(outputFile != "")
		{
			print!"-dpdf"(outputFile~"_pressure");
		}
		Thread.sleep(100.dur!"msecs");
	}
	
	if(plotRho)
	{
		auto rho = getDensity(mesh);
		
		figure;
		
		contourf(meshgrid.X, meshgrid.Y, rho, 350, `LineStyle`, `none`);
		colorbar;
		title(format("Density, t = %4.4f", t));
		xlabel("x");
		ylabel("y");
		if(outputFile != "")
		{
			print!"-dpdf"(outputFile~"_density");
		}
		Thread.sleep(100.dur!"msecs");
	}
	
	if(plotE)
	{
		auto e = getEnergy(mesh);
		
		figure;
		
		contourf(meshgrid.X, meshgrid.Y, e, 350, `LineStyle`, `none`);
		colorbar;
		title(format("Energy, t = %4.4f", t));
		xlabel("x");
		ylabel("y");
		if(outputFile != "")
		{
			print!"-dpdf"(outputFile~"_energy");
		}
		Thread.sleep(100.dur!"msecs");
	}
	
	if(plotM && !movie)
	{
		auto M = getMach(mesh);
		
		figure;
		
		contourf(meshgrid.X, meshgrid.Y, M, 350, `LineStyle`, `none`);
		colorbar;
		title(format("Mach number, t = %4.4f", t));
		xlabel("x");
		ylabel("y");
		if(outputFile != "")
		{
			print!"-dpdf"(outputFile~"_mach");
		}
		Thread.sleep(100.dur!"msecs");
	}
	
	if(plotSpeed)
	{
		auto s = getSpeed(mesh);
		
		figure;
		
		contourf(meshgrid.X, meshgrid.Y, s, 350, `LineStyle`, `none`);
		colorbar;
		title(format("Flow speed, t = %4.4f", t));
		xlabel("x");
		ylabel("y");
		if(outputFile != "")
		{
			print!"-dpdf"(outputFile~"_speed");
		}
		Thread.sleep(100.dur!"msecs");
	}
	
	if(movie)
	{
		mkdir("./output");
		auto dir = array(dirEntries("./", SpanMode.shallow));
		dir.sort!"a.timeLastModified < b.timeLastModified";
		
		import std.string : rightJustify;
		import std.conv : to;
		foreach(int i, dirFile; dir)
		{
			mesh = loadMesh(dirFile.name, dt, t);
			//auto p = getPressure(mesh);
			double[][] p;
			if(plotP)
			{
				p = getPressure(mesh);
			}
			else if(plotM)
			{
				p = getMach(mesh);
			}
		
			//figure;
			
			contourf(meshgrid.X, meshgrid.Y, p, 350, `LineStyle`, `none`);
			colorbar;
			//caxis([0.9, 1.2]);
			if(plotP)
			{
				caxis([0.5, 6]);
			}
			else if(plotM)
			{
				caxis([0, 3]);
			}
			
			axis("equal");
			title(format("Pressure, t = %4.4f", t));
			xlabel("x");
			ylabel("y");
			string strFrame = i.to!string.rightJustify(4, '0');
			print!"-dpng"(format("./output/%s.png", strFrame));
			Thread.sleep(100.dur!"msecs");
		}
	}
	writeln(file);
}