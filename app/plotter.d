/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.plotter;

import core.thread : Thread;

import std.conv;
import std.datetime : dur;
import std.format : format;
import std.getopt;
import std.meta : aliasSeqOf;
import std.stdio;
import std.string;

import rpp.client.rpc;

import numd.utility;

import ebb.mesh;
import ebb.euler;

static this()
{
	initRPP("127.0.0.1", 54000);
}

void plotFunc(alias func)(Mesh mesh, ref Meshgrid!double meshgrid, double t, double[] caxisLims, string titleStr, string outputFile, string suffix)
{
	double[][] value = func(mesh);
	figure;
	
	contourf(meshgrid.X, meshgrid.Y, value, 350, `LineStyle`, `none`);
	colorbar;
	
	if(caxisLims.length > 0)
	{
		caxis(caxisLims);
	}

	title(format("%s, t = %4.4f", titleStr, t));
	xlabel("x");
	ylabel("y");
	axis("equal");

	if(outputFile != "")
	{
		print!"-dpdf"(outputFile~"_"~suffix);
	}

	Thread.sleep(100.dur!"msecs");
}

double[][] parseCaxis(string caxisStr)
{
	double[][] caxisLimits;
	auto clims = caxisStr.split(';');

	foreach(clim; clims)
	{
		auto parsedLims = clim.split(',');
		caxisLimits ~= [parsedLims[0].strip.chompPrefix("[").to!double, parsedLims[1].strip.chomp("]").to!double];
	}

	return caxisLimits;
}

void main(string[] args)
{
	string file, outputFile, caxisString;
	bool plotU = false, plotV = false, plotP = false, plotRho = false, plotE = false, plotM = false, plotSpeed = false, movie = false;
	
	auto res = getopt(args, std.getopt.config.caseSensitive, std.getopt.config.bundling,
							"file|f", "file to plot", &file, "uvel|u", "plot u velocity component", &plotU, "vvel|v", "plot v velocity component", &plotV,
							"pressure|p", "plot pressure", &plotP, "density|d", "plot density", &plotRho, "energy|e", "plot energy", &plotE,
							"Mach|M", "plot Mach number", &plotM, "speed|s", "plot flow speed", &plotSpeed,
							"movie|m", "generate output images for movie", &movie, "output|o", "output file name", &outputFile,
							"caxis|c", "Set the colorbar limits. Use the following format: [p1low, p1high];[p2low, p2high]", &caxisString);

	auto caxisLims = parseCaxis(caxisString);
	uint caIdx = 0;
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
	
	if(file == "" && !movie)
	{
		writeln("no input file, exiting");
		return;
	}

	Mesh mesh;
	double dt, t;
	mesh = loadMesh(file, dt, t);
	Meshgrid!double meshgrid = buildMeshgrid(mesh);
	
	if(plotU)
	{
		double[] caLim;
		if(caIdx < caxisLims.length)
		{
			caLim = caxisLims[caIdx];
			caIdx++;
		}
		plotFunc!(getVelocity!0)(mesh, meshgrid, t, caLim, "U velocity", outputFile, "uvel");
	}
	
	if(plotV)
	{
		double[] caLim;
		if(caIdx < caxisLims.length)
		{
			caLim = caxisLims[caIdx];
			caIdx++;
		}
		plotFunc!(getVelocity!1)(mesh, meshgrid, t, caLim, "V velocity", outputFile, "vvel");
	}
	
	if(plotP && !movie)
	{
		double[] caLim;
		if(caIdx < caxisLims.length)
		{
			caLim = caxisLims[caIdx];
			caIdx++;
		}
		plotFunc!(getPressure)(mesh, meshgrid, t, caLim, "Pressure", outputFile, "pressure");
	}
	
	if(plotRho)
	{
		double[] caLim;
		if(caIdx < caxisLims.length)
		{
			caLim = caxisLims[caIdx];
			caIdx++;
		}
		plotFunc!(getDensity)(mesh, meshgrid, t, caLim, "Density", outputFile, "density");
	}
	
	if(plotE)
	{
		double[] caLim;
		if(caIdx < caxisLims.length)
		{
			caLim = caxisLims[caIdx];
			caIdx++;
		}
		plotFunc!(getEnergy)(mesh, meshgrid, t, caLim, "Energy", outputFile, "energy");
	}
	
	if(plotM && !movie)
	{
		double[] caLim;
		if(caIdx < caxisLims.length)
		{
			caLim = caxisLims[caIdx];
			caIdx++;
		}
		plotFunc!(getMach)(mesh, meshgrid, t, caLim, "Mach number", outputFile, "mach");
	}
	
	if(plotSpeed && !movie)
	{
		double[] caLim;
		if(caIdx < caxisLims.length)
		{
			caLim = caxisLims[caIdx];
			caIdx++;
		}
		plotFunc!(getSpeed)(mesh, meshgrid, t, caLim, "Flow speed", outputFile, "speed");
	}
	
	if(movie)
	{
		auto dir = array(dirEntries("./", SpanMode.shallow));
		mkdir("./output");
		dir.sort!"a.name < b.name";
		
		string titleStr = "";
		import std.string : rightJustify;
		import std.conv : to;
		foreach(int i, dirFile; dir)
		{
			mesh = loadMesh(dirFile.name, dt, t);
			double[][] value;
			if(plotP)
			{
				value = getPressure(mesh);
				titleStr = "Pressure";
			}
			else if(plotM)
			{
				value = getMach(mesh);
				titleStr = "Mach Number";
			}
			else if(plotSpeed)
			{
				value = getSpeed(mesh);
				titleStr = "Flow Speed";
			}
			
			contourf(meshgrid.X, meshgrid.Y, value, 350, `LineStyle`, `none`);
			colorbar;

			caxis(caxisLims[0]);

			axis("equal");
			title(format("%s, t = %4.4f", titleStr, t));
			xlabel("x");
			ylabel("y");
			string strFrame = i.to!string.rightJustify(7, '0');
			print!"-dpng"(format("./output/%s.png", strFrame));
			Thread.sleep(100.dur!"msecs");
		}
	}
}