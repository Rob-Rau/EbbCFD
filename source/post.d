/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.post;

import std.algorithm : canFind, filter, sort;
import std.array : array, split;
import std.conv;
import std.file : dirEntries, DirEntry, mkdir, SpanMode;
import std.getopt;
import std.math;
import std.stdio;
import std.string;

import numd.linearalgebra.matrix;

import rpp.client.rpc;

import ebb.mesh;

static this()
{
	initRPP("127.0.0.1", 54000);
}

double[][] parseLimits(string limitsStr)
{
	double[][] axisLimits;
	auto lims = limitsStr.split(';');

	foreach(lim; lims)
	{
		auto parsedLims = lim.split(',');
		axisLimits ~= [parsedLims[0].strip.chompPrefix("[").to!double, parsedLims[1].strip.chomp("]").to!double];
	}

	return axisLimits;
}

int main(string[] args)
{
	string titleStr;
	string saveFile;
	string legendStr;
	string otherDataSet;
	string meshFile;
	string tag;
	string limits;
	double aoa;
	auto res = getopt(args, std.getopt.config.caseSensitive, std.getopt.config.bundling,
					  "t|title", "Plot title", &titleStr, "s|save", "File to save", &saveFile,
					  "l|legend", "Legend entries. separate entries with a semicolon", &legendStr,
					  "d|data2", "A second data set to plot", &otherDataSet,
					  "m|mesh", "Mesh file to load", &meshFile,
					  "a|tag", "Boundary tag to compute forces on", &tag,
					  "A|AOA", "Angle of attack of the flow", &aoa,
					  "i|limits", "Plot with axis limits. Use the following format: [xlow, xhigh];[ylow, yhigh]", &limits);

	if(res.helpWanted)
	{
		writeln("plotter options:");
		foreach(opt; res.options)
		{
			writeln(opt.optShort, " | ", opt.optLong, "\t\t", opt.help);
		}
		return 1;
	}

	if(meshFile == "")
	{
		writeln("Mesh file not supplied, exiting");
		return 1;
	}

	if(aoa.isNaN)
	{
		writeln("No supplied angle of attack, exiting");
		return 1;
	}

	UMesh2 mesh;
	if(meshFile.canFind(".gri"))
	{
		mesh = parseXflowMesh(meshFile);
	}
	else
	{
		writeln("Unsupported mesh type, exiting");
		return 1;
	}

	if(tag == "")
	{
		writeln("Boundary tag not supplied, exiting");
		return 1;
	}

	auto dir = array(dirEntries("./", SpanMode.shallow));
	dir = dir.filter!(a => a.name.canFind(".msln")).array.sort!("a.name.chompPrefix(\"./save_\").chomp(\".msln\").to!uint < b.name.chompPrefix(\"./save_\").chomp(\".msln\").to!uint").array;

	//auto mesh = loadMatlabSolution(meshFile, dir[0]);

	double[] lift = new double[dir.length];
	double[] time = new double[dir.length];

	double[] lift2, time2;
	DirEntry[] dir2;

	if(otherDataSet != "")
	{
		dir2 = array(dirEntries(otherDataSet, SpanMode.shallow));
		//dir2 = dir2.sort!"a.name < b.name".filter!(a => a.name.canFind(".mesh")).array;
		dir2 = dir2.filter!(a => a.name.canFind(".msln")).array.sort!("a.name.chompPrefix(\"./save_\").chomp(\".msln\").to!uint < b.name.chompPrefix(\"./save_\").chomp(\".msln\").to!uint").array;

		lift2 = new double[dir2.length];
		time2 = new double[dir2.length];
	}

	//double aoa = -45 * (PI/180);
	aoa *= (PI/180);

	auto rotMat = Matrix!(2, 2)(cos(aoa), -sin(aoa), sin(aoa), cos(aoa));

	foreach(int i, dirFile; dir)
	{
		double dt, t;
		//auto mesh = loadMesh(dirFile.name, dt, t);
		loadMatlabSolution(mesh, dirFile.name);

		Vector!2 f = rotMat*mesh.computeBoundaryForces(tag);
		lift[i] = f[1];
		time[i] = i;
	}

	if(otherDataSet != "")
	{
		foreach(int i, dirFile; dir2)
		{
			double dt, t;
			//auto mesh = loadMesh(dirFile.name, dt, t);
			loadMatlabSolution(mesh, dirFile.name);

			Vector!2 f = rotMat*mesh.computeBoundaryForces(tag);
			lift2[i] = f[1];
			time2[i] = i;
		}
	}

	figure;
	if(otherDataSet == "")
	{
		plot(time, lift);
	}
	else
	{
		plot(time, lift, time2, lift2);
	}

	if(limits != "")
	{
		auto lims = parseLimits(limits);
		axis([lims[0][0], lims[0][1], lims[1][0], lims[1][1]]);
	}
/+
	xlabel("save iteration");
	ylabel("lift");

	if(legendStr != "")
	{
		legend(legendStr.split(';'), "interpreter", "latex");
	}
+/
	if(legendStr != "")
	{
		setupPlot("save iteration", "lift", legendStr.split(';'), 12, "");
	}
	else
	{
		setupPlot("save iteration", "lift", ["data"], 12, "");
	}

	if(titleStr != "")
	{
		title(titleStr);
	}

	if(saveFile != "")
	{
		print!"-dpdf"(saveFile);
	}

	return 0;
}