/+ Copyright (c) 2016 Robert F. Rau II +/
module ebb.post;

import std.algorithm : canFind, filter, sort;
import std.array : array, split;
import std.file : dirEntries, DirEntry, mkdir, SpanMode;
import std.getopt;
import std.math;
import std.stdio;

import numd.linearalgebra.matrix;

import rpp.client.rpc;

import ebb.mesh;

static this()
{
	initRPP("127.0.0.1", 54000);
}

void main(string[] args)
{
	string titleStr;
	string saveFile;
	string legendStr;
	string otherDataSet;
	auto res = getopt(args, std.getopt.config.caseSensitive, std.getopt.config.bundling,
					  "t|title", "Plot title", &titleStr, "s|save", "File to save", &saveFile,
					  "l|legend", "Legend entries. separate entries with a semicolon", &legendStr,
					  "d|data2", "A second data set to plot", &otherDataSet);

	if(res.helpWanted)
	{
		writeln("plotter options:");
		foreach(opt; res.options)
		{
			writeln(opt.optShort, " | ", opt.optLong, "\t\t", opt.help);
		}
		return;
	}

	
	auto dir = array(dirEntries("./", SpanMode.shallow));
	dir = dir.sort!"a.name < b.name".filter!(a => a.name.canFind(".mesh")).array;

	double[] lift = new double[dir.length];
	double[] time = new double[dir.length];

	double[] lift2, time2;
	DirEntry[] dir2;

	if(otherDataSet != "")
	{
		dir2 = array(dirEntries(otherDataSet, SpanMode.shallow));
		dir2 = dir2.sort!"a.name < b.name".filter!(a => a.name.canFind(".mesh")).array;

		lift2 = new double[dir2.length];
		time2 = new double[dir2.length];
	}

	double aoa = -45 * (PI/180);

	auto rotMat = Matrix!(2, 2)(cos(aoa), -sin(aoa), sin(aoa), cos(aoa));

	foreach(int i, dirFile; dir)
	{
		double dt, t;
		auto mesh = loadMesh(dirFile.name, dt, t);

		Vector!2 f = rotMat*mesh.computeBodyForces();
		lift[i] = f[0];
		time[i] = t;
	}

	if(otherDataSet != "")
	{
		foreach(int i, dirFile; dir2)
		{
			double dt, t;
			auto mesh = loadMesh(dirFile.name, dt, t);

			Vector!2 f = rotMat*mesh.computeBodyForces();
			lift2[i] = f[0];
			time2[i] = t;
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

	xlabel("t [s]");
	ylabel("lift");

	if(legendStr != "")
	{
		legend(legendStr.split(';'), "interpreter", "latex");
	}

	if(titleStr != "")
	{
		title(titleStr);
	}

	if(saveFile != "")
	{
		print!"-dpdf"(saveFile);
	}
}