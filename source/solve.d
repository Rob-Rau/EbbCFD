/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.solve;

import core.atomic;
import core.stdc.stdio : fopen, fwrite, fopen, printf, snprintf;
import core.sys.posix.signal;

import std.algorithm : canFind, countUntil, min, map, max, reduce, sum;
import std.conv;
import std.file;
import std.math;
import std.meta : aliasSeqOf;
import std.stdio;
import std.typecons;

import numd.utility;
import numd.linearalgebra.matrix;

import ebb.config;
import ebb.euler;
import ebb.exception;
public import ebb.solvers.finitevolume;
import ebb.flux;
import ebb.integrate;
import ebb.mesh;
import ebb.io;
import ebb.mpid;

alias solverList = aliasSeqOf!(["FiniteVolume"]);

enum PhysicalDimensions: size_t
{
	TwoD = 2,
	ThreeD = 3
}

alias SolverParams = Tuple!(size_t, "size", PhysicalDimensions, "dims");
alias solverParams = tuple!("size", "dims");

alias FluxResult(size_t size) = Tuple!(Vector!size, "flux", double, "waveSpeed");
alias fluxResult = tuple!("flux", "waveSpeed");
