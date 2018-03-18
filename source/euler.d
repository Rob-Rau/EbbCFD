/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.euler;

import std.math;

import numd.linearalgebra.matrix;

import ebb.mesh;

@nogc double getPressure(ref inout Vector!4 q, inout double gamma)
{
	return (gamma - 1)*(q[3] - 0.5*q[0]*((q[1]/q[0])^^2.0 + (q[2]/q[0])^^2.0));
}
