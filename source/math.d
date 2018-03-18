/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.math;

import std.math;

import ebb.config;
import ebb.exception;
import ebb.mesh;
import ebb.mpid;

import numd.utility;
import numd.linearalgebra.matrix;

class Solution(size_t size)
{
	double[] dts;
	double dt;
	Vector!size[] dQdt;
	Vector!size[] Q;
	ulong iteration;
}

@nogc Vector!4 abs(Vector!4 data)
{
	return Vector!4(std.math.abs(data[0]), std.math.abs(data[1]), std.math.abs(data[2]), std.math.abs(data[3]));
}

@nogc S sum(string member = "", string op = "", S = double, T)(T[] data)
{
	if(data.length <= 6)
	{
		mixin("S s = "~op~"(data[0]"~member~");");
		for(uint i = 1; i < data.length; i++)
		{
			mixin("s += "~op~"(data[i]"~member~");");
		}
		return s;
	}
	else
	{
		uint m = cast(uint)floor(cast(double)data.length/2.0);
		return sum!(member, op, S, T)(data[0..m-1]) + sum!(member, op, S, T)(data[m..$-1]);
	}
}
