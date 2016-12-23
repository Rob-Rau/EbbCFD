module ebb.exception;

import std.algorithm;
import std.conv;
import std.exception;
import std.stdio;
import std.meta;
import std.string;
import std.typecons;
import std.traits;

import ebb.mesh;

import numd.linearalgebra.matrix;

class CellException : Exception
{
	this(string msg)
	{
		super(msg);
	}

	UCell2 cell;	
}

class EdgeException : Exception
{
	this(string msg)
	{
		super(msg);
	}

	Edge edge;
}

private alias Exceptions = Tuple!(string, Exception)[];
private Exceptions exceptions;
private Exceptions addExceptions(alias mod)()
{
	Exceptions newExceptions;
	foreach(i, member; __traits(allMembers, mod))
    {
        static if(__traits(compiles, __traits(classInstanceSize, mixin(member))))
        {
            foreach(base; BaseClassesTuple!(mixin(member)))
            {
                static if(is(base : Exception))
                {
                    pragma(msg, "Class "~member~" is an exception");
                    mixin("newExceptions ~= tuple("~member~".stringof, cast(Exception)new "~member~"(\"\"));");
					break;
                }
            }
        }
    }
	return newExceptions;
}

static this()
{
	exceptions ~= tuple(Exception.stringof, new Exception(""));
	exceptions ~= addExceptions!(ebb.exception);
	exceptions ~= addExceptions!(std.exception);
	exceptions ~= addExceptions!(std.conv);
}

@nogc void enforce(T = Exception)(bool cond, string msg, string file = __FILE__, int line = __LINE__)
{
	if(!cond)
	{
		auto idx = exceptions.countUntil!"a[0].equal(b)"(T.stringof);
		if(idx >= 0)
		{
			auto ex = cast(T)exceptions[idx][1];
			ex.msg = msg;
			throw ex;
		}
		else
		{
			auto ex = cast(Exception)exceptions[0][1];
			ex.msg = msg;
			ex.line = line;
			ex.file = file;
			throw ex;
		}
	}
}

@nogc void enforce(T = Exception)(bool cond, string msg, UCell2 cell, string file = __FILE__, int line = __LINE__)
	if(is(T : CellException))
{
	if(!cond)
	{
		auto idx = exceptions.countUntil!"a[0].equal(b)"(T.stringof);
		if(idx >= 0)
		{
			auto ex = cast(T)exceptions[idx][1];
			ex.msg = msg;
			ex.cell = cell;
			ex.line = line;
			ex.file = file;
			throw ex;
		}
		else
		{
			auto ex = cast(Exception)exceptions[0][1];
			ex.msg = msg;
			ex.line = line;
			ex.file = file;
			throw ex;
		}
	}
}

@nogc void enforce(T = Exception)(bool cond, string msg, Edge edge, string file = __FILE__, int line = __LINE__)
	if(is(T : EdgeException))
{
	if(!cond)
	{
		auto idx = exceptions.countUntil!"a[0].equal(b)"(T.stringof);
		if(idx >= 0)
		{
			auto ex = cast(T)exceptions[idx][1];
			ex.msg = msg;
			ex.edge = edge;
			ex.line = line;
			ex.file = file;
			throw ex;
		}
		else
		{
			auto ex = cast(Exception)exceptions[0][1];
			ex.msg = msg;
			ex.line = line;
			ex.file = file;
			throw ex;
		}
	}
}