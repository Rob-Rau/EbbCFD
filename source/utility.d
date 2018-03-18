/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.utility;

import std.algorithm;
import std.conv;
import std.exception;
import std.math;
import std.meta;
import std.range;
import std.stdio;
import std.string;
import std.typecons;

string switchBuilder(int level, string switchVar, Args...)(string statement)
{
	alias list = AliasSeq!Args;

	string fillPlaceHolder(int level)(in string statement, string arg)
	{
		auto strSlice = statement;
		ptrdiff_t searchIdx = 0;
		string newStatement = statement;
		while(searchIdx < statement.length)
		{
			auto idxStart = newStatement.indexOf('{', searchIdx);
			if(idxStart == -1)
			{
				break;
			}
			
			auto idxEnd = newStatement.indexOf('}', idxStart);
			assert(idxEnd > idxStart);
			
			auto sliceLen = idxEnd - (idxStart + 1);
			if(newStatement[idxStart+1..idxStart+1 + sliceLen].isNumeric)
			{
				auto strLevel = newStatement[idxStart+1..idxStart+1 + sliceLen].to!int;
				if(strLevel == level)
				{
					newStatement = newStatement[0..idxStart] ~ arg.strip('\"') ~ newStatement[idxEnd + 1..$];
				}
			}

			searchIdx = idxStart + 1;
		}
		
		return newStatement;
	}

	string switchStatement =
	"final switch("~switchVar~")\n"~
	"{\n";
	foreach(arg; list)
	{
		auto thisStatement = fillPlaceHolder!level(statement, arg.stringof);
		switchStatement ~= 
			"\tcase "~arg.stringof~":\n"~
				"\t\t"~thisStatement~"\n"~
				"\t\tbreak;";
	}
		
	switchStatement ~=
	"}";
	return switchStatement;
}
