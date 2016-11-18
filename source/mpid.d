module ebb.mpid;

import std.conv;
import std.exception;
import std.meta;
import std.stdio;
import std.traits;

import numd.linearalgebra.matrix;
import numd.utility;

public import mpi;
public import mpi.util;

MPI_Datatype[string] customDataTypes;

int MPI_Shutdown()
{
	foreach(customDataType; customDataTypes.byValue)
	{
		MPI_Type_free(&customDataType);
	}
	return MPI_Finalize();
}

void logln(S...)(S args)
{
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	writeln("Proc ", id, ": ", args);
}

MPI_Datatype toMPIType(T)()
{
	static if(!is(T : char))
	{
		static if(is(T : int) && !is(T: uint))
		{
			return MPI_INT32_T;
		}
		else static if(is(T : uint))
		{
			return MPI_UINT32_T;
		}
		else static if(is(T : long) && !is(T : ulong))
		{
			return MPI_INT64_T;
		}
		else static if(is(T : ulong))
		{
			return MPI_UINT64_T;
		}
		else static if(is(T : float) && !is(T : double))
		{
			return MPI_FLOAT;
		}
		else static if(is(T : double))
		{
			return MPI_DOUBLE;
		}
		else
		{
			auto customType = (T.stringof in customDataTypes); 
			if(customType is null)
			{
				T inst;

				MPI_Datatype newDataType;

				long[] offsets;
				MPI_Datatype[] dataTypes;
				int[] blocklengths;

				foreach(uint i, member; FieldNameTuple!T)
				{
					static assert(!isDynamicArray!(Fields!T[i]), "Dynamic arrays not supported");
					static if(isArray!(Fields!T[i]))
					{
						mixin("blocklengths ~= inst."~member~".length;");
						mixin("dataTypes ~= toMPIType!("~ForeachType!(Fields!T[i]).stringof~");");
						mixin("offsets ~= cast(uint)(cast(void*)&inst."~member~"[0] - cast(void*)&inst);");
					}
					else
					{
						blocklengths ~= 1;
						mixin("dataTypes ~= toMPIType!("~Fields!T[i].stringof~");");
						mixin("offsets ~= cast(uint)(cast(void*)&inst."~member~" - cast(void*)&inst);");
					}
				}

				MPI_Type_create_struct(cast(int)FieldNameTuple!T.length, blocklengths.ptr, offsets.ptr, dataTypes.ptr, &newDataType);
				MPI_Type_commit(&newDataType);

				customDataTypes[T.stringof] = newDataType;
				return newDataType;
			}
			else
			{
				return customDataTypes[T.stringof];
			}
		}
	}
	else
	{
		return MPI_CHAR;
	}
}

void send(T)(MPI_Comm comm, T data, uint proc, int tag)
{
	MPI_Send(&data, 1, toMPIType!T, proc, tag, comm);
}

T recv(T)(MPI_Comm comm, uint from, int tag)
{
	T data;
	MPI_Status status;
	MPI_Recv(&data, 1, toMPIType!T, from, tag, comm, &status);
	return data;
}

void sendArray(T)(MPI_Comm comm, T[] data, uint proc, int tag)
{
	uint len = cast(uint)data.length;
	MPI_Send(&len, 1, MPI_UINT32_T, proc, tag, comm);
	MPI_Send(data.ptr, cast(uint)data.length, toMPIType!T, proc, tag + 1, comm);
}

T[] recvArray(T)(MPI_Comm comm, uint from, const int tag)
{
	uint nElems;
	T[] data;

	MPI_Status status;
	MPI_Recv(&nElems, 1, MPI_UINT32_T, from, tag, comm, &status);
	data = new T[nElems];
	
	enforce(status.MPI_ERROR == MPI_SUCCESS, "Error receiving array length. MPI Error: "~status.MPI_ERROR.to!string);
	
	MPI_Recv(data.ptr, nElems, toMPIType!T, 0, tag + 1, comm, &status);

	return data;
}
