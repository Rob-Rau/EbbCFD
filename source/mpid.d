module ebb.mpid;

import std.conv;
import std.meta;
import std.stdio;
import std.traits;
import std.typecons;

import ebb.exception;

import numd.linearalgebra.matrix;
import numd.utility;

public import mpi;
public import mpi.util;

MPI_Datatype[string] customDataTypes;

int shutdown()
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

//MPI_Datatype toMPIType(T)
template toMPIType(T)
{
	MPI_Datatype toMPIType()
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
}

@nogc void barrier(MPI_Comm comm)
{
	auto ret = MPI_Barrier(comm);
	enforce(ret == MPI_SUCCESS, "MPI_Barrier failed");
}

void sendInit(T)(MPI_Comm comm, T[] buffer, int to, int tag, ref MPI_Request request)
{
	auto ret = MPI_Send_init(buffer.ptr, cast(int)buffer.length, toMPIType!(T), to, tag, comm, &request);
	enforce(ret == MPI_SUCCESS, "MPI_Send_init failed with error: "~ret.to!string);
}

@nogc void startall(MPI_Request[] requests)
{
	if(requests.length > 0)
	{
		auto ret = MPI_Startall(cast(int)requests.length, requests.ptr);
		enforce(ret == MPI_SUCCESS, "MPI_Startall failed");
	}
}

@nogc void waitall(MPI_Request[] requests, MPI_Status[] statuses)
{
	if(requests.length > 0)
	{
		auto ret = MPI_Waitall(cast(int)requests.length, requests.ptr, statuses.ptr);
		enforce(ret == MPI_SUCCESS, "MPI_Waitall failed");
	}
}

@nogc void allreduce(T)(MPI_Comm comm, T send, T recv, MPI_Op op)
{
	static if(__traits(compiles, ForeachType!(Unqual!T)))
	{
		alias uT = ForeachType!(Unqual!T);
	}
	else
	{
		alias uT = Unqual!T;
	}

	auto ret = MPI_Allreduce(&send, &recv, 1, MPI_DOUBLE, op, comm);
}

void recvInit(T)(MPI_Comm comm, T[] buffer, int from, int tag, ref MPI_Request request)
{
	auto ret = MPI_Recv_init(buffer.ptr, cast(int)buffer.length, toMPIType!(T), from, tag, comm, &request);
	enforce(ret == MPI_SUCCESS, "MPI_Recv_init failed with error: "~ret.to!string);
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
