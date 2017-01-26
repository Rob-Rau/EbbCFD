module ebb.mpid;

import std.conv;
import std.meta;
import std.stdio;
import std.traits;
import std.typecons;

import ebb.exception;

import numd.linearalgebra.matrix;
import numd.utility;

version(Have_mpi)
{
	public import mpi;
	public import mpi.util;

	alias Datatype = MPI_Datatype;
	alias Comm = MPI_Comm;
	alias Request = MPI_Request;
	alias Status = MPI_Status;
	alias Op = MPI_Op;
}
else
{
	struct Datatype {}
	struct Comm {}
	struct Request {}
	struct Status
	{
		int MPI_SOURCE;
        int MPI_TAG;
        int MPI_ERROR;
	}
	struct Op{}

	public Comm MPI_COMM_WORLD;
	public Comm MPI_COMM_SELF;
	public Op MPI_MIN;
	public Op MPI_MAX;
	public Op MPI_SUM;

	enum MPI_ANY_SOURCE = -1;
	import std.datetime;
	static StopWatch sw;
	static this()
	{
		sw.reset;
		sw.start;
	}
}

Datatype[string] customDataTypes;

void abort(Comm comm, int reason)
{
	version(Have_mpi)
	{
		MPI_Abort(comm, reason);
	}
}

void init(string[] args)
{
	version(Have_mpi)
	{
		int argc = cast(int)args.length;
		auto argv = args.toArgv;
		MPI_Init(&argc, &argv);
	}
}

int size(Comm comm)
{
	version(Have_mpi)
	{
		int p;
		MPI_Comm_size(comm, &p);
		return p;
	}
	else
	{
		return 1;
	}
}

double wtime()
{
	version(Have_mpi)
	{
		return MPI_Wtime();
	}
	else
	{
		return sw.peek.msecs/1000.0;
	}
}

int rank(Comm comm)
{
	version(Have_mpi)
	{
		int id;
		MPI_Comm_rank(MPI_COMM_WORLD, &id);
		return id;
	}
	else
	{
		return 0;
	}
}

int shutdown()
{
	version(Have_mpi)
	{
		pragma(msg, "have mpi");
		foreach(customDataType; customDataTypes.byValue)
		{
			MPI_Type_free(&customDataType);
		}
		return MPI_Finalize();
	}
	else
	{
		pragma(msg, "no mpi");
		return 0;
	}
}

void logln(S...)(S args)
{
	version(Have_mpi)
	{
		int id;
		MPI_Comm_rank(MPI_COMM_WORLD, &id);
		writeln("Proc ", id, ": ", args);
	}
	else
	{
		writeln("Proc 0: ", args);
	}
}

//MPI_Datatype toMPIType(T)
template toMPIType(T)
{
	Datatype toMPIType()
	{
		version(Have_mpi)
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

						Datatype newDataType;

						ptrdiff_t[] offsets;
						Datatype[] dataTypes;
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
		else
		{
			return Datatype();
		}
	}
}

@nogc void barrier(Comm comm)
{
	version(Have_mpi)
	{
		auto ret = MPI_Barrier(comm);
		enforce(ret == MPI_SUCCESS, "MPI_Barrier failed");
	}
}

void sendInit(T)(Comm comm, T[] buffer, int to, int tag, ref Request request)
{
	version(Have_mpi)
	{
		auto ret = MPI_Send_init(buffer.ptr, cast(int)buffer.length, toMPIType!(T), to, tag, comm, &request);
		enforce(ret == MPI_SUCCESS, "MPI_Send_init failed with error: "~ret.to!string);
	}
}

@nogc void startall(Request[] requests)
{
	version(Have_mpi)
	{
		if(requests.length > 0)
		{
			auto ret = MPI_Startall(cast(int)requests.length, requests.ptr);
			enforce(ret == MPI_SUCCESS, "MPI_Startall failed");
		}
	}
}

@nogc void waitall(Request[] requests, Status[] statuses)
{
	version(Have_mpi)
	{
		if(requests.length > 0)
		{
			auto ret = MPI_Waitall(cast(int)requests.length, requests.ptr, statuses.ptr);
			enforce(ret == MPI_SUCCESS, "MPI_Waitall failed");
		}
	}
}

@nogc void allreduce(T)(Comm comm, T send, ref T recv, Op op)
	if(is(T: double) || is(T: uint) || is(T: int))
{
	version(Have_mpi)
	{
		static if(is(T: double))
		{
			Datatype type = MPI_DOUBLE;
		}
		else static if(is(T: uint))
		{
			Datatype type = MPI_UINT32_T;
		}
		else static if(is(T: int))
		{
			Datatype type = MPI_INT32_T;
		}

		auto ret = MPI_Allreduce(&send, &recv, 1, type, op, comm);
		enforce(ret == MPI_SUCCESS, "MPI_AllReduce failed");
	}
}

@nogc void reduce(T)(Comm comm, T[] send, T[] recv, Op op, int root)
	if(is(T: double) || is(T: uint) || is(T: int))
{
	version(Have_mpi)
	{
		static if(is(T: double))
		{
			Datatype type = MPI_DOUBLE;
		}
		else static if(is(T: uint))
		{
			Datatype type = MPI_UINT32_T;
		}
		else static if(is(T: int))
		{
			Datatype type = MPI_INT32_T;
		}

		auto ret = MPI_Reduce(send.ptr, recv.ptr, cast(int)send.length, type, op, root, comm);
		enforce(ret == MPI_SUCCESS, "MPI_AllReduce failed");
	}
}

void recvInit(T)(Comm comm, T[] buffer, int from, int tag, ref Request request)
{
	version(Have_mpi)
	{
		auto ret = MPI_Recv_init(buffer.ptr, cast(int)buffer.length, toMPIType!(T), from, tag, comm, &request);
		enforce(ret == MPI_SUCCESS, "MPI_Recv_init failed with error: "~ret.to!string);
	}
}

void send(T)(Comm comm, T data, uint proc, int tag)
	if(!isDynamicArray!T)
{
	version(Have_mpi)
	{
		MPI_Send(&data, 1, toMPIType!T, proc, tag, comm);
	}
}

void send(T)(Comm comm, T[] data, uint proc, int tag)
{
	version(Have_mpi)
	{
		MPI_Send(data.ptr, cast(uint)data.length, toMPIType!T, proc, tag, comm);
	}
}

void recv(T)(Comm comm, T[] data, uint from, int tag)
{
	version(Have_mpi)
	{
		Status status;
		MPI_Recv(data.ptr, cast(int)data.length, toMPIType!T, from, tag, comm, &status);
	}
}

T recv(T)(Comm comm, uint from, int tag)
	if(!isDynamicArray!T)
{
	version(Have_mpi)
	{
		T data;
		Status status;
		MPI_Recv(&data, 1, toMPIType!T, from, tag, comm, &status);
		return data;
	}
	else
	{
		auto t = T;
		return t;
	}
}

Request irecv(T)(Comm comm, T[] data, uint from, uint tag)
{
	Request request;
	version(Have_mpi)
	{
		auto ret = MPI_Irecv(data.ptr, cast(int)data.length, toMPIType!T, from, tag, comm, &request);
		enforce(ret == MPI_SUCCESS, "MPI_Irecv failed");
	}
	return request;
}

void probe(Comm comm, int source, int tag, ref Status status)
{
	version(Have_mpi)
	{
		MPI_Probe(source, tag, comm, &status);
	}
}

void sendArray(T)(Comm comm, T[] data, uint proc, int tag)
{
	version(Have_mpi)
	{
		uint len = cast(uint)data.length;
		MPI_Send(&len, 1, MPI_UINT32_T, proc, tag, comm);
		MPI_Send(data.ptr, cast(uint)data.length, toMPIType!T, proc, tag + 1, comm);
	}
}

T[] recvArray(T)(Comm comm, uint from, const int tag)
{
	version(Have_mpi)
	{
		uint nElems;
		T[] data;

		Status status;
		MPI_Recv(&nElems, 1, MPI_UINT32_T, from, tag, comm, &status);
		data = new T[nElems];
		
		enforce(status.MPI_ERROR == MPI_SUCCESS, "Error receiving array length. MPI Error: "~status.MPI_ERROR.to!string);
		
		MPI_Recv(data.ptr, nElems, toMPIType!T, 0, tag + 1, comm, &status);

		return data;
	}
	else
	{
		T[] t;
		return t;
	}
}
