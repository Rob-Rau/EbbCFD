/+ Copyright (c) 2018 Robert F. Rau II +/
module ebb.mpid;

import std.conv;
import std.meta;
import std.stdio;
import std.traits;
import std.typecons;

import ebb.exception;
import ebb.solve;

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
mixin template GenPrintf(S...)
{
	private string genPrintf()
	{
		alias pfArgs = AliasSeq!S;
		string pfString = "@nogc void gpPrintf(S...)(S args)\n{\n\tprintf(\"";
		string argsString;
		foreach(i, pfArg; pfArgs)
		{
			pragma(msg, pfArg.stringof);
			static if(is(pfArg: int) && !is(pfArg: uint))
			{
				pfString ~= `%d`;
				argsString ~= `,args[`~i.to!string~`]`;
			}
			else static if(is(pfArg: uint))
			{
				pfString ~= `%u`;
				argsString ~= `,args[`~i.to!string~`]`;
			}
			else static if(is(pfArg: long) && !is(pfArg: ulong))
			{
				pfString ~= `%ld`;
				argsString ~= `,args[`~i.to!string~`]`;
			}
			else static if(is(pfArg: ulong))
			{
				pfString ~= `%llu`;
				argsString ~= `,args[`~i.to!string~`]`;
			}
			else static if(is(pfArg: string))
			{
				pfString ~= `%s`;
				argsString ~= `,args[`~i.to!string~`].ptr`;
			}
			else static if(isDynamicArray!pfArg)
			{
				pragma(msg, "dyn arr");
			}
		}

		pfString ~= `\n"`~argsString~`);
}`;
		return pfString;
	}
	pragma(msg, genPrintf);
	mixin(genPrintf);
}

void logln(S...)(S args)
{
	import core.stdc.stdio : fopen, fwrite, fopen, printf, snprintf;
	version(Have_mpi)
	{
		int id;
		MPI_Comm_rank(MPI_COMM_WORLD, &id);
		//string pfString = GenPrintf!(args);
		//mixin GenPrintf!(string, int, string, S);
		//gpPrintf("Proc ", id, ": ", args);
		//mixin GenPrintf!S;
		writeln("Proc ", id, ": ", args);
	}
	else
	{
		writeln("Proc 0: ", args);
	}
}

static @nogc void startall(Request[] requests)
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

static @nogc void waitall(Request[] requests, Status[] statuses)
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

struct Communication(SolverParams params)
{
	immutable vecSize = params.size;
	immutable dims = params.dims;

	immutable int commTag = 2000;

	private Comm comm;

	private Request[] sendRequests;
	private Request[] sendGradRequests;
	private Request[] recvRequests;
	private Request[] recvGradRequests;
	private Status[] statuses;

	private Vector!vecSize[][] sendStateBuffers;
	private Vector!vecSize[][] recvStateBuffers;

	private Matrix!(vecSize, dims)[][] sendGradBuffers;
	private Matrix!(vecSize, dims)[][] recvGradBuffers;

	@nogc void startSendRecvState(Vector!vecSize[] q, immutable(size_t[][]) idxMapping)
	{
		foreach(commIdx, idxArray; idxMapping)
		{
			foreach(i, idx; idxArray)
			{
				sendStateBuffers[commIdx][i] = q[idx];
			}
		}
		sendRequests.startall;
		recvRequests.startall;
	}

	@nogc void finishRecvState(ref Vector!vecSize[] q, immutable(size_t[][]) idxMapping)
	{
		recvRequests.waitall(statuses);
		foreach(commIdx, idxArray; idxMapping)
		{
			foreach(i, idx; idxArray)
			{
				q[idx] = recvStateBuffers[commIdx][i];
			}
		}
	}

	@nogc void startSendRecvStateGradient(Matrix!(vecSize, dims)[] dQ, immutable(size_t[][]) idxMapping)
	{
		foreach(commIdx, idxArray; idxMapping)
		{
			foreach(i, idx; idxArray)
			{
				sendGradBuffers[commIdx][i] = dQ[idx];
			}
		}
		
		sendGradRequests.startall;
		recvGradRequests.startall;
	}

	@nogc void finishRecvStateGradient(ref Matrix!(vecSize, dims)[] dQ, immutable(size_t[][]) idxMapping)
	{
		recvGradRequests.waitall(statuses);
		foreach(commIdx, idxArray; idxMapping)
		{
			foreach(i, idx; idxArray)
			{
				dQ[idx] = recvGradBuffers[commIdx][i];
			}
		}
	}

	static Datatype[string] customDataTypes;

	void abort(int reason)
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

	int size()
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

	int rank()
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

	template toMPIType(T)
	{
		static Datatype toMPIType()
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

	@nogc void barrier()
	{
		version(Have_mpi)
		{
			auto ret = MPI_Barrier(comm);
			enforce(ret == MPI_SUCCESS, "MPI_Barrier failed");
		}
	}

	void sendInit(T)(T[] buffer, int to, int tag, ref Request request)
	{
		version(Have_mpi)
		{
			auto ret = MPI_Send_init(buffer.ptr, cast(int)buffer.length, toMPIType!(T), to, tag, comm, &request);
			enforce(ret == MPI_SUCCESS, "MPI_Send_init failed with error: "~ret.to!string);
		}
	}

	@nogc void allreduce(T)(T send, ref T recv, Op op)
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

	@nogc void reduce(T)(T[] send, T[] recv, Op op, int root)
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

	void recvInit(T)(T[] buffer, int from, int tag, ref Request request)
	{
		version(Have_mpi)
		{
			auto ret = MPI_Recv_init(buffer.ptr, cast(int)buffer.length, toMPIType!(T), from, tag, comm, &request);
			enforce(ret == MPI_SUCCESS, "MPI_Recv_init failed with error: "~ret.to!string);
		}
	}

	void send(T, P)(T data, P proc, int tag)
		if(!isDynamicArray!T)
	{
		version(Have_mpi)
		{
			MPI_Send(&data, 1, toMPIType!T, proc.to!int, tag, comm);
		}
	}

	void send(T, P)(T[] data, P proc, int tag)
	{
		version(Have_mpi)
		{
			MPI_Send(data.ptr, data.length.to!int, toMPIType!T, proc.to!uint, tag, comm);
		}
	}

	void recv(T)(T[] data, uint from, int tag)
	{
		version(Have_mpi)
		{
			Status status;
			MPI_Recv(data.ptr, cast(int)data.length, toMPIType!T, from, tag, comm, &status);
		}
	}

	T recv(T, P)(P from, int tag)
		if(!isDynamicArray!T)
	{
		version(Have_mpi)
		{
			T data;
			Status status;
			MPI_Recv(&data, 1, toMPIType!T, from.to!int, tag, comm, &status);
			return data;
		}
		else
		{
			auto t = T;
			return t;
		}
	}

	Request irecv(T)(T[] data, uint from, uint tag)
	{
		Request request;
		version(Have_mpi)
		{
			auto ret = MPI_Irecv(data.ptr, cast(int)data.length, toMPIType!T, from, tag, comm, &request);
			enforce(ret == MPI_SUCCESS, "MPI_Irecv failed");
		}
		return request;
	}

	void probe(int source, int tag, ref Status status)
	{
		version(Have_mpi)
		{
			MPI_Probe(source, tag, comm, &status);
		}
	}

	void sendArray(T, P)(T[] data, P proc, int tag)
	{
		version(Have_mpi)
		{
			uint len = data.length.to!uint;
			MPI_Send(&len, 1, MPI_UINT32_T, proc.to!int, tag, comm);
			MPI_Send(data.ptr, data.length.to!int, toMPIType!T, proc.to!int, tag + 1, comm);
		}
	}

	T[] recvArray(T, P)(P from, const int tag)
	{
		version(Have_mpi)
		{
			int nElems;
			T[] data;

			Status status;
			MPI_Recv(&nElems, 1, MPI_UINT32_T, from.to!int, tag, comm, &status);
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

}
