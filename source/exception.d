module ebb.exception;

import ebb.mesh;

import numd.linearalgebra.matrix;

class SolverException : Exception
{
	this(string msg)
	{
		super(msg);
	}

	enum SExceptionType
	{
		EdgeException,
		CellException
	}

	@nogc struct EdgeException
	{
		double pL, pR;
		Vector!4 flux, qL, qR;
		Vector!2 normal;
		uint cellL, cellR;

		@nogc this(double pL, double pR, Vector!4 flux, Vector!4 qL, Vector!4 qR, Vector!2 n, uint cellL, uint cellR)
		{
			this.pL = pL;
			this.pR = pR;
			this.flux = flux;
			this.qL = qL;
			this.qR = qR;
			this.normal = n;
			this.cellL = cellL;
			this.cellR = cellR;
		}
	}

	@nogc struct CellException
	{
		double p;
		Vector!4 q;
		uint cellIdx;

		@nogc this(double p, Vector!4 q, uint cellIdx)
		{
			this.p = p;
			this.q = q;
			this.cellIdx = cellIdx;
		}
	}

	SExceptionType exceptionType;
	CellException cExcept;
	EdgeException eExcept;

	@nogc void SetException(SExceptionType type, string msg, EdgeException ex, string file = __FILE__, size_t line = __LINE__)
	{
		exceptionType = type;
		eExcept = ex;
		this.msg = msg;
		this.file = file;
		this.line = line;
	}

	@nogc void SetException(SExceptionType type, string msg, CellException ex, string file = __FILE__, size_t line = __LINE__)
	{
		exceptionType = type;
		cExcept = ex;
		this.msg = msg;
		this.file = file;
		this.line = line;
	}

	/+
	printf("pL = %f\n", getPressure(mesh.q[i]));
	printf("qL = [%f, %f, %f, %f]\n", mesh.q[i][0], mesh.q[i][1], mesh.q[i][2], mesh.q[i][3]);
	printf("cell = %d\n", i);
	+/
	/+
	printf("pL = %f\n", getPressure(mesh.edges[i].q[0]));
	printf("pR = %f\n", getPressure(mesh.edges[i].q[1]));
	printf("Flux = [%f, %f, %f, %f]\n", mesh.edges[i].flux[0], mesh.edges[i].flux[1], mesh.edges[i].flux[2], mesh.edges[i].flux[3]);
	printf("qL = [%f, %f, %f, %f]\n", mesh.edges[i].q[0][0], mesh.edges[i].q[0][1], mesh.edges[i].q[0][2], mesh.edges[i].q[0][3]);
	printf("qR = [%f, %f, %f, %f]\n", mesh.edges[i].q[1][0], mesh.edges[i].q[1][1], mesh.edges[i].q[1][2], mesh.edges[i].q[1][3]);
	printf("cell L = %d\n", mesh.edges[i].cellIdx[0]);
	printf("cell R = %d\n", mesh.edges[i].cellIdx[1]);
	printf("normal = [%f, %f]\n", mesh.edges[i].normal[0], mesh.edges[i].normal[1]);
	+/
}