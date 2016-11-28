/* Converted to D from parmetis.h by htod */
module parmetis;
/*
 * Copyright 1997-2003, Regents of the University of Minnesota
 *
 * parmetis.h
 *
 * This file contains function prototypes and constrant definitions for 
 * ParMETIS
 *
 * Started 7/21/03
 * George
 *
 */

extern (C):

import mpi;

/*************************************************************************
* Data-structures
**************************************************************************/
/* Undefine the following #define in order to use short int as the idxtype */

/* Indexes are as long as integers for now */
alias int idxtype;

/*************************************************************************
* Constants 
**************************************************************************/
const PARMETIS_MAJOR_VERSION = 3;
const PARMETIS_MINOR_VERSION = 1;

/*************************************************************************
* Function prototypes
**************************************************************************/

/*-------------------------------------------------------------------
* API Introduced with Release 3.0 (current API) 
*--------------------------------------------------------------------*/
void  ParMETIS_V3_AdaptiveRepart(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, float *tpwgts, float *ubvec, float *ipc2redist, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

void  ParMETIS_V3_PartGeomKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *ncon, int *nparts, float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

void  ParMETIS_V3_PartGeom(idxtype *vtxdist, int *ndims, float *xyz, idxtype *part, MPI_Comm *comm);

void  ParMETIS_V3_PartKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

void  ParMETIS_V3_Mesh2Dual(idxtype *elmdist, idxtype *eptr, idxtype *eind, int *numflag, int *ncommonnodes, idxtype **xadj, idxtype **adjncy, MPI_Comm *comm);

void  ParMETIS_V3_PartMeshKway(idxtype *elmdist, idxtype *eptr, idxtype *eind, idxtype *elmwgt, int *wgtflag, int *numflag, int *ncon, int *ncommonnodes, int *nparts, float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

void  ParMETIS_V3_NodeND(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm);

void  ParMETIS_V3_RefineKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);



/*------------------------------------------------------------------
* Backward compatibility routines with Release 2.0
*-------------------------------------------------------------------*/
void  ParMETIS_PartKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

void  ParMETIS_PartGeomKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *nparts, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

void  ParMETIS_PartGeom(idxtype *vtxdist, int *ndims, float *xyz, idxtype *part, MPI_Comm *comm);

void  ParMETIS_PartGeomRefine(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

void  ParMETIS_RefineKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

void  ParMETIS_RepartLDiffusion(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

void  ParMETIS_RepartGDiffusion(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

void  ParMETIS_RepartRemap(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

void  ParMETIS_RepartMLRemap(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

void  ParMETIS_NodeND(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm);

void  ParMETIS_SerialNodeND(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm);




/*-------------------------------------------------------------------
* Backward compatibility routines with Release 1.0 
*--------------------------------------------------------------------*/
void  PARKMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, idxtype *part, int *options, MPI_Comm comm);

void  PARGKMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, int ndims, float *xyz, idxtype *part, int *options, MPI_Comm comm);

void  PARGRMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, int ndims, float *xyz, idxtype *part, int *options, MPI_Comm comm);

void  PARGMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int ndims, float *xyz, idxtype *part, int *options, MPI_Comm comm);

void  PARRMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, idxtype *part, int *options, MPI_Comm comm);

void  PARUAMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, idxtype *part, int *options, MPI_Comm comm);

void  PARDAMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, idxtype *part, int *options, MPI_Comm comm);
