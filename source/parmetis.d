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

//C     #ifndef __parmetis_h__
//C     #define __parmetis_h__

//#include <mpi.h>
//#define MPI_Comm void
//C     typedef void MPI_Comm;
extern (C):

import mpi;
//C     #ifndef _MSC_VER
//C     #define __cdecl
//C     #endif


/*************************************************************************
* Data-structures
**************************************************************************/
/* Undefine the following #define in order to use short int as the idxtype */
//C     #define IDXTYPE_INT

/* Indexes are as long as integers for now */
//C     #ifdef IDXTYPE_INT
//C     typedef int idxtype;
alias int idxtype;
//C     #else
//C     typedef short idxtype;
//C     #endif


/*************************************************************************
* Constants 
**************************************************************************/
//C     #define PARMETIS_MAJOR_VERSION        3
//C     #define PARMETIS_MINOR_VERSION        1
const PARMETIS_MAJOR_VERSION = 3;

const PARMETIS_MINOR_VERSION = 1;

/*************************************************************************
* Function prototypes
**************************************************************************/
//C     #ifdef __cplusplus
//C     extern "C" {
//C     #endif

/*-------------------------------------------------------------------
* API Introduced with Release 3.0 (current API) 
*--------------------------------------------------------------------*/
//C     void __cdecl ParMETIS_V3_AdaptiveRepart(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
//C     	     idxtype *vsize, idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, 
//C     	     int *nparts, float *tpwgts, float *ubvec, float *ipc2redist, 
//C     	     int *options, int *edgecut, idxtype *part, MPI_Comm *comm);
void  ParMETIS_V3_AdaptiveRepart(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, float *tpwgts, float *ubvec, float *ipc2redist, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

//C     void __cdecl ParMETIS_V3_PartGeomKway(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
//C     	     idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, 
//C     	     int *ncon, int *nparts, float *tpwgts, float *ubvec, int *options, 
//C     	     int *edgecut, idxtype *part, MPI_Comm *comm);
void  ParMETIS_V3_PartGeomKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *ncon, int *nparts, float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

//C     void __cdecl ParMETIS_V3_PartGeom(
//C                  idxtype *vtxdist, int *ndims, float *xyz, idxtype *part, MPI_Comm *comm);
void  ParMETIS_V3_PartGeom(idxtype *vtxdist, int *ndims, float *xyz, idxtype *part, MPI_Comm *comm);

//C     void __cdecl ParMETIS_V3_PartKway(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
//C     	     idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, 
//C     	     float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, 
//C     	     MPI_Comm *comm);
void  ParMETIS_V3_PartKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

//C     void __cdecl ParMETIS_V3_Mesh2Dual(
//C                  idxtype *elmdist, idxtype *eptr, idxtype *eind, int *numflag, 
//C     	     int *ncommonnodes, idxtype **xadj, idxtype **adjncy, MPI_Comm *comm);
void  ParMETIS_V3_Mesh2Dual(idxtype *elmdist, idxtype *eptr, idxtype *eind, int *numflag, int *ncommonnodes, idxtype **xadj, idxtype **adjncy, MPI_Comm *comm);

//C     void __cdecl ParMETIS_V3_PartMeshKway(
//C                  idxtype *elmdist, idxtype *eptr, idxtype *eind, idxtype *elmwgt, 
//C     	     int *wgtflag, int *numflag, int *ncon, int *ncommonnodes, int *nparts, 
//C     	     float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, 
//C     	     MPI_Comm *comm);
void  ParMETIS_V3_PartMeshKway(idxtype *elmdist, idxtype *eptr, idxtype *eind, idxtype *elmwgt, int *wgtflag, int *numflag, int *ncon, int *ncommonnodes, int *nparts, float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

//C     void __cdecl ParMETIS_V3_NodeND(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, 
//C     	     int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm);
void  ParMETIS_V3_NodeND(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm);

//C     void __cdecl ParMETIS_V3_RefineKway(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
//C     	     idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, 
//C     	     float *tpwgts, float *ubvec, int *options, int *edgecut, 
//C     	     idxtype *part, MPI_Comm *comm);
void  ParMETIS_V3_RefineKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);



/*------------------------------------------------------------------
* Backward compatibility routines with Release 2.0
*-------------------------------------------------------------------*/
//C     void __cdecl ParMETIS_PartKway(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
//C     	     idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, 
//C                  int *edgecut, idxtype *part, MPI_Comm *comm);
void  ParMETIS_PartKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

//C     void __cdecl ParMETIS_PartGeomKway(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, 
//C     	     int *wgtflag, int *numflag, int *ndims, float *xyz, int *nparts, int *options, 
//C     	     int *edgecut, idxtype *part, MPI_Comm *comm);
void  ParMETIS_PartGeomKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *nparts, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

//C     void __cdecl ParMETIS_PartGeom(
//C                  idxtype *vtxdist, int *ndims, float *xyz, idxtype *part, MPI_Comm *comm);
void  ParMETIS_PartGeom(idxtype *vtxdist, int *ndims, float *xyz, idxtype *part, MPI_Comm *comm);

//C     void __cdecl ParMETIS_PartGeomRefine(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
//C     	     idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, 
//C     	     int *options, int *edgecut, idxtype *part, MPI_Comm *comm);
void  ParMETIS_PartGeomRefine(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

//C     void __cdecl ParMETIS_RefineKway(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
//C     	     idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, 
//C     	     idxtype *part, MPI_Comm *comm);
void  ParMETIS_RefineKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

//C     void __cdecl ParMETIS_RepartLDiffusion(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
//C     	     idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, 
//C     	     idxtype *part, MPI_Comm *comm);
void  ParMETIS_RepartLDiffusion(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

//C     void __cdecl ParMETIS_RepartGDiffusion(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
//C     	     idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, 
//C     	     idxtype *part, MPI_Comm *comm);
void  ParMETIS_RepartGDiffusion(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

//C     void __cdecl ParMETIS_RepartRemap(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, 
//C     	     int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, 
//C     	     MPI_Comm *comm);
void  ParMETIS_RepartRemap(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

//C     void __cdecl ParMETIS_RepartMLRemap(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, 
//C     	     int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, 
//C     	     MPI_Comm *comm);
void  ParMETIS_RepartMLRemap(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm);

//C     void __cdecl ParMETIS_NodeND(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, 
//C     	     idxtype *order, idxtype *sizes, MPI_Comm *comm);
void  ParMETIS_NodeND(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm);

//C     void __cdecl ParMETIS_SerialNodeND(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, 
//C     	     idxtype *order, idxtype *sizes, MPI_Comm *comm);
void  ParMETIS_SerialNodeND(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm);




/*-------------------------------------------------------------------
* Backward compatibility routines with Release 1.0 
*--------------------------------------------------------------------*/
//C     void __cdecl PARKMETIS(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
//C     	     idxtype *part, int *options, MPI_Comm comm);
void  PARKMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, idxtype *part, int *options, MPI_Comm comm);

//C     void __cdecl PARGKMETIS(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt,
//C                  int ndims, float *xyz, idxtype *part, int *options, MPI_Comm comm);
void  PARGKMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, int ndims, float *xyz, idxtype *part, int *options, MPI_Comm comm);

//C     void __cdecl PARGRMETIS(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt,
//C                  int ndims, float *xyz, idxtype *part, int *options, MPI_Comm comm);
void  PARGRMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, int ndims, float *xyz, idxtype *part, int *options, MPI_Comm comm);

//C     void __cdecl PARGMETIS(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int ndims, float *xyz,
//C                  idxtype *part, int *options, MPI_Comm comm);
void  PARGMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int ndims, float *xyz, idxtype *part, int *options, MPI_Comm comm);

//C     void __cdecl PARRMETIS(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, 
//C     	     idxtype *adjwgt, idxtype *part, int *options, MPI_Comm comm);
void  PARRMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, idxtype *part, int *options, MPI_Comm comm);

//C     void __cdecl PARUAMETIS(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, 
//C     	     idxtype *adjwgt, idxtype *part, int *options, MPI_Comm comm);
void  PARUAMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, idxtype *part, int *options, MPI_Comm comm);

//C     void __cdecl PARDAMETIS(
//C                  idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt,
//C                  idxtype *part, int *options, MPI_Comm comm);
void  PARDAMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, idxtype *part, int *options, MPI_Comm comm);

//C     #ifdef __cplusplus
//C     }
//C     #endif




//C     #endif 
