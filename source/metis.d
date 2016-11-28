/* Converted to D from metis.h by htod */
module metis;
/*!
\file metis.h 
\brief This file contains function prototypes and constant definitions for METIS
 *
\author George
\date   Started 8/9/02
\version\verbatim $Id$\endverbatim
*/

//C     #ifndef _METIS_H_
//C     #define _METIS_H_ 

/****************************************************************************
* A set of defines that can be modified by the user
*****************************************************************************/

/*--------------------------------------------------------------------------
 Specifies the width of the elementary data type that will hold information
 about vertices and their adjacency lists.

 Possible values:
   32 : Use 32 bit signed integers
   64 : Use 64 bit signed integers

 A width of 64 should be specified if the number of vertices or the total
 number of edges in the graph exceed the limits of a 32 bit signed integer
 i.e., 2^31-1.
 Proper use of 64 bit integers requires that the c99 standard datatypes
 int32_t and int64_t are supported by the compiler.
 GCC does provides these definitions in stdint.h, but it may require some
 modifications on other architectures.
--------------------------------------------------------------------------*/
//C     #define IDXTYPEWIDTH 32

const IDXTYPEWIDTH = 32;

/*--------------------------------------------------------------------------
 Specifies the data type that will hold floating-point style information.

 Possible values:
   32 : single precission floating point (float)
   64 : double precission floating point (double)
--------------------------------------------------------------------------*/
//C     #define REALTYPEWIDTH 32

const REALTYPEWIDTH = 32;


/****************************************************************************
* In principle, nothing needs to be changed beyond this point, unless the
* int32_t and int64_t cannot be found in the normal places.
*****************************************************************************/

/* Uniform definitions for various compilers */
//C     #if defined(_MSC_VER)
//C       #define COMPILER_MSC
//C     #endif
//C     #if defined(__ICC)
//C       #define COMPILER_ICC
//C     #endif
//C     #if defined(__GNUC__)
//C       #define COMPILER_GCC
//C     #endif

/* Include c99 int definitions and need constants. When building the library,
 * these are already defined by GKlib; hence the test for _GKLIB_H_ */
//C     #ifndef _GKLIB_H_
//C     #ifdef COMPILER_MSC
//#include <limits.h>

//C     typedef __int32 int32_t;
//C     typedef __int64 int64_t;
//C     #define PRId32       "I32d"
//C     #define PRId64       "I64d"
//C     #define SCNd32       "ld"
//C     #define SCNd64       "I64d"
//C     #define INT32_MIN    ((int32_t)_I32_MIN)
//C     #define INT32_MAX    _I32_MAX
//C     #define INT64_MIN    ((int64_t)_I64_MIN)
//C     #define INT64_MAX    _I64_MAX
//C     #else
//#include <inttypes.h>
//C     #endif
//C     #endif


/*------------------------------------------------------------------------
* Setup the basic datatypes
*-------------------------------------------------------------------------*/
//C     #if IDXTYPEWIDTH == 32

//C       typedef int idx_t;
extern (C):

alias int idx_t;

import core.stdc.stdint;
import core.stdc.math;
import core.stdc.stdlib;
import core.stdc.float_;
import core.stdc.inttypes;

//C       #define IDX_MAX   INT32_MAX
//C       #define IDX_MIN   INT32_MIN
alias INT32_MAX IDX_MAX;

alias INT32_MIN IDX_MIN;
//C       #define SCIDX  SCNd32
//C       #define PRIDX  PRId32
alias SCNd32 SCIDX;

alias PRId32 PRIDX;
//C       #define strtoidx      strtol
//C       #define iabs          abs
alias strtol strtoidx;
//C     #elif IDXTYPEWIDTH == 64
alias abs iabs;
//C       typedef int64_t idx_t;

//C       #define IDX_MAX   INT64_MAX
//C       #define IDX_MIN   INT64_MIN

//C       #define SCIDX  SCNd64
//C       #define PRIDX  PRId64

//C     #ifdef COMPILER_MSC
//C       #define strtoidx      _strtoi64
//C     #else
//C       #define strtoidx      strtoll
//C     #endif
//C       #define iabs          labs
//C     #else
//C       #error "Incorrect user-supplied value fo IDXTYPEWIDTH"
//C     #endif


//C     #if REALTYPEWIDTH == 32
//C       typedef float real_t;
alias float real_t;

//C       #define SCREAL         "f"
//C       #define PRREAL         "f"
//C       #define REAL_MAX       FLT_MAX
//C       #define REAL_MIN       FLT_MIN
alias FLT_MAX REAL_MAX;
//C       #define REAL_EPSILON   FLT_EPSILON
alias FLT_MIN REAL_MIN;

alias FLT_EPSILON REAL_EPSILON;
//C       #define rabs          fabsf
//C       #define REALEQ(x,y) ((rabs((x)-(y)) <= FLT_EPSILON))
alias fabsf rabs;

//C     #ifdef COMPILER_MSC
//C       #define strtoreal     (float)strtod
//C     #else
//C       #define strtoreal     strtof
//C     #endif
alias strtof strtoreal;
//C     #elif REALTYPEWIDTH == 64
//C       typedef double real_t;

//C       #define SCREAL         "lf"
//C       #define PRREAL         "lf"
//C       #define REAL_MAX       DBL_MAX
//C       #define REAL_MIN       DBL_MIN
//C       #define REAL_EPSILON   DBL_EPSILON

//C       #define rabs          fabs
//C       #define REALEQ(x,y) ((rabs((x)-(y)) <= DBL_EPSILON))

//C       #define strtoreal     strtod
//C     #else
//C       #error "Incorrect user-supplied value for REALTYPEWIDTH"
//C     #endif


/*------------------------------------------------------------------------
* Constant definitions 
*-------------------------------------------------------------------------*/
/* Metis's version number */
//C     #define METIS_VER_MAJOR         5
//C     #define METIS_VER_MINOR         1
const METIS_VER_MAJOR = 5;
//C     #define METIS_VER_SUBMINOR      0
const METIS_VER_MINOR = 1;

const METIS_VER_SUBMINOR = 0;
/* The maximum length of the options[] array */
//C     #define METIS_NOPTIONS          40

const METIS_NOPTIONS = 40;


/*------------------------------------------------------------------------
* Function prototypes 
*-------------------------------------------------------------------------*/

//C     #ifdef _WINDLL
//C     #define METIS_API(type) __declspec(dllexport) type __cdecl
//C     #elif defined(__cdecl)
//C     #define METIS_API(type) type __cdecl
//C     #else
//C     #define METIS_API(type) type
//C     #endif



//C     #ifdef __cplusplus
//C     extern "C" {
//C     #endif

//C     METIS_API(int) METIS_PartGraphRecursive(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, 
//C                       idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, 
//C                       idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
//C                       idx_t *edgecut, idx_t *part);
int  METIS_PartGraphRecursive(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part);

//C     METIS_API(int) METIS_PartGraphKway(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, 
//C                       idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, 
//C                       idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
//C                       idx_t *edgecut, idx_t *part);
int  METIS_PartGraphKway(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part);

//C     METIS_API(int) METIS_MeshToDual(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, 
//C                       idx_t *ncommon, idx_t *numflag, idx_t **r_xadj, idx_t **r_adjncy);
int  METIS_MeshToDual(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, idx_t *ncommon, idx_t *numflag, idx_t **r_xadj, idx_t **r_adjncy);

//C     METIS_API(int) METIS_MeshToNodal(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, 
//C                       idx_t *numflag, idx_t **r_xadj, idx_t **r_adjncy);
int  METIS_MeshToNodal(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, idx_t *numflag, idx_t **r_xadj, idx_t **r_adjncy);

//C     METIS_API(int) METIS_PartMeshNodal(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind,
//C                       idx_t *vwgt, idx_t *vsize, idx_t *nparts, real_t *tpwgts, 
//C                       idx_t *options, idx_t *objval, idx_t *epart, idx_t *npart);
int  METIS_PartMeshNodal(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, idx_t *vwgt, idx_t *vsize, idx_t *nparts, real_t *tpwgts, idx_t *options, idx_t *objval, idx_t *epart, idx_t *npart);

//C     METIS_API(int) METIS_PartMeshDual(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind,
//C                       idx_t *vwgt, idx_t *vsize, idx_t *ncommon, idx_t *nparts, 
//C                       real_t *tpwgts, idx_t *options, idx_t *objval, idx_t *epart, 
//C                       idx_t *npart);
int  METIS_PartMeshDual(idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, idx_t *vwgt, idx_t *vsize, idx_t *ncommon, idx_t *nparts, real_t *tpwgts, idx_t *options, idx_t *objval, idx_t *epart, idx_t *npart);

//C     METIS_API(int) METIS_NodeND(idx_t *nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
//C                       idx_t *options, idx_t *perm, idx_t *iperm);
int  METIS_NodeND(idx_t *nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *options, idx_t *perm, idx_t *iperm);

//C     METIS_API(int) METIS_Free(void *ptr);
int  METIS_Free(void *ptr);

//C     METIS_API(int) METIS_SetDefaultOptions(idx_t *options);
int  METIS_SetDefaultOptions(idx_t *options);


/* These functions are used by ParMETIS */

//C     METIS_API(int) METIS_NodeNDP(idx_t nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
//C                        idx_t npes, idx_t *options, idx_t *perm, idx_t *iperm, 
//C                        idx_t *sizes);
int  METIS_NodeNDP(idx_t nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t npes, idx_t *options, idx_t *perm, idx_t *iperm, idx_t *sizes);

//C     METIS_API(int) METIS_ComputeVertexSeparator(idx_t *nvtxs, idx_t *xadj, idx_t *adjncy, 
//C                        idx_t *vwgt, idx_t *options, idx_t *sepsize, idx_t *part);
int  METIS_ComputeVertexSeparator(idx_t *nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *options, idx_t *sepsize, idx_t *part);

//C     METIS_API(int) METIS_NodeRefine(idx_t nvtxs, idx_t *xadj, idx_t *vwgt, idx_t *adjncy,
//C                        idx_t *where, idx_t *hmarker, real_t ubfactor);
int  METIS_NodeRefine(idx_t nvtxs, idx_t *xadj, idx_t *vwgt, idx_t *adjncy, idx_t *where, idx_t *hmarker, real_t ubfactor);


//C     #ifdef __cplusplus
//C     }
//C     #endif



/*------------------------------------------------------------------------
* Enum type definitions 
*-------------------------------------------------------------------------*/
/*! Return codes */
//C     typedef enum {
//C       METIS_OK              = 1,    /*!< Returned normally */
//C       METIS_ERROR_INPUT     = -2,   /*!< Returned due to erroneous inputs and/or options */
//C       METIS_ERROR_MEMORY    = -3,   /*!< Returned due to insufficient memory */
//C       METIS_ERROR           = -4    /*!< Some other errors */
//C     } rstatus_et; 
enum
{
    METIS_OK = 1,
    METIS_ERROR_INPUT = -2,
    METIS_ERROR_MEMORY = -3,
    METIS_ERROR = -4,
}
alias int rstatus_et;


/*! Operation type codes */
//C     typedef enum {
//C       METIS_OP_PMETIS,       
//C       METIS_OP_KMETIS,
//C       METIS_OP_OMETIS
//C     } moptype_et;
enum
{
    METIS_OP_PMETIS,
    METIS_OP_KMETIS,
    METIS_OP_OMETIS,
}
alias int moptype_et;


/*! Options codes (i.e., options[]) */
//C     typedef enum {
//C       METIS_OPTION_PTYPE,
//C       METIS_OPTION_OBJTYPE,
//C       METIS_OPTION_CTYPE,
//C       METIS_OPTION_IPTYPE,
//C       METIS_OPTION_RTYPE,
//C       METIS_OPTION_DBGLVL,
//C       METIS_OPTION_NITER,
//C       METIS_OPTION_NCUTS,
//C       METIS_OPTION_SEED,
//C       METIS_OPTION_NO2HOP,
//C       METIS_OPTION_MINCONN,
//C       METIS_OPTION_CONTIG,
//C       METIS_OPTION_COMPRESS,
//C       METIS_OPTION_CCORDER,
//C       METIS_OPTION_PFACTOR,
//C       METIS_OPTION_NSEPS,
//C       METIS_OPTION_UFACTOR,
//C       METIS_OPTION_NUMBERING,

  /* Used for command-line parameter purposes */
//C       METIS_OPTION_HELP,
//C       METIS_OPTION_TPWGTS,
//C       METIS_OPTION_NCOMMON,
//C       METIS_OPTION_NOOUTPUT,
//C       METIS_OPTION_BALANCE,
//C       METIS_OPTION_GTYPE,
//C       METIS_OPTION_UBVEC
//C     } moptions_et;
enum
{
    METIS_OPTION_PTYPE,
    METIS_OPTION_OBJTYPE,
    METIS_OPTION_CTYPE,
    METIS_OPTION_IPTYPE,
    METIS_OPTION_RTYPE,
    METIS_OPTION_DBGLVL,
    METIS_OPTION_NITER,
    METIS_OPTION_NCUTS,
    METIS_OPTION_SEED,
    METIS_OPTION_NO2HOP,
    METIS_OPTION_MINCONN,
    METIS_OPTION_CONTIG,
    METIS_OPTION_COMPRESS,
    METIS_OPTION_CCORDER,
    METIS_OPTION_PFACTOR,
    METIS_OPTION_NSEPS,
    METIS_OPTION_UFACTOR,
    METIS_OPTION_NUMBERING,
    METIS_OPTION_HELP,
    METIS_OPTION_TPWGTS,
    METIS_OPTION_NCOMMON,
    METIS_OPTION_NOOUTPUT,
    METIS_OPTION_BALANCE,
    METIS_OPTION_GTYPE,
    METIS_OPTION_UBVEC,
}
alias int moptions_et;


/*! Partitioning Schemes */
//C     typedef enum {
//C       METIS_PTYPE_RB, 
//C       METIS_PTYPE_KWAY                
//C     } mptype_et;
enum
{
    METIS_PTYPE_RB,
    METIS_PTYPE_KWAY,
}
alias int mptype_et;

/*! Graph types for meshes */
//C     typedef enum {
//C       METIS_GTYPE_DUAL,
//C       METIS_GTYPE_NODAL               
//C     } mgtype_et;
enum
{
    METIS_GTYPE_DUAL,
    METIS_GTYPE_NODAL,
}
alias int mgtype_et;

/*! Coarsening Schemes */
//C     typedef enum {
//C       METIS_CTYPE_RM,
//C       METIS_CTYPE_SHEM
//C     } mctype_et;
enum
{
    METIS_CTYPE_RM,
    METIS_CTYPE_SHEM,
}
alias int mctype_et;

/*! Initial partitioning schemes */
//C     typedef enum {
//C       METIS_IPTYPE_GROW,
//C       METIS_IPTYPE_RANDOM,
//C       METIS_IPTYPE_EDGE,
//C       METIS_IPTYPE_NODE,
//C       METIS_IPTYPE_METISRB
//C     } miptype_et;
enum
{
    METIS_IPTYPE_GROW,
    METIS_IPTYPE_RANDOM,
    METIS_IPTYPE_EDGE,
    METIS_IPTYPE_NODE,
    METIS_IPTYPE_METISRB,
}
alias int miptype_et;


/*! Refinement schemes */
//C     typedef enum {
//C       METIS_RTYPE_FM,
//C       METIS_RTYPE_GREEDY,
//C       METIS_RTYPE_SEP2SIDED,
//C       METIS_RTYPE_SEP1SIDED
//C     } mrtype_et;
enum
{
    METIS_RTYPE_FM,
    METIS_RTYPE_GREEDY,
    METIS_RTYPE_SEP2SIDED,
    METIS_RTYPE_SEP1SIDED,
}
alias int mrtype_et;


/*! Debug Levels */
//C     typedef enum {
//C       METIS_DBG_INFO       = 1,       /*!< Shows various diagnostic messages */
//C       METIS_DBG_TIME       = 2,       /*!< Perform timing analysis */
//C       METIS_DBG_COARSEN    = 4,	  /*!< Show the coarsening progress */
//C       METIS_DBG_REFINE     = 8,	  /*!< Show the refinement progress */
//C       METIS_DBG_IPART      = 16, 	  /*!< Show info on initial partitioning */
//C       METIS_DBG_MOVEINFO   = 32, 	  /*!< Show info on vertex moves during refinement */
//C       METIS_DBG_SEPINFO    = 64, 	  /*!< Show info on vertex moves during sep refinement */
//C       METIS_DBG_CONNINFO   = 128,     /*!< Show info on minimization of subdomain connectivity */
//C       METIS_DBG_CONTIGINFO = 256,     /*!< Show info on elimination of connected components */ 
//C       METIS_DBG_MEMORY     = 2048,    /*!< Show info related to wspace allocation */
//C     } mdbglvl_et;
enum
{
    METIS_DBG_INFO = 1,
    METIS_DBG_TIME,
    METIS_DBG_COARSEN = 4,
    METIS_DBG_REFINE = 8,
    METIS_DBG_IPART = 16,
    METIS_DBG_MOVEINFO = 32,
    METIS_DBG_SEPINFO = 64,
    METIS_DBG_CONNINFO = 128,
    METIS_DBG_CONTIGINFO = 256,
    METIS_DBG_MEMORY = 2048,
}
alias int mdbglvl_et;


/* Types of objectives */
//C     typedef enum {
//C       METIS_OBJTYPE_CUT,
//C       METIS_OBJTYPE_VOL,
//C       METIS_OBJTYPE_NODE
//C     } mobjtype_et;
enum
{
    METIS_OBJTYPE_CUT,
    METIS_OBJTYPE_VOL,
    METIS_OBJTYPE_NODE,
}
alias int mobjtype_et;



//C     #endif  /* _METIS_H_ */
