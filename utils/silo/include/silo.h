/*
Copyright (c) 1994 - 2010, Lawrence Livermore National Security, LLC.
LLNL-CODE-425250.
All rights reserved.

This file is part of Silo. For details, see silo.llnl.gov.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the disclaimer below.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the disclaimer (as noted
     below) in the documentation and/or other materials provided with
     the distribution.
   * Neither the name of the LLNS/LLNL nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

THIS SOFTWARE  IS PROVIDED BY  THE COPYRIGHT HOLDERS  AND CONTRIBUTORS
"AS  IS" AND  ANY EXPRESS  OR IMPLIED  WARRANTIES, INCLUDING,  BUT NOT
LIMITED TO, THE IMPLIED  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A  PARTICULAR  PURPOSE ARE  DISCLAIMED.  IN  NO  EVENT SHALL  LAWRENCE
LIVERMORE  NATIONAL SECURITY, LLC,  THE U.S.  DEPARTMENT OF  ENERGY OR
CONTRIBUTORS BE LIABLE FOR  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR  CONSEQUENTIAL DAMAGES  (INCLUDING, BUT NOT  LIMITED TO,
PROCUREMENT OF  SUBSTITUTE GOODS  OR SERVICES; LOSS  OF USE,  DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER  IN CONTRACT, STRICT LIABILITY,  OR TORT (INCLUDING
NEGLIGENCE OR  OTHERWISE) ARISING IN  ANY WAY OUT  OF THE USE  OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

This work was produced at Lawrence Livermore National Laboratory under
Contract No.  DE-AC52-07NA27344 with the DOE.

Neither the  United States Government nor  Lawrence Livermore National
Security, LLC nor any of  their employees, makes any warranty, express
or  implied,  or  assumes  any  liability or  responsibility  for  the
accuracy, completeness,  or usefulness of  any information, apparatus,
product, or  process disclosed, or  represents that its use  would not
infringe privately-owned rights.

Any reference herein to  any specific commercial products, process, or
services by trade name,  trademark, manufacturer or otherwise does not
necessarily  constitute or imply  its endorsement,  recommendation, or
favoring  by  the  United  States  Government  or  Lawrence  Livermore
National Security,  LLC. The views  and opinions of  authors expressed
herein do not necessarily state  or reflect those of the United States
Government or Lawrence Livermore National Security, LLC, and shall not
be used for advertising or product endorsement purposes.
*/

/*
 * SILO Public header file.
 *
 * This header file defines public constants and public prototypes.
 * Before including this file, the application should define
 * which file formats will be used.
 *
 * WARNING: The `#define' statements in this file are used when
 *      generating the Fortran include file `silo.inc'.  Any
 *     such symbol that should not be an integer parameter
 *     in the Fortran include file should have the text 
 *     `NO_FORTRAN_DEFINE' on the same line.  #define statements
 *     that define macros (or any value not beginning with
 *     one of [a-zA-Z0-9_]) are ignored.
 */
#ifndef SILO_H
#define SILO_H

#ifdef __cplusplus
extern "C" {
#endif

/* Set the base type for datatype'd pointers (that is pointers whose
   ultimate type is deteremined by an additional 'int datatype' function
   argument or struct member) as float (legacy) and void (modern). The
   DB_DTPTR is the base type. The '1' and '2' variants are for singley
   subscripted and doubley subscripted arrays, respectively. If the
   definitions of DB_DTPTR below reference 'float', then this silo.h
   header file was configured with --enable-legacy-datatyped-pointers
   and it represents the legacy (float) pointers that the silo
   library has always had since its original writing. If, instead,
   you see 'void' (the default configuration), then this silo.h header
   file is using the modern (void) pointers. In that case, note also
   that because C compiler's often do not handle correctly nor
   distinguish between a void* and a void**, both the singley and
   doubley subscripted variants will have only a single star. Rest
   assured they are still treated as doubley subscripted in the
   implementation. */
#define DB_DTPTR  void  /* NO_FORTRAN_DEFINE */
#define DB_DTPTR1 void* /* NO_FORTRAN_DEFINE */
#define DB_DTPTR2 void* /* NO_FORTRAN_DEFINE */

/* Permit client to explicitly require the legacy mode
   for datatyped pointers */
#ifdef DB_USE_LEGACY_DTPTR
#ifdef DB_USE_MODERN_DTPTR
#error cannot specify BOTH legacy and modern datatyped pointers 
#endif
#undef DB_DTPTR  /* NO_FORTRAN_DEFINE */
#undef DB_DTPTR1 /* NO_FORTRAN_DEFINE */
#undef DB_DTPTR2 /* NO_FORTRAN_DEFINE */
#define DB_DTPTR  float   /* NO_FORTRAN_DEFINE */
#define DB_DTPTR1 float*  /* NO_FORTRAN_DEFINE */
#define DB_DTPTR2 float** /* NO_FORTRAN_DEFINE */
#endif

/* Permit client to explicitly require the modern mode
   for datatyped pointers */
#ifdef DB_USE_MODERN_DTPTR
#undef DB_DTPTR  /* NO_FORTRAN_DEFINE */
#undef DB_DTPTR1 /* NO_FORTRAN_DEFINE */
#undef DB_DTPTR2 /* NO_FORTRAN_DEFINE */
#define DB_DTPTR  void  /* NO_FORTRAN_DEFINE */
#define DB_DTPTR1 void* /* NO_FORTRAN_DEFINE */
#define DB_DTPTR2 void* /* NO_FORTRAN_DEFINE */
#endif

#include <stdio.h>

#ifndef FALSE
#define FALSE   0
#endif
#ifndef TRUE
#define TRUE    1
#endif

/* In the definitions for different parts of the version number, below,
   we use leading '0x0' to deal with possible blank minor and/or patch
   version number but still allow base-10, numeric comparison in the GE
   macro. */

/* Major release number of silo library. */
#define SILO_VERS_MAJ 4

/* Minor release number of silo library. Can be empty. */
#define SILO_VERS_MIN 0x08

/* Patch release number of silo library.  Can be empty. */
#define SILO_VERS_PAT 0x0

/* Pre-release release number of silo library.  Can be empty. */
#define SILO_VERS_PRE 

/* The symbol Silo uses to enforce link-time
   header/object version compatibility */
#define SILO_VERS_TAG Silo_version_4_8

/* Useful macro for comparing Silo versions (and DB_ alias) */
#define SILO_VERSION_GE(Maj,Min,Pat)  \
        (((SILO_VERS_MAJ==Maj) && (SILO_VERS_MIN==0x0 ## Min) && (SILO_VERS_PAT>=0x0 ## Pat)) || \
         ((SILO_VERS_MAJ==Maj) && (SILO_VERS_MIN>0x0 ## Min)) || \
         (SILO_VERS_MAJ>Maj))
#define DB_VERSION_GE(Maj,Min,Pat) SILO_VERSION_GE(Maj,Min,Pat)

/*-------------------------------------------------------------------------
 * Drivers.  This is a list of every driver that a user could use.  Not all of
 * them are necessarily compiled into the library.  However, users are free
 * to try without getting compilation errors.  They are listed here so that
 * silo.h doesn't have to be generated every time the library is recompiled.
 *--------------------------------------------------------------------------*/
#define DB_NETCDF 0
#define DB_PDB 2 /* PDB Lite */
#define DB_TAURUS 3
#define DB_UNKNOWN 5
#define DB_DEBUG 6
#define DB_HDF5X 7
#define DB_PDBP 1 /* PDB Proper */

/* DO NOT USE. Obsoleted ways of specifying different HDF5 vfds */
#define DB_HDF5_SEC2_OBSOLETE 0x100
#define DB_HDF5_STDIO_OBSOLETE 0x200
#define DB_HDF5_CORE_OBSOLETE 0x300
#define DB_HDF5_MPIO_OBSOLETE 0x400
#define DB_HDF5_MPIOP_OBSOLETE 0x500

/* symbols for various HDF5 vfds */
#define DB_H5VFD_DEFAULT 0
#define DB_H5VFD_SEC2    1
#define DB_H5VFD_STDIO   2
#define DB_H5VFD_CORE    3
#define DB_H5VFD_LOG     4
#define DB_H5VFD_SPLIT   5
#define DB_H5VFD_DIRECT  6
#define DB_H5VFD_FAMILY  7
#define DB_H5VFD_MPIO    8
#define DB_H5VFD_MPIP    9
#define DB_H5VFD_SILO    10

/* Macro for defining various HDF5 vfds as 'type' arg in create/open.
   The 11 bit shift is to avoid possible collision with older versions
   of Silo header file where VFDs where specified in bits 8-11. Their
   obsoleted values are listed above. */ 
#define DB_HDF5_OPTS(OptsId) (DB_HDF5X|((OptsId&0x3F)<<11))

/* Monikers for default file options sets */
/* We just make the default options sets the same as the vfd is */
#define DB_FILE_OPTS_H5_DEFAULT_DEFAULT DB_H5VFD_DEFAULT 
#define DB_FILE_OPTS_H5_DEFAULT_SEC2    DB_H5VFD_SEC2 
#define DB_FILE_OPTS_H5_DEFAULT_STDIO   DB_H5VFD_STDIO 
#define DB_FILE_OPTS_H5_DEFAULT_CORE    DB_H5VFD_CORE 
#define DB_FILE_OPTS_H5_DEFAULT_LOG     DB_H5VFD_LOG 
#define DB_FILE_OPTS_H5_DEFAULT_SPLIT   DB_H5VFD_SPLIT 
#define DB_FILE_OPTS_H5_DEFAULT_DIRECT  DB_H5VFD_DIRECT 
#define DB_FILE_OPTS_H5_DEFAULT_FAMILY  DB_H5VFD_FAMILY 
#define DB_FILE_OPTS_H5_DEFAULT_MPIO    DB_H5VFD_MPIO
#define DB_FILE_OPTS_H5_DEFAULT_MPIP    DB_H5VFD_MPIP
#define DB_FILE_OPTS_H5_DEFAULT_SILO    DB_H5VFD_SILO 
#define DB_FILE_OPTS_LAST DB_FILE_OPTS_H5_DEFAULT_SILO

/* Various default HDF5 driver options. Users can define their own
   sets of options using DBRegisterFileOptionsSets(). */
#define DB_HDF5 DB_HDF5_OPTS(DB_FILE_OPTS_H5_DEFAULT_DEFAULT)
#define DB_HDF5_SEC2 DB_HDF5_OPTS(DB_FILE_OPTS_H5_DEFAULT_SEC2)
#define DB_HDF5_STDIO DB_HDF5_OPTS(DB_FILE_OPTS_H5_DEFAULT_STDIO)
#define DB_HDF5_CORE DB_HDF5_OPTS(DB_FILE_OPTS_H5_DEFAULT_CORE)
#define DB_HDF5_LOG DB_HDF5_OPTS(DB_FILE_OPTS_H5_DEFAULT_LOG)
#define DB_HDF5_SPLIT DB_HDF5_OPTS(DB_FILE_OPTS_H5_DEFAULT_SPLIT)
#define DB_HDF5_DIRECT DB_HDF5_OPTS(DB_FILE_OPTS_H5_DEFAULT_DIRECT)
#define DB_HDF5_FAMILY DB_HDF5_OPTS(DB_FILE_OPTS_H5_DEFAULT_FAMILY)
#define DB_HDF5_MPIO DB_HDF5_OPTS(DB_FILE_OPTS_H5_DEFAULT_MPIO)
#define DB_HDF5_MPIOP DB_HDF5_OPTS(DB_FILE_OPTS_H5_DEFAULT_MPIP)
#define DB_HDF5_MPIP DB_HDF5_OPTS(DB_FILE_OPTS_H5_DEFAULT_MPIP)
#define DB_HDF5_SILO DB_HDF5_OPTS(DB_FILE_OPTS_H5_DEFAULT_SILO)

/*-------------------------------------------------------------------------
 * Other library-wide constants.
 *-------------------------------------------------------------------------*/
#define DB_NFILES       256         /*Max simultaneously open files */
#define DB_NFILTERS     32          /*Number of filters defined */

/*-------------------------------------------------------------------------
 * Constants.  All of these constants are always defined in the application.
 * Each group of constants defined here are small integers used as an index
 * into an array.  Many of the groups have a total count of items in the
 * group that will be used for array allocation and error checking--don't
 * forget to increment the value when adding a new item to a constant group.
 *-------------------------------------------------------------------------
 */

/* The following identifiers are for use with the DBDataReadMask() call.  They
 * specify what portions of the data beyond the metadata is allocated
 * and read.  Note that since these values are only necessary when reading
 * and since the Fortran interface is primarily a write interface, it is not
 * necessary for these symbols to appear in the silo.inc file. However, the
 * reason they DO NOT APPEAR there inspite of the fact that DO NOT HAVE the
 * 'NO_FORTRAN_DEFINE' text appearing on each line is that the mkinc script
 * requires an underscore in the symbol name for it to appear in silo.inc. */
#define DBAll                0xffffffff
#define DBNone               0x00000000
#define DBCalc               0x00000001
#define DBMatMatnos          0x00000002
#define DBMatMatlist         0x00000004
#define DBMatMixList         0x00000008
#define DBCurveArrays        0x00000010
#define DBPMCoords           0x00000020
#define DBPVData             0x00000040
#define DBQMCoords           0x00000080
#define DBQVData             0x00000100
#define DBUMCoords           0x00000200
#define DBUMFacelist         0x00000400
#define DBUMZonelist         0x00000800
#define DBUVData             0x00001000
#define DBFacelistInfo       0x00002000
#define DBZonelistInfo       0x00004000
#define DBMatMatnames        0x00008000
#define DBUMGlobNodeNo       0x00010000
#define DBZonelistGlobZoneNo 0x00020000
#define DBMatMatcolors       0x00040000
#define DBCSGMBoundaryInfo   0x00080000
#define DBCSGMZonelist       0x00100000
#define DBCSGMBoundaryNames  0x00200000
#define DBCSGVData           0x00400000
#define DBCSGZonelistZoneNames 0x00800000
#define DBCSGZonelistRegNames  0x01000000
#define DBMMADJNodelists     0x02000000
#define DBMMADJZonelists     0x04000000
#define DBPMGlobNodeNo       0x08000000

/* Macros used for exporting symbols on Win32 systems. */
#ifndef SILO_API
#ifdef _WIN32
/* Make Silo a DLL by default. */
#ifdef SILO_STATIC_LIBRARY
#define SILO_API
#else
#ifdef SILO_EXPORTS
#define SILO_API __declspec(dllexport)
#else
#define SILO_API __declspec(dllimport)
#endif
#endif
#else
#define SILO_API
#endif
#endif

/* Definitions for COORD_TYPE */
/* Placed before DBObjectType enum because the
   DB_QUAD_CURV and DB_QUAD_RECT symbols are
   sometimes used as DBObjectType */

#define  DB_COLLINEAR           130
#define  DB_NONCOLLINEAR        131
#define  DB_QUAD_RECT           DB_COLLINEAR
#define  DB_QUAD_CURV           DB_NONCOLLINEAR

/* Objects that can be stored in a data file */
typedef enum {
    DB_INVALID_OBJECT= -1,       /*causes enum to be signed, do not remove,
                                   space before minus sign necessary for lint*/
    DB_QUADRECT = DB_QUAD_RECT,
    DB_QUADCURV = DB_QUAD_CURV,
    DB_QUADMESH=500,
    DB_QUADVAR=501,
    DB_UCDMESH=510,
    DB_UCDVAR=511,
    DB_MULTIMESH=520,
    DB_MULTIVAR=521,
    DB_MULTIMAT=522,
    DB_MULTIMATSPECIES=523,
    DB_MULTIBLOCKMESH=DB_MULTIMESH,
    DB_MULTIBLOCKVAR=DB_MULTIVAR,
    DB_MULTIMESHADJ=524,
    DB_MATERIAL=530,
    DB_MATSPECIES=531,
    DB_FACELIST=550,
    DB_ZONELIST=551,
    DB_EDGELIST=552,
    DB_PHZONELIST=553,
    DB_CSGZONELIST=554,
    DB_CSGMESH=555,
    DB_CSGVAR=556,
    DB_CURVE=560,
    DB_DEFVARS=565,
    DB_POINTMESH=570,
    DB_POINTVAR=571,
    DB_ARRAY=580,
    DB_DIR=600,
    DB_VARIABLE=610,
    DB_MRGTREE=611,
    DB_GROUPELMAP=612,
    DB_MRGVAR=613,
    DB_USERDEF=700
} DBObjectType;

/* Data types */
typedef enum {
    DB_INT=16,
    DB_SHORT=17,
    DB_LONG=18,
    DB_FLOAT=19,
    DB_DOUBLE=20,
    DB_CHAR=21,
    DB_LONG_LONG=22,
    DB_NOTYPE=25           /*unknown type */
} DBdatatype;

/* Flags for DBCreate */
#define         DB_CLOBBER      0
#define         DB_NOCLOBBER    1

/* Flags for DBOpen */
#define         DB_READ         1
#define         DB_APPEND       2

/* Target machine for DBCreate */
#define         DB_LOCAL        0
#define         DB_SUN3         10
#define         DB_SUN4         11
#define         DB_SGI          12
#define         DB_RS6000       13
#define         DB_CRAY         14
#define         DB_INTEL        15

/* Options */
#define DBOPT_FIRST             260
#define DBOPT_ALIGN             260
#define DBOPT_COORDSYS          262
#define DBOPT_CYCLE             263
#define DBOPT_FACETYPE          264
#define DBOPT_HI_OFFSET         265
#define DBOPT_LO_OFFSET         266
#define DBOPT_LABEL             267
#define DBOPT_XLABEL            268
#define DBOPT_YLABEL            269
#define DBOPT_ZLABEL            270
#define DBOPT_MAJORORDER        271
#define DBOPT_NSPACE            272
#define DBOPT_ORIGIN            273
#define DBOPT_PLANAR            274
#define DBOPT_TIME              275
#define DBOPT_UNITS             276
#define DBOPT_XUNITS            277
#define DBOPT_YUNITS            278
#define DBOPT_ZUNITS            279
#define DBOPT_DTIME             280
#define DBOPT_USESPECMF         281
#define DBOPT_XVARNAME          282
#define DBOPT_YVARNAME          283
#define DBOPT_ZVARNAME          284
#define DBOPT_ASCII_LABEL       285
#define DBOPT_MATNOS            286
#define DBOPT_NMATNOS           287
#define DBOPT_MATNAME           288
#define DBOPT_NMAT              289
#define DBOPT_NMATSPEC          290
#define DBOPT_BASEINDEX         291 /* quad meshes for node and zone */
#define DBOPT_ZONENUM           292 /* ucd meshes for zone */
#define DBOPT_NODENUM           293 /* ucd/point meshes for node */
#define DBOPT_BLOCKORIGIN       294
#define DBOPT_GROUPNUM          295
#define DBOPT_GROUPORIGIN       296
#define DBOPT_NGROUPS           297
#define DBOPT_MATNAMES          298
#define DBOPT_EXTENTS_SIZE      299
#define DBOPT_EXTENTS           300
#define DBOPT_MATCOUNTS         301
#define DBOPT_MATLISTS          302
#define DBOPT_MIXLENS           303
#define DBOPT_ZONECOUNTS        304
#define DBOPT_HAS_EXTERNAL_ZONES 305
#define DBOPT_PHZONELIST        306
#define DBOPT_MATCOLORS         307
#define DBOPT_BNDNAMES          308
#define DBOPT_REGNAMES          309
#define DBOPT_ZONENAMES         310
#define DBOPT_HIDE_FROM_GUI     311
#define DBOPT_TOPO_DIM          312 /* TOPOlogical DIMension */
#define DBOPT_REFERENCE         313 /* reference object */
#define DBOPT_GROUPINGS_SIZE    314 /* size of grouping array */
#define DBOPT_GROUPINGS         315 /* groupings array */
#define DBOPT_GROUPINGNAMES     316 /* array of size determined by 
                                       number of groups of names of groups. */
#define DBOPT_ALLOWMAT0         317 /* Turn off material numer "0" warnings*/
#define DBOPT_MRGTREE_NAME      318
#define DBOPT_REGION_PNAMES     319
#define DBOPT_TENSOR_RANK       320
#define DBOPT_MMESH_NAME        321
#define DBOPT_TV_CONNECTIVITY   322
#define DBOPT_DISJOINT_MODE     323
#define DBOPT_MRGV_ONAMES       324
#define DBOPT_MRGV_RNAMES       325
#define DBOPT_SPECNAMES         326
#define DBOPT_SPECCOLORS        327
#define DBOPT_LLONGNZNUM        328
#define DBOPT_CONSERVED         329
#define DBOPT_EXTENSIVE         330
#define DBOPT_MB_FILE_NS        331
#define DBOPT_MB_BLOCK_NS       332
#define DBOPT_MB_BLOCK_TYPE     333
#define DBOPT_MB_EMPTY_LIST     334
#define DBOPT_MB_EMPTY_COUNT    335
#define DBOPT_LAST              499 

/* Options relating to virtual file drivers */
#define DBOPT_H5_FIRST              500
#define DBOPT_H5_VFD                500
#define DBOPT_H5_RAW_FILE_OPTS      501
#define DBOPT_H5_RAW_EXTENSION      502
#define DBOPT_H5_META_FILE_OPTS     503
#define DBOPT_H5_META_EXTENSION     504
#define DBOPT_H5_CORE_ALLOC_INC     505
#define DBOPT_H5_CORE_NO_BACK_STORE 506
#define DBOPT_H5_META_BLOCK_SIZE    507
#define DBOPT_H5_SMALL_RAW_SIZE     508
#define DBOPT_H5_ALIGN_MIN          509
#define DBOPT_H5_ALIGN_VAL          510
#define DBOPT_H5_DIRECT_MEM_ALIGN   511
#define DBOPT_H5_DIRECT_BLOCK_SIZE  512
#define DBOPT_H5_DIRECT_BUF_SIZE    513
#define DBOPT_H5_LOG_NAME           514
#define DBOPT_H5_LOG_BUF_SIZE       515
#define DBOPT_H5_MPIO_COMM          516
#define DBOPT_H5_MPIO_INFO          517
#define DBOPT_H5_MPIP_NO_GPFS_HINTS 518
#define DBOPT_H5_SIEVE_BUF_SIZE     519
#define DBOPT_H5_CACHE_NELMTS       520
#define DBOPT_H5_CACHE_NBYTES       521
#define DBOPT_H5_CACHE_POLICY       522
#define DBOPT_H5_FAM_SIZE           523
#define DBOPT_H5_FAM_FILE_OPTS      524
#define DBOPT_H5_USER_DRIVER_ID     525
#define DBOPT_H5_USER_DRIVER_INFO   526
#define DBOPT_H5_SILO_BLOCK_SIZE    527
#define DBOPT_H5_SILO_BLOCK_COUNT   528
#define DBOPT_H5_SILO_LOG_STATS     529
#define DBOPT_H5_SILO_USE_DIRECT    530
#define DBOPT_H5_LAST               599

/* Error trapping method */
#define         DB_TOP          0 /*default--API traps  */
#define         DB_NONE         1 /*no errors trapped  */
#define         DB_ALL          2 /*all levels trap (traceback) */
#define         DB_ABORT        3 /*abort() is called  */
#define         DB_SUSPEND      4 /*suspend error reporting temporarily */
#define         DB_RESUME       5 /*resume normal error reporting */
#define         DB_ALL_AND_DRVR 6 /*DB_ALL + driver error reporting */

/* Errors */
#define     E_NOERROR   0       /*No error   */
#define     E_BADFTYPE  1       /*Bad file type   */
#define     E_NOTIMP    2       /*Callback not implemented */
#define     E_NOFILE    3       /*No data file specified    */
#define     E_INTERNAL  5       /*Internal error        */
#define     E_NOMEM     6       /*Not enough memory     */
#define     E_BADARGS   7       /*Bad argument to function  */
#define     E_CALLFAIL  8       /*Low-level function failure    */
#define     E_NOTFOUND  9       /*Object not found      */
#define     E_TAURSTATE 10      /*Taurus: database state error  */
#define     E_MSERVER   11      /*SDX: too many connections */
#define     E_PROTO     12      /*SDX: protocol error       */
#define     E_NOTDIR    13      /*Not a directory       */
#define     E_MAXOPEN   14      /*Too many open files  */
#define     E_NOTFILTER 15      /*Filter(s) not found  */
#define     E_MAXFILTERS    16  /*Too many filters  */
#define     E_FEXIST    17      /*File already exists  */
#define     E_FILEISDIR 18      /*File is actually a directory */
#define     E_FILENOREAD    19  /*File lacks read permission. */
#define     E_SYSTEMERR 20      /*System level error occured. */
#define     E_FILENOWRITE 21    /*File lacks write permission. */
#define     E_INVALIDNAME 22    /* Variable name is invalid */
#define     E_NOOVERWRITE 23    /*Overwrite not permitted */
#define     E_CHECKSUM  24      /*Checksum failed */
#define     E_COMPRESSION  25   /*Compression failed */
#define     E_GRABBED   26      /*Low level driver enabled */
#define     E_NOTREG    27      /*The dbfile pointer is not resitered. */
#define     E_CONCURRENT 28     /*File is opened multiply and not all read-only. */
#define     E_DRVRCANTOPEN 29   /*Driver cannot open the file. */
#define     E_BADOPTCLASS 30    /*Optlist contains options for wrong class */
#define     E_NOTENABLEDINBUILD 31 /*feature not enabled in this build */
#define     E_MAXFILEOPTSETS 32    /*Too many file options sets */
#define     E_NOHDF5 33         /*HDF5 driver not available */
#define     E_NERRORS   50

/* Definitions for MAJOR_ORDER */
#define  DB_ROWMAJOR            0
#define  DB_COLMAJOR            1

/* Definitions for CENTERING */
#define  DB_NOTCENT             0
#define  DB_NODECENT            110
#define  DB_ZONECENT            111
#define  DB_FACECENT            112
#define  DB_BNDCENT             113 /* for CSG meshes only */
#define  DB_EDGECENT            114
#define  DB_BLOCKCENT           115 /* for 'block-centered' data on multimeshs */

/* Definitions for COORD_SYSTEM */
#define  DB_CARTESIAN           120
#define  DB_CYLINDRICAL         121
#define  DB_SPHERICAL           122
#define  DB_NUMERICAL           123
#define  DB_OTHER               124

/* Definitions for ZONE FACE_TYPE */
#define  DB_RECTILINEAR         100
#define  DB_CURVILINEAR         101

/* Definitions for PLANAR */
#define  DB_AREA                140
#define  DB_VOLUME              141

/* Definitions for flag values */
#define DB_ON                    1000
#define DB_OFF                  -1000

/* Definitions for disjoint flag */
#define DB_ABUTTING              142
#define DB_FLOATING              143

/* Definitions for derived variable types */
#define DB_VARTYPE_SCALAR               200
#define DB_VARTYPE_VECTOR               201
#define DB_VARTYPE_TENSOR               202
#define DB_VARTYPE_SYMTENSOR            203
#define DB_VARTYPE_ARRAY                204
#define DB_VARTYPE_MATERIAL             205
#define DB_VARTYPE_SPECIES              206
#define DB_VARTYPE_LABEL                207

/* Definitions for CSG boundary types 
   Designed so low-order 16 bits are unused.

   The last few characters of the symbol are intended
   to indicate the representational form of the surface type

   G = generalized form  (n values, depends on surface type)
   P = point (3 values, x,y,z in 3D, 2 values in 2D x,y)
   N = normal (3 values, Nx,Ny,Nz in 3D, 2 values in 2D Nx,Ny)
   R = radius (1 value)
   A = angle (1 value in degrees)
   L = length (1 value)
   X = x-intercept (1 value)
   Y = y-intercept (1 value)
   Z = z-intercept (1 value)
   K = arbitrary integer
   F = planar face defined by point-normal pair (6 values)
   */
#define DBCSG_QUADRIC_G         0x01000000
#define DBCSG_SPHERE_PR         0x02010000
#define DBCSG_ELLIPSOID_PRRR    0x02020000
#define DBCSG_PLANE_G           0x03000000
#define DBCSG_PLANE_X           0x03010000
#define DBCSG_PLANE_Y           0x03020000
#define DBCSG_PLANE_Z           0x03030000
#define DBCSG_PLANE_PN          0x03040000
#define DBCSG_PLANE_PPP         0x03050000
#define DBCSG_CYLINDER_PNLR     0x04000000
#define DBCSG_CYLINDER_PPR      0x04010000
#define DBCSG_BOX_XYZXYZ        0x05000000
#define DBCSG_CONE_PNLA         0x06000000
#define DBCSG_CONE_PPA          0x06010000
#define DBCSG_POLYHEDRON_KF     0x07000000
#define DBCSG_HEX_6F            0x07010000
#define DBCSG_TET_4F            0x07020000
#define DBCSG_PYRAMID_5F        0x07030000
#define DBCSG_PRISM_5F          0x07040000

/* Definitions for 2D CSG boundary types */
#define DBCSG_QUADRATIC_G       0x08000000
#define DBCSG_CIRCLE_PR         0x09000000
#define DBCSG_ELLIPSE_PRR       0x09010000
#define DBCSG_LINE_G            0x0A000000
#define DBCSG_LINE_X            0x0A010000
#define DBCSG_LINE_Y            0x0A020000
#define DBCSG_LINE_PN           0x0A030000
#define DBCSG_LINE_PP           0x0A040000
#define DBCSG_BOX_XYXY          0x0B000000
#define DBCSG_ANGLE_PNLA        0x0C000000
#define DBCSG_ANGLE_PPA         0x0C010000
#define DBCSG_POLYGON_KP        0x0D000000
#define DBCSG_TRI_3P            0x0D010000
#define DBCSG_QUAD_4P           0x0D020000

/* Definitions for CSG Region operators */
#define DBCSG_INNER             0x7F000000
#define DBCSG_OUTER             0x7F010000
#define DBCSG_ON                0x7F020000
#define DBCSG_UNION             0x7F030000
#define DBCSG_INTERSECT         0x7F040000
#define DBCSG_DIFF              0x7F050000
#define DBCSG_COMPLIMENT        0x7F060000
#define DBCSG_XFORM             0x7F070000
#define DBCSG_SWEEP             0x7F080000

/* definitions for MRG Tree traversal flags */
#define DB_PREORDER             0x00000001
#define DB_POSTORDER            0x00000002
#define DB_FROMCWR              0x00000004

/* Miscellaneous constants */
#define     DB_F77NULL  (-99)   /*Fortran NULL pointer      */
#define     DB_F77NULLSTRING  "NULLSTRING"  /* FORTRAN STRING */

/*-------------------------------------------------------------------------
 * Index selection macros
 *-------------------------------------------------------------------------
 */
#define I4D(s,i,j,k,l) (l)*s[3]+(k)*s[2]+(j)*s[1]+(i)*s[0]
#define I3D(s,i,j,k)   (k)*s[2]+(j)*s[1]+(i)*s[0]
#define I2D(s,i,j)     (j)*s[1]+(i)*s[0]

/*-------------------------------------------------------------------------
 * Structures (just the public parts).
 *-------------------------------------------------------------------------
 */

/*
 * Database table of contents for the current directory only.
 */
typedef struct DBtoc_ {

    char         **curve_names;
    int            ncurve;

    char         **multimesh_names;
    int            nmultimesh;

    char         **multimeshadj_names;
    int            nmultimeshadj;

    char         **multivar_names;
    int            nmultivar;

    char         **multimat_names;
    int            nmultimat;

    char         **multimatspecies_names;
    int            nmultimatspecies;

    char         **csgmesh_names;
    int            ncsgmesh;

    char         **csgvar_names;
    int            ncsgvar;

    char         **defvars_names;
    int            ndefvars;

    char         **qmesh_names;
    int            nqmesh;

    char         **qvar_names;
    int            nqvar;

    char         **ucdmesh_names;
    int            nucdmesh;

    char         **ucdvar_names;
    int            nucdvar;

    char         **ptmesh_names;
    int            nptmesh;

    char         **ptvar_names;
    int            nptvar;

    char         **mat_names;
    int            nmat;

    char         **matspecies_names;
    int            nmatspecies;

    char         **var_names;
    int            nvar;

    char         **obj_names;
    int            nobj;

    char         **dir_names;
    int            ndir;

    char         **array_names;
    int            narrays;

    char         **mrgtree_names;
    int            nmrgtrees;

    char         **groupelmap_names;
    int            ngroupelmaps;

    char         **mrgvar_names;
    int            nmrgvars;

} DBtoc;

/*----------------------------------------------------------------------------
 * Database Curve Object
 *--------------------------------------------------------------------------
 */
typedef struct DBcurve_ {
/*----------- X vs. Y (Curve) Data -----------*/
    int            id;          /* Identifier for this object */
    int            datatype;    /* Datatype for x and y (float, double) */
    int            origin;      /* '0' or '1' */
    char          *title;       /* Title for curve */
    char          *xvarname;    /* Name of domain (x) variable */
    char          *yvarname;    /* Name of range  (y) variable */
    char          *xlabel;      /* Label for x-axis */
    char          *ylabel;      /* Label for y-axis */
    char          *xunits;      /* Units for domain */
    char          *yunits;      /* Units for range  */
    DB_DTPTR      *x;           /* Domain values for curve */
    DB_DTPTR      *y;           /* Range  values for curve */
    int            npts;        /* Number of points in curve */
    int            guihide;     /* Flag to hide from post-processor's GUI */
    char          *reference;   /* Label to reference object */
} DBcurve;

typedef struct DBdefvars_ {
    int            ndefs;       /* number of definitions */
    char         **names;       /* [ndefs] derived variable names */
    int           *types;       /* [ndefs] derived variable types */
    char         **defns;       /* [ndefs] derived variable definitions */
    int        *guihides;       /* [ndefs] flags to hide from
                                   post-processor's GUI */
} DBdefvars;

typedef struct DBpointmesh_ {
/*----------- Point Mesh -----------*/
    int            id;          /* Identifier for this object */
    int            block_no;    /* Block number for this mesh */
    int            group_no;    /* Block group number for this mesh */
    char          *name;        /* Name associated with this mesh */
    int            cycle;       /* Problem cycle number */
    char          *units[3];    /* Units for each axis */
    char          *labels[3];   /* Labels for each axis */
    char          *title;       /* Title for curve */

    DB_DTPTR      *coords[3];   /* Coordinate values */
    float          time;        /* Problem time */
    double         dtime;       /* Problem time, double data type */
   /*
    * The following two fields really only contain 3 elements.  However, silo
    * contains a bug in PJ_ReadVariable() as called by DBGetPointmesh() which
    * can cause three doubles to be stored there instead of three floats.
    */
    float          min_extents[6];  /* Min mesh extents [ndims] */
    float          max_extents[6];  /* Max mesh extents [ndims] */

    int            datatype;    /* Datatype for coords (float, double) */
    int            ndims;       /* Number of computational dimensions */
    int            nels;        /* Number of elements in mesh */
    int            origin;      /* '0' or '1' */
    int            guihide;     /* Flag to hide from post-processor's GUI */
    void          *gnodeno;     /* global node ids */
    char          *mrgtree_name; /* optional name of assoc. mrgtree object */
    int            gnznodtype;  /* datatype for global node/zone ids */
} DBpointmesh;

/*----------------------------------------------------------------------------
 * Multi-Block Mesh Object
 *--------------------------------------------------------------------------
 */
typedef struct DBmultimesh_ {
/*----------- Multi-Block Mesh -----------*/
    int            id;          /* Identifier for this object */
    int            nblocks;     /* Number of blocks in mesh */
    int            ngroups;     /* Number of block groups in mesh */
    int           *meshids;     /* Array of mesh-ids which comprise mesh */
    char         **meshnames;   /* Array of mesh-names for meshids */
    int           *meshtypes;   /* Array of mesh-type indicators [nblocks] */
    int           *dirids;      /* Array of directory ID's which contain blk */
    int            blockorigin; /* Origin (0 or 1) of block numbers */
    int            grouporigin; /* Origin (0 or 1) of group numbers */
    int            extentssize; /* size of each extent tuple */
    double        *extents;     /* min/max extents of coords of each block */
    int           *zonecounts;  /* array of zone counts for each block */
    int           *has_external_zones;  /* external flags for each block */
    int            guihide;     /* Flag to hide from post-processor's GUI */
    int            lgroupings;  /* size of groupings array */
    int            *groupings;  /* Array of mesh-ids, group-ids, and counts */
    char          **groupnames; /* Array of group-names for groupings  */
    char          *mrgtree_name;/* optional name of assoc. mrgtree object */
    int            tv_connectivity;
    int            disjoint_mode;
    int            topo_dim;    /* Topological dimension; max of all blocks. */ 
    char          *file_ns;     /* namescheme for files (in lieu of meshnames) */
    char          *block_ns;    /* namescheme for block objects (in lieu of meshnames) */
    int            block_type;  /* constant block type for all blocks (in lieu of meshtypes) */
    int           *empty_list;  /* list of empty block #'s (option for namescheme) */
    int            empty_cnt;   /* size of empty list */
} DBmultimesh;

/*----------------------------------------------------------------------------
 * Multi-Block Mesh Adjacency Object
 *--------------------------------------------------------------------------
 */
typedef struct DBmultimeshadj_ {
/*----------- Multi-Block Mesh Adjacency -----------*/
    int            nblocks;     /* Number of blocks in mesh */
    int            blockorigin; /* Origin (0 or 1) of block numbers */
    int           *meshtypes;   /* Array of mesh-type indicators [nblocks] */
    int           *nneighbors;  /* Array [nblocks] neighbor counts */

    int           lneighbors;
    int           *neighbors;   /* Array [lneighbors] neighbor block numbers */
    int           *back;        /* Array [lneighbors] neighbor block back */

    int            totlnodelists;
    int           *lnodelists;  /* Array [lneighbors] of node counts shared */
    int          **nodelists;   /* Array [lneighbors] nodelists shared */

    int            totlzonelists;
    int           *lzonelists;  /* Array [lneighbors] of zone counts adjacent */
    int          **zonelists;   /* Array [lneighbors] zonelists adjacent */
} DBmultimeshadj;

/*----------------------------------------------------------------------------
 * Multi-Block Variable Object
 *--------------------------------------------------------------------------
 */
typedef struct DBmultivar_ {
/*----------- Multi-Block Variable -----------*/
    int            id;          /* Identifier for this object  */
    int            nvars;       /* Number of variables   */
    int            ngroups;     /* Number of block groups in mesh */
    char         **varnames;    /* Variable names   */
    int           *vartypes;    /* variable types   */
    int            blockorigin; /* Origin (0 or 1) of block numbers */
    int            grouporigin; /* Origin (0 or 1) of group numbers */
    int            extentssize; /* size of each extent tuple */
    double        *extents;     /* min/max extents of each block */
    int            guihide;     /* Flag to hide from post-processor's GUI */
    char         **region_pnames;
    char          *mmesh_name;
    int            tensor_rank;    /* DB_VARTYPE_XXX */
    int            conserved;   /* indicates if the variable should be conserved
                                   under various operations such as interp. */
    int            extensive;   /* indicates if the variable reprsents an extensiv
                                   physical property (as opposed to intensive) */
    char          *file_ns;     /* namescheme for files (in lieu of meshnames) */
    char          *block_ns;    /* namescheme for block objects (in lieu of meshnames) */
    int            block_type;  /* constant block type for all blocks (in lieu of meshtypes) */
    int           *empty_list;  /* list of empty block #'s (option for namescheme) */
    int            empty_cnt;   /* size of empty list */
} DBmultivar;

/*-------------------------------------------------------------------------
 * Multi-material
 *-------------------------------------------------------------------------
 */
typedef struct DBmultimat_ {
    int            id;          /* Identifier for this object  */
    int            nmats;       /* Number of materials   */
    int            ngroups;     /* Number of block groups in mesh */
    char         **matnames;    /* names of constiuent DBmaterial objects */
    int            blockorigin; /* Origin (0 or 1) of block numbers */
    int            grouporigin; /* Origin (0 or 1) of group numbers */
    int           *mixlens;     /* array of mixlen values in each mat */
    int           *matcounts;   /* counts of unique materials in each block */
    int           *matlists;    /* list of materials in each block */
    int            guihide;     /* Flag to hide from post-processor's GUI */
    int            nmatnos;     /* global number of materials over all pieces */
    int           *matnos;      /* global list of material numbers */
    char         **matcolors;   /* optional colors for materials */
    char         **material_names; /* optional names of the materials */
    int            allowmat0;   /* Flag to allow material "0" */
    char          *mmesh_name;
    char          *file_ns;     /* namescheme for files (in lieu of meshnames) */
    char          *block_ns;    /* namescheme for block objects (in lieu of meshnames) */
    int           *empty_list;  /* list of empty block #'s (option for namescheme) */
    int            empty_cnt;   /* size of empty list */
} DBmultimat;

/*-------------------------------------------------------------------------
 * Multi-species
 *-------------------------------------------------------------------------
 */
typedef struct DBmultimatspecies_ {
    int            id;          /* Identifier for this object  */
    int            nspec;       /* Number of species   */
    int            ngroups;     /* Number of block groups in mesh */
    char         **specnames;   /* Species object names   */    
    int            blockorigin; /* Origin (0 or 1) of block numbers */
    int            grouporigin; /* Origin (0 or 1) of group numbers */
    int            guihide;     /* Flag to hide from post-processor's GUI */
    int            nmat;        /* equiv. to nmatnos of a DBmultimat */
    int           *nmatspec;    /* equiv. to matnos of a DBmultimat */
    char         **species_names; /* optional names of the species */
    char         **speccolors;  /* optional colors for species */
    char          *file_ns;     /* namescheme for files (in lieu of meshnames) */
    char          *block_ns;    /* namescheme for block objects (in lieu of meshnames) */
    int           *empty_list;  /* list of empty block #'s (option for namescheme) */
    int            empty_cnt;   /* size of empty list */
} DBmultimatspecies;

/*----------------------------------------------------------------------
 *  Definitions for the FaceList, ZoneList, and EdgeList structures
 *  used for describing UCD meshes.
 *----------------------------------------------------------------------
 */

#define DB_ZONETYPE_BEAM        10

#define DB_ZONETYPE_POLYGON     20
#define DB_ZONETYPE_TRIANGLE    23
#define DB_ZONETYPE_QUAD        24

#define DB_ZONETYPE_POLYHEDRON  30
#define DB_ZONETYPE_TET         34
#define DB_ZONETYPE_PYRAMID     35
#define DB_ZONETYPE_PRISM       36
#define DB_ZONETYPE_HEX         38

typedef struct DBzonelist_ {
    int            ndims;       /* Number of dimensions (2,3) */
    int            nzones;      /* Number of zones in list */
    int            nshapes;     /* Number of zone shapes */
    int           *shapecnt;    /* [nshapes] occurences of each shape */
    int           *shapesize;   /* [nshapes] Number of nodes per shape */
    int           *shapetype;   /* [nshapes] Type of shape */
    int           *nodelist;    /* Sequent lst of nodes which comprise zones */
    int            lnodelist;   /* Number of nodes in nodelist */
    int            origin;      /* '0' or '1' */
    int            min_index;   /* Index of first real zone */
    int            max_index;   /* Index of last real zone */

/*--------- Optional zone attributes ---------*/
    int           *zoneno;      /* [nzones] zone number of each zone */
    void          *gzoneno;     /* [nzones] global zone number of each zone */
    int            gnznodtype;  /* datatype for global node/zone ids */
} DBzonelist;

typedef struct DBphzonelist_ {
    int            nfaces;      /* Number of faces in facelist (aka "facetable") */
    int           *nodecnt;     /* Count of nodes in each face */
    int            lnodelist;   /* Length of nodelist used to construct faces */
    int           *nodelist;    /* List of nodes used in all faces */
    char          *extface;     /* boolean flag indicating if a face is external */
    int            nzones;      /* Number of zones in this zonelist */
    int           *facecnt;     /* Count of faces in each zone */
    int            lfacelist;   /* Length of facelist used to construct zones */
    int           *facelist;    /* List of faces used in all zones */
    int            origin;      /* '0' or '1' */
    int            lo_offset;   /* Index of first non-ghost zone */
    int            hi_offset;   /* Index of last non-ghost zone */

/*--------- Optional zone attributes ---------*/
    int           *zoneno;      /* [nzones] zone number of each zone */
    void          *gzoneno;     /* [nzones] global zone number of each zone */
    int            gnznodtype;  /* datatype for global node/zone ids */
} DBphzonelist;

typedef struct DBfacelist_ {
/*----------- Required components ------------*/
    int            ndims;       /* Number of dimensions (2,3) */
    int            nfaces;      /* Number of faces in list */
    int            origin;      /* '0' or '1' */
    int           *nodelist;    /* Sequent list of nodes comprise faces */
    int            lnodelist;   /* Number of nodes in nodelist */

/*----------- 3D components ------------------*/
    int            nshapes;     /* Number of face shapes */
    int           *shapecnt;    /* [nshapes] Num of occurences of each shape */
    int           *shapesize;   /* [nshapes] Number of nodes per shape */

/*----------- Optional type component---------*/
    int            ntypes;      /* Number of face types */
    int           *typelist;    /* [ntypes] Type ID for each type */
    int           *types;       /* [nfaces] Type info for each face */

/*--------- Optional node attributes ---------*/
    int           *nodeno;      /* [lnodelist] node number of each node */

/*----------- Optional zone-reference component---------*/
    int           *zoneno;      /* [nfaces] Zone number for each face */
} DBfacelist;

typedef struct DBedgelist_ {
    int            ndims;       /* Number of dimensions (2,3) */
    int            nedges;      /* Number of edges */
    int           *edge_beg;    /* [nedges] */
    int           *edge_end;    /* [nedges] */
    int            origin;      /* '0' or '1' */
} DBedgelist;

typedef struct DBquadmesh_ {
/*----------- Quad Mesh -----------*/
    int            id;          /* Identifier for this object */
    int            block_no;    /* Block number for this mesh */
    int            group_no;    /* Block group number for this mesh */
    char          *name;        /* Name associated with mesh */
    int            cycle;       /* Problem cycle number */
    int            coord_sys;   /* Cartesian, cylindrical, spherical */
    int            major_order; /* 1 indicates row-major for multi-d arrays */
    int            stride[3];   /* Offsets to adjacent elements  */
    int            coordtype;   /* Coord array type: collinear,
                                 * non-collinear */
    int            facetype;    /* Zone face type: rect, curv */
    int            planar;      /* Sentinel: zones represent area or volume? */

    DB_DTPTR      *coords[3];   /* Mesh node coordinate ptrs [ndims] */
    int            datatype;    /* Type of coordinate arrays (double,float) */
    float          time;        /* Problem time */
    double         dtime;       /* Problem time, double data type */
   /*
    * The following two fields really only contain 3 elements.  However, silo
    * contains a bug in PJ_ReadVariable() as called by DBGetQuadmesh() which
    * can cause three doubles to be stored there instead of three floats.
    */
    float          min_extents[6];  /* Min mesh extents [ndims] */
    float          max_extents[6];  /* Max mesh extents [ndims] */

    char          *labels[3];   /* Label associated with each dimension */
    char          *units[3];    /* Units for variable, e.g, 'mm/ms' */
    int            ndims;       /* Number of computational dimensions */
    int            nspace;      /* Number of physical dimensions */
    int            nnodes;      /* Total number of nodes */

    int            dims[3];     /* Number of nodes per dimension */
    int            origin;      /* '0' or '1' */
    int            min_index[3];   /* Index in each dimension of 1st
                                    * non-phoney */
    int            max_index[3];   /* Index in each dimension of last
                                    * non-phoney */
    int            base_index[3];  /* Lowest real i,j,k value for this block */
    int            start_index[3]; /* i,j,k values corresponding to original
                                    * mesh */
    int            size_index[3];  /* Number of nodes per dimension for 
                                    * original mesh */
    int            guihide;     /* Flag to hide from post-processor's GUI */
    char          *mrgtree_name; /* optional name of assoc. mrgtree object */
} DBquadmesh;

typedef struct DBucdmesh_ {
/*----------- Unstructured Cell Data (UCD) Mesh -----------*/
    int            id;          /* Identifier for this object */
    int            block_no;    /* Block number for this mesh */
    int            group_no;    /* Block group number for this mesh */
    char          *name;        /* Name associated with mesh */
    int            cycle;       /* Problem cycle number */
    int            coord_sys;   /* Coordinate system */
    int            topo_dim;    /* Topological dimension. */ 
    char          *units[3];    /* Units for variable, e.g, 'mm/ms' */
    char          *labels[3];   /* Label associated with each dimension */

    DB_DTPTR      *coords[3];   /* Mesh node coordinates */
    int            datatype;    /* Type of coordinate arrays (double,float) */
    float          time;        /* Problem time */
    double         dtime;       /* Problem time, double data type */
   /*
    * The following two fields really only contain 3 elements.  However, silo
    * contains a bug in PJ_ReadVariable() as called by DBGetUcdmesh() which
    * can cause three doubles to be stored there instead of three floats.
    */
    float          min_extents[6];  /* Min mesh extents [ndims] */
    float          max_extents[6];  /* Max mesh extents [ndims] */

    int            ndims;       /* Number of computational dimensions */
    int            nnodes;      /* Total number of nodes */
    int            origin;      /* '0' or '1' */

    DBfacelist    *faces;       /* Data structure describing mesh faces */
    DBzonelist    *zones;       /* Data structure describing mesh zones */
    DBedgelist    *edges;       /* Data struct describing mesh edges
                                 * (option) */

/*--------- Optional node attributes ---------*/
    void          *gnodeno;     /* [nnodes] global node number of each node */

/*--------- Optional zone attributes ---------*/
    int           *nodeno;      /* [nnodes] node number of each node */

/*--------- Optional polyhedral zonelist ---------*/
    DBphzonelist  *phzones;     /* Data structure describing mesh zones */

    int            guihide;     /* Flag to hide from post-processor's GUI */
    char          *mrgtree_name; /* optional name of assoc. mrgtree object */
    int            tv_connectivity;
    int            disjoint_mode;
    int            gnznodtype;  /* datatype for global node/zone ids */
} DBucdmesh;

/*----------------------------------------------------------------------------
 * Database Mesh-Variable Object
 *---------------------------------------------------------------------------
 */
typedef struct DBquadvar_ {
/*----------- Quad Variable -----------*/
    int            id;          /* Identifier for this object */
    char          *name;        /* Name of variable */
    char          *units;       /* Units for variable, e.g, 'mm/ms' */
    char          *label;       /* Label (perhaps for editing purposes) */
    int            cycle;       /* Problem cycle number */
    int            meshid;      /* Identifier for associated mesh (Deprecated Sep2005) */

    DB_DTPTR     **vals;        /* Array of pointers to data arrays */
    int            datatype;    /* Type of data pointed to by 'val' */
    int            nels;        /* Number of elements in each array */
    int            nvals;       /* Number of arrays pointed to by 'vals' */
    int            ndims;       /* Rank of variable */
    int            dims[3];     /* Number of elements in each dimension */

    int            major_order; /* 1 indicates row-major for multi-d arrays */
    int            stride[3];   /* Offsets to adjacent elements  */
    int            min_index[3];  /* Index in each dimension of 1st
                                   * non-phoney */
    int            max_index[3];  /* Index in each dimension of last
                                   * non-phoney */
    int            origin;      /* '0' or '1' */
    float          time;        /* Problem time */
    double         dtime;       /* Problem time, double data type */
   /*
    * The following field really only contains 3 elements.  However, silo
    * contains a bug in PJ_ReadVariable() as called by DBGetQuadvar() which
    * can cause three doubles to be stored there instead of three floats.
    */
    float          align[6];    /* Centering and alignment per dimension */

    DB_DTPTR     **mixvals;     /* nvals ptrs to data arrays for mixed zones */
    int            mixlen;      /* Num of elmts in each mixed zone data
                                 * array */

    int            use_specmf;  /* Flag indicating whether to apply species
                                 * mass fractions to the variable. */

    int            ascii_labels;/* Treat variable values as ASCII values
                                   by rounding to the nearest integer in
                                   the range [0, 255] */
    char          *meshname;    /* Name of associated mesh */
    int            guihide;     /* Flag to hide from post-processor's GUI */
    char         **region_pnames;
    int            conserved;   /* indicates if the variable should be conserved
                                   under various operations such as interp. */
    int            extensive;   /* indicates if the variable reprsents an extensiv
                                   physical property (as opposed to intensive) */
    int            centering;   /* explicit centering knowledge; should agree
                                   with alignment. */
} DBquadvar;

typedef struct DBucdvar_ {
/*----------- Unstructured Cell Data (UCD) Variable -----------*/
    int            id;          /* Identifier for this object */
    char          *name;        /* Name of variable */
    int            cycle;       /* Problem cycle number */
    char          *units;       /* Units for variable, e.g, 'mm/ms' */
    char          *label;       /* Label (perhaps for editing purposes) */
    float          time;        /* Problem time */
    double         dtime;       /* Problem time, double data type */
    int            meshid;      /* Identifier for associated mesh (Deprecated Sep2005) */

    DB_DTPTR     **vals;        /* Array of pointers to data arrays */
    int            datatype;    /* Type of data pointed to by 'vals' */
    int            nels;        /* Number of elements in each array */
    int            nvals;       /* Number of arrays pointed to by 'vals' */
    int            ndims;       /* Rank of variable */
    int            origin;      /* '0' or '1' */

    int            centering;   /* Centering within mesh (nodal or zonal) */
    DB_DTPTR     **mixvals;     /* nvals ptrs to data arrays for mixed zones */
    int            mixlen;      /* Num of elmts in each mixed zone data
                                 * array */

    int            use_specmf;  /* Flag indicating whether to apply species
                                 * mass fractions to the variable. */
    int            ascii_labels;/* Treat variable values as ASCII values
                                   by rounding to the nearest integer in
                                   the range [0, 255] */
    char          *meshname;    /* Name of associated mesh */
    int            guihide;     /* Flag to hide from post-processor's GUI */
    char         **region_pnames;
    int            conserved;   /* indicates if the variable should be conserved
                                   under various operations such as interp. */
    int            extensive;   /* indicates if the variable reprsents an extensiv
                                   physical property (as opposed to intensive) */
} DBucdvar;

typedef struct DBmeshvar_ {
/*----------- Generic Mesh-Data Variable -----------*/
    int            id;          /* Identifier for this object */
    char          *name;        /* Name of variable */
    char          *units;       /* Units for variable, e.g, 'mm/ms' */
    char          *label;       /* Label (perhaps for editing purposes) */
    int            cycle;       /* Problem cycle number */
    int            meshid;      /* Identifier for associated mesh (Deprecated Sep2005) */

    DB_DTPTR     **vals;        /* Array of pointers to data arrays */
    int            datatype;    /* Type of data pointed to by 'val' */
    int            nels;        /* Number of elements in each array */
    int            nvals;       /* Number of arrays pointed to by 'vals' */
    int            nspace;      /* Spatial rank of variable */
    int            ndims;       /* Rank of 'vals' array(s) (computatnl rank) */

    int            origin;      /* '0' or '1' */
    int            centering;   /* Centering within mesh (nodal,zonal,other) */
    float          time;        /* Problem time */
    double         dtime;       /* Problem time, double data type */
   /*
    * The following field really only contains 3 elements.  However, silo
    * contains a bug in PJ_ReadVariable() as called by DBGetPointvar() which
    * can cause three doubles to be stored there instead of three floats.
    */
    float          align[6];    /* Alignmnt per dimension if
                                 * centering==other */

    /* Stuff for multi-dimensional arrays (ndims > 1) */
    int            dims[3];     /* Number of elements in each dimension */
    int            major_order; /* 1 indicates row-major for multi-d arrays */
    int            stride[3];   /* Offsets to adjacent elements  */
   /*
    * The following two fields really only contain 3 elements.  However, silo
    * contains a bug in PJ_ReadVariable() as called by DBGetUcdmesh() which
    * can cause three doubles to be stored there instead of three floats.
    */
    int            min_index[6];  /* Index in each dimension of 1st
                                   * non-phoney */
    int            max_index[6];  /* Index in each dimension of last
                                    non-phoney */

    int            ascii_labels;/* Treat variable values as ASCII values
                                   by rounding to the nearest integer in
                                   the range [0, 255] */
    char          *meshname;      /* Name of associated mesh */
    int            guihide;     /* Flag to hide from post-processor's GUI */
    char         **region_pnames;
    int            conserved;   /* indicates if the variable should be conserved
                                   under various operations such as interp. */
    int            extensive;   /* indicates if the variable reprsents an extensiv
                                   physical property (as opposed to intensive) */
} DBmeshvar;

typedef struct DBmaterial_ {
/*----------- Material Information -----------*/
    int            id;          /* Identifier */
    char          *name;        /* Name of this material information block */
    int            ndims;       /* Rank of 'matlist' variable */
    int            origin;      /* '0' or '1' */
    int            dims[3];     /* Number of elements in each dimension */
    int            major_order; /* 1 indicates row-major for multi-d arrays */
    int            stride[3];   /* Offsets to adjacent elements in matlist */

    int            nmat;        /* Number of materials */
    int           *matnos;      /* Array [nmat] of valid material numbers */
    char         **matnames;    /* Array of material names   */
    int           *matlist;     /* Array[nzone] w/ mat. number or mix index */
    int            mixlen;      /* Length of mixed data arrays (mix_xxx) */
    int            datatype;    /* Type of volume-fractions (double,float) */
    DB_DTPTR      *mix_vf;      /* Array [mixlen] of volume fractions */
    int           *mix_next;    /* Array [mixlen] of mixed data indeces */
    int           *mix_mat;     /* Array [mixlen] of material numbers */
    int           *mix_zone;    /* Array [mixlen] of back pointers to mesh */

    char         **matcolors;   /* Array of material colors */
    char          *meshname;    /* Name of associated mesh */
    int            allowmat0;   /* Flag to allow material "0" */
    int            guihide;     /* Flag to hide from post-processor's GUI */
} DBmaterial;

typedef struct DBmatspecies_ {
/*----------- Species Information -----------*/
    int            id;          /* Identifier */
    char          *name;        /* Name of this matspecies information block */
    char          *matname;     /* Name of material object with which the
                                 * material species object is associated. */
    int            nmat;        /* Number of materials */
    int           *nmatspec;    /* Array of lngth nmat of the num of material
                                 * species associated with each material. */
    int            ndims;       /* Rank of 'speclist' variable */
    int            dims[3];     /* Number of elements in each dimension of the
                                 * 'speclist' variable. */
    int            major_order; /* 1 indicates row-major for multi-d arrays */
    int            stride[3];   /* Offsts to adjacent elmts in 'speclist'  */

    int            nspecies_mf; /* Total number of species mass fractions. */
    DB_DTPTR      *species_mf;  /* Array of length nspecies_mf of mass
                                 * frations of the material species. */
    int           *speclist;    /* Zone array of dimensions described by ndims
                                 * and dims.  Each element of the array is an
                                 * index into one of the species mass fraction
                                 * arrays.  A positive value is the index in
                                 * the species_mf array of the mass fractions
                                 * of the clean zone's material species:
                                 * species_mf[speclist[i]] is the mass fraction
                                 * of the first species of material matlist[i]
                                 * in zone i. A negative value means that the
                                 * zone is a mixed zone and that the array
                                 * mix_speclist contains the index to the
                                 * species mas fractions: -speclist[i] is the
                                 * index in the 'mix_speclist' array for zone
                                 * i. */
    int            mixlen;      /* Length of 'mix_speclist' array. */
    int           *mix_speclist;  /* Array of lgth mixlen of 1-orig indices
                                   * into the 'species_mf' array.
                                   * species_mf[mix_speclist[j]] is the index
                                   * in array species_mf' of the first of the
                                   * mass fractions for material
                                   * mix_mat[j]. */

    int            datatype;    /* Datatype of mass fraction data. */
    int            guihide;     /* Flag to hide from post-processor's GUI */
    char         **specnames;   /* Array of species names; length is sum of nmatspec   */
    char         **speccolors;  /* Array of species colors; length is sum of nmatspec */
} DBmatspecies;

typedef struct DBcsgzonelist_ {
/*----------- CSG Zonelist -----------*/
    int            nregs;       /* Number of regions in regionlist */
    int            origin;      /* '0' or '1' */

    int           *typeflags;   /* [nregs] type info about each region */
    int           *leftids;     /* [nregs] left operand region refs */
    int           *rightids;    /* [nregs] right operand region refs */
    void          *xform;       /* [lxforms] transformation coefficients */
    int            lxform;      /* length of xforms array */
    int            datatype;    /* type of data in xforms array */

    int            nzones;      /* number of zones */
    int           *zonelist;    /* [nzones] region ids (complete regions) */
    int            min_index;   /* Index of first real zone */
    int            max_index;   /* Index of last real zone */

/*--------- Optional zone attributes ---------*/
    char         **regnames;   /* [nregs] names of each region */
    char         **zonenames;  /* [nzones] names of each zone */
} DBcsgzonelist;

typedef struct DBcsgmesh_ {
/*----------- CSG Mesh -----------*/
    int            block_no;    /* Block number for this mesh */
    int            group_no;    /* Block group number for this mesh */
    char          *name;        /* Name associated with mesh */
    int            cycle;       /* Problem cycle number */
    char          *units[3];    /* Units for variable, e.g, 'mm/ms' */
    char          *labels[3];   /* Label associated with each dimension */

    int            nbounds;     /* Total number of boundaries */
    int           *typeflags;   /* nbounds boundary type info flags */
    int           *bndids;      /* optional, nbounds explicit ids */

    void          *coeffs;      /* coefficients in the representation of
                                   each boundary */
    int            lcoeffs;     /* length of coeffs array */
    int           *coeffidx;    /* array of nbounds offsets into coeffs
                                   for each boundary's coefficients */
    int            datatype;    /* data type of coeffs data */

    float          time;        /* Problem time */
    double         dtime;       /* Problem time, double data type */
    double         min_extents[3];  /* Min mesh extents [ndims] */
    double         max_extents[3];  /* Max mesh extents [ndims] */

    int            ndims;       /* Number of spatial & topological dimensions */
    int            origin;      /* '0' or '1' */

    DBcsgzonelist *zones;       /* Data structure describing mesh zones */

/*--------- Optional boundary attributes ---------*/
    char         **bndnames;     /* [nbounds] boundary names */
    int            guihide;     /* Flag to hide from post-processor's GUI */
    char          *mrgtree_name; /* optional name of assoc. mrgtree object */
    int            tv_connectivity;
    int            disjoint_mode;
} DBcsgmesh;

typedef struct DBcsgvar_ {
/*----------- CSG Variable -----------*/
    char          *name;        /* Name of variable */
    int            cycle;       /* Problem cycle number */
    char          *units;       /* Units for variable, e.g, 'mm/ms' */
    char          *label;       /* Label (perhaps for editing purposes) */
    float          time;        /* Problem time */
    double         dtime;       /* Problem time, double data type */

    void         **vals;        /* Array of pointers to data arrays */
    int            datatype;    /* Type of data pointed to by 'vals' */
    int            nels;        /* Number of elements in each array */
    int            nvals;       /* Number of arrays pointed to by 'vals' */

    int            centering;   /* Centering within mesh (nodal or zonal) */

    int            use_specmf;  /* Flag indicating whether to apply species
                                 * mass fractions to the variable. */
    int            ascii_labels;/* Treat variable values as ASCII values
                                   by rounding to the nearest integer in
                                   the range [0, 255] */
    char          *meshname;    /* Name of associated mesh */
    int            guihide;     /* Flag to hide from post-processor's GUI */
    char         **region_pnames;
    int            conserved;   /* indicates if the variable should be conserved
                                   under various operations such as interp. */
    int            extensive;   /* indicates if the variable reprsents an extensiv
                                   physical property (as opposed to intensive) */
} DBcsgvar;

/*-------------------------------------------------------------------------
 * A compound array is an array whose elements are simple arrays. A simple
 * array is an array whose elements are all of the same primitive data
 * type: float, double, integer, long...  All of the simple arrays of
 * a compound array have elements of the same data type.
 *-------------------------------------------------------------------------
 */
typedef struct DBcompoundarray_ {
    int            id;          /*identifier of the compound array */
    char          *name;        /*name of te compound array  */
    char         **elemnames;   /*names of the simple array elements */
    int           *elemlengths; /*lengths of the simple arrays  */
    int            nelems;      /*number of simple arrays  */
    void          *values;      /*simple array values   */
    int            nvalues;     /*sum reduction of `elemlengths' vector */
    int            datatype;    /*simple array element data type */
} DBcompoundarray;

typedef struct DBoptlist_ {

    int           *options;     /* Vector of option identifiers */
    void         **values;      /* Vector of pointers to option values */
    int            numopts;     /* Number of options defined */
    int            maxopts;     /* Total length of option/value arrays */

} DBoptlist;

typedef struct DBobject_ {

    char          *name;
    char          *type;        /* Type of group/object */
    char         **comp_names;  /* Array of component names */
    char         **pdb_names;   /* Array of internal (PDB) variable names */
    int            ncomponents; /* Number of components */
    int            maxcomponents;  /* Max number of components */

} DBobject;

typedef struct _DBmrgtnode {
    char *name;
    int  narray;
    char **names;
    int type_info_bits;
    int max_children;
    char *maps_name;
    int nsegs;
    int *seg_ids;
    int *seg_lens;
    int *seg_types;
    int num_children;
    struct _DBmrgtnode **children;

    /* internal stuff to support updates, i/o, etc. */
    int walk_order;
    struct _DBmrgtnode  *parent;
} DBmrgtnode;

typedef void (*DBmrgwalkcb)(DBmrgtnode *tnode, int nat_node_num, void *data);

typedef struct _DBmrgtree {
    char *name;
    char *src_mesh_name;
    int src_mesh_type;
    int type_info_bits;
    int num_nodes;
    DBmrgtnode *root;
    DBmrgtnode *cwr;

    char **mrgvar_onames;
    char **mrgvar_rnames;
} DBmrgtree;

typedef struct _DBmrgvar {
    char *name;
    char *mrgt_name;
    int ncomps;
    char **compnames;
    int nregns;
    char **reg_pnames;
    int datatype;
    void **data;
} DBmrgvar ;

typedef struct _DBgroupelmap {
    char *name;
    int num_segments;
    int *groupel_types;
    int *segment_lengths;
    int *segment_ids;
    int **segment_data;
    void **segment_fracs;
    int fracs_data_type;
} DBgroupelmap;

#if !defined(DB_MAX_EXPSTRS) /* NO_FORTRAN_DEFINE */
#define DB_MAX_EXPSTRS 8 /* NO_FORTRAN_DEFINE */
#endif

typedef struct _DBnamescheme
{
    char *fmt;              /* orig. format string */
    const char **fmtptrs;   /* ptrs into first (printf) part of fmt for each conversion spec. */
    int fmtlen;             /* len of first part of fmt */
    int ncspecs;            /* # of conversion specs in first part of fmt */
    char delim;             /* delimiter char used for parts of fmt */
    int nembed;             /* number of last embedded string encountered */
    char *embedstrs[DB_MAX_EXPSTRS]; /* ptrs to copies of embedded strings */
    int narrefs;            /* number of array refs in conversion specs */
    char **arrnames;        /* array names used by array refs */
    const int **arrvals;    /* pointer to actual array data assoc. with each name */
    char **exprstrs;        /* expressions to be evaluated for each conv. spec. */
} DBnamescheme;

typedef struct DBfile *___DUMMY_TYPE;  /* Satisfy ANSI scope rules */

/*
 * All file formats are now anonymous except for the public properties
 * and public methods.
 */
typedef struct DBfile_pub {

    /* Public Properties */
    char          *name;        /*name of file    */
    int            type;        /*file type    */
    DBtoc         *toc;         /*table of contents   */
    int            dirid;       /*directory ID    */
    int            fileid;      /*unique file id [0,DB_NFILES-1] */
    int            pathok;      /*driver handles paths in names */
    int            Grab;        /*drive has access to low-level interface */
    void          *GrabId;      /*pointer to low-level driver descriptor */
    char          *file_lib_version; /* version of lib file was created with */

    /* Public Methods */
    int            (*close)(struct DBfile *);
    int            (*exist)(struct DBfile *, char *);
    int            (*pause)(struct DBfile *);
    int            (*cont)(struct DBfile *);
    int            (*newtoc)(struct DBfile *);
    DBObjectType   (*inqvartype)(struct DBfile *, char*);
    int            (*uninstall)(struct DBfile *);
    DBobject      *(*g_obj)(struct DBfile *, char *);
    int            (*c_obj)(struct DBfile *, DBobject *, int);
    int            (*w_obj)(struct DBfile *, DBobject *, int);
    void          *(*g_comp)(struct DBfile *, char *, char *);
    int            (*g_comptyp)(struct DBfile *, char *, char *);
    int            (*w_comp)(struct DBfile *, DBobject *, char *, char *,
                             char *, const void *, int, long *);
    int            (*write) (struct DBfile *, char *, void *, int *, int, int);
    int            (*writeslice)(struct DBfile *, char *, void *, int,
                                 int[], int[], int[], int[], int);
    void          *(*g_attr)(struct DBfile *, char *, char *);
    int            (*g_dir)(struct DBfile *, char *);
    int            (*mkdir)(struct DBfile *, char *);
    int            (*cd)(struct DBfile *, char *);
    int            (*cdid)(struct DBfile *, int);
    int            (*r_att)(struct DBfile *, char *, char *, void *);
    int            (*r_var)(struct DBfile *, char *, void *);
    int            (*r_var1)(struct DBfile *, char *, int, void *);
    int            (*module)(struct DBfile *, FILE *);
    int            (*r_varslice)(struct DBfile *, char *, int *, int *, int *,
                                 int, void *);
    int            (*g_compnames)(struct DBfile *, char *, char ***, char ***);
    DBcompoundarray *(*g_ca)(struct DBfile *, char *);
    DBcurve       *(*g_cu)(struct DBfile *, char *);
    DBdefvars     *(*g_defv)(struct DBfile *, const char *);
    DBmaterial    *(*g_ma)(struct DBfile *, char *);
    DBmatspecies  *(*g_ms)(struct DBfile *, char *);
    DBmultimesh   *(*g_mm)(struct DBfile *, char *);
    DBmultivar    *(*g_mv)(struct DBfile *, char *);
    DBmultimat    *(*g_mt)(struct DBfile *, char *);
    DBmultimatspecies *(*g_mms)(struct DBfile *, char *);
    DBpointmesh   *(*g_pm)(struct DBfile *, char *);
    DBmeshvar     *(*g_pv)(struct DBfile *, char *);
    DBquadmesh    *(*g_qm)(struct DBfile *, char *);
    DBquadvar     *(*g_qv)(struct DBfile *, char *);
    DBucdmesh     *(*g_um)(struct DBfile *, char *);
    DBucdvar      *(*g_uv)(struct DBfile *, char *);
    DBfacelist    *(*g_fl)(struct DBfile *, char *);
    DBzonelist    *(*g_zl)(struct DBfile *, char *);
    void          *(*g_var)(struct DBfile *, char *);
    int            (*g_varbl)(struct DBfile *, char *);  /*byte length */
    int            (*g_varlen)(struct DBfile *, char *);  /*nelems */
    int            (*g_vardims)(struct DBfile*, char*, int, int*); /*dims*/
    int            (*g_vartype)(struct DBfile *, char *);
    int            (*i_meshname)(struct DBfile *, char *, char *);
    int            (*i_meshtype)(struct DBfile *, char *);
    int            (*p_ca)(struct DBfile *, char *, char **, int *, int,
                           void *, int, int, DBoptlist *);
    int            (*p_cu)(struct DBfile *, char *, void *, void *, int, int,
                           DBoptlist *);
    int            (*p_defv)(struct DBfile *, const char *, int, 
                           char **, const int *, char **,
                           DBoptlist **);
    int            (*p_fl)(struct DBfile *, char *, int, int, int *, int, int,
                           int *, int *, int *, int, int *, int *, int);
    int            (*p_ma)(struct DBfile *, char *, char *, int, int *, int *,
                           int *, int, int *, int *, int *, DB_DTPTR1, int, int,
                           DBoptlist *);
    int            (*p_ms)(struct DBfile *, char *, char *, int, int *, int *,
                           int *, int, int, DB_DTPTR1, int *, int, int,
                           DBoptlist *);
    int            (*p_mm)(struct DBfile *, char *, int, char **, int *,
                           DBoptlist *);
    int            (*p_mv)(struct DBfile *, char *, int, char **, int *,
                           DBoptlist *);
    int            (*p_mt)(struct DBfile *, char *, int, char **, DBoptlist *);
    int            (*p_mms)(struct DBfile *, char *, int, char **, DBoptlist *);
    int            (*p_pm)(struct DBfile *, char *, int, DB_DTPTR2, int, int,
                           DBoptlist *);
    int            (*p_pv)(struct DBfile *, char *, char *, int, DB_DTPTR2, int,
                           int, DBoptlist *);
    int            (*p_qm)(struct DBfile *, char *, char **, DB_DTPTR2, int *,
                           int, int, int, DBoptlist *);
    int            (*p_qv)(struct DBfile *, char *, char *, int, char **,
                           DB_DTPTR2, int *, int, DB_DTPTR2, int, int, int,
                           DBoptlist *);
    int            (*p_um)(struct DBfile *, char *, int, char **, DB_DTPTR2,
                           int, int, char *, char *, int, DBoptlist *);
    int            (*p_sm)(struct DBfile *, char *, char *,
                           int, char *, char *, DBoptlist *);
    int            (*p_uv)(struct DBfile *, char *, char *, int, char **,
                           DB_DTPTR2, int, DB_DTPTR2, int, int, int,
                           DBoptlist *);
    int            (*p_zl)(struct DBfile *, char *, int, int, int *, int, int,
                           int *, int *, int);
    int            (*p_zl2)(struct DBfile *, char *, int, int, int *, int, int,
                            int, int, int *, int *, int *, int, DBoptlist *);
    /* MCM-27Jul04: We added these to the end to avert potential
       link-time compatibility issues with older versions of the
       library. Some user's of Silo circumvent its version check
       which would ordinarily keep different versions from being
       mixed by defining an appropriate global symbol. */
    DBphzonelist  *(*g_phzl)(struct DBfile *, char *);
    int            (*p_phzl)(struct DBfile *, char *,
                             int, int *, int, int *, char *,
                             int, int *, int, int *,
                             int, int, int,
                             DBoptlist *);
    int            (*p_csgzl)(struct DBfile *, const char *, int, const int *,
                              const int *, const int *, const void *, int, int,
                              int, const int *, DBoptlist *);
    DBcsgzonelist *(*g_csgzl)(struct DBfile *, const char *);
    int            (*p_csgm)(struct DBfile *, const char *, int, int,
                             const int *, const int *,
                             const void *, int, int, const double *,
                             const char *, DBoptlist *);
    DBcsgmesh     *(*g_csgm)(struct DBfile *, const char *);
    int            (*p_csgv)(struct DBfile *, const char *, const char *, int,
                             char **, void **, int, int, int,
                             DBoptlist *);
    DBcsgvar      *(*g_csgv)(struct DBfile *, const char *);
    DBmultimeshadj *(*g_mmadj)(struct DBfile *, const char *, int, const int *);
    int            (*p_mmadj)(struct DBfile *, const char *, int, const int *,
                              const int *, const int *, const int *, const int *,
                              int **, const int *, int **,
                              DBoptlist *optlist);
    int            (*p_mrgt)(struct DBfile *dbfile, const char *name, const char *mesh_name,
                             DBmrgtree *tree, DBoptlist *opts);
    DBmrgtree     *(*g_mrgt)(struct DBfile *, const char *name);
    int            (*p_grplm)(struct DBfile *dbfile, const char *map_name,
                             int num_segments, int *groupel_types,
			     int *segment_lengths, int *segment_ids,
			     int **segment_data, void **segment_fracs,
                             int fracs_data_type, DBoptlist *opts);
    DBgroupelmap  *(*g_grplm)(struct DBfile *dbfile, const char *name);
    int            (*p_mrgv)(struct DBfile *dbfile, const char *name,
                             const char *mrgt_name,
                             int ncomps, char **compnames,
                             int nregns, char **reg_pnames,
                             int datatype, void **data, DBoptlist *opts);
    DBmrgvar      *(*g_mrgv)(struct DBfile *dbfile, const char *name);
    int            (*free_z)(struct DBfile *, const char *);
    int            (*cpdir)(struct DBfile *, const char *,
                            struct DBfile *, const char *);

    int          (*sort_obo)(struct DBfile *dbfile, int nobjs,
                             const char *const *const obj_names, int *ranks);
} DBfile_pub;

typedef struct DBfile {
    DBfile_pub     pub;
    /*private part follows per device driver */
} DBfile;

typedef void (*DBErrFunc_t)(char*);


/* The first prototypes here are the functions by which client code first
 * gets into Silo.  They are separated out because they do a version number
 * check for us.  Client code doesn't actually use these functions.
 * Instead, it uses macros like DBOpen, DBCreate, etc.
 *
 * If any functions are added that provide first-call access to Silo, they
 * should be set up as macro/function pairs, just as these are.  The way to
 * determine if a function is a "first-call" function is to ask whether
 * there are any Silo calls that must happen before it.  If there are not,
 * then the function is a "first-call" function and should have this
 * macro/function pair.  */

SILO_API extern DBfile  *DBOpenReal(const char *, int, int);
SILO_API extern DBfile  *DBCreateReal(const char *, int, int, const char *, int);
SILO_API extern int      DBInqFileReal(const char *);

#define SILO_VSTRING_NAME "_silolibinfo"
#define SILO_VSTRING PACKAGE_VERSION
SILO_API extern int SILO_VERS_TAG;
#define CheckVersion SILO_VERS_TAG = 1

#define DBOpen(name, target, mode) \
    (CheckVersion, DBOpenReal(name, target, mode))

#define DBCreate(name, mode, target, info, type) \
    (CheckVersion, DBCreateReal(name, mode, target, info, type))

#define DBInqFile(name) \
    (CheckVersion, DBInqFileReal(name))

/* Prototypes for regular API functions. */
SILO_API extern DBcompoundarray *DBAllocCompoundarray(void);
SILO_API extern DBcurve *DBAllocCurve(void);
SILO_API extern DBdefvars *DBAllocDefvars(int);
SILO_API extern DBmultimesh *DBAllocMultimesh(int);
SILO_API extern DBmultimeshadj *DBAllocMultimeshadj(int);
SILO_API extern DBmultivar *DBAllocMultivar(int);
SILO_API extern DBmultimat *DBAllocMultimat(int);
SILO_API extern DBmultimatspecies *DBAllocMultimatspecies(int);
SILO_API extern DBcsgmesh *DBAllocCsgmesh(void);
SILO_API extern DBquadmesh *DBAllocQuadmesh(void);
SILO_API extern DBpointmesh *DBAllocPointmesh(void);
SILO_API extern DBmeshvar *DBAllocMeshvar(void);
SILO_API extern DBucdmesh *DBAllocUcdmesh(void);
SILO_API extern DBcsgvar *DBAllocCsgvar(void);
SILO_API extern DBquadvar *DBAllocQuadvar(void);
SILO_API extern DBucdvar *DBAllocUcdvar(void);
SILO_API extern DBzonelist *DBAllocZonelist(void);
SILO_API extern DBphzonelist *DBAllocPHZonelist(void);
SILO_API extern DBcsgzonelist *DBAllocCSGZonelist(void);
SILO_API extern DBedgelist *DBAllocEdgelist(void);
SILO_API extern DBfacelist *DBAllocFacelist(void);
SILO_API extern DBmaterial *DBAllocMaterial(void);
SILO_API extern DBmatspecies *DBAllocMatspecies(void);
SILO_API extern DBnamescheme *DBAllocNamescheme(void);
SILO_API extern DBgroupelmap *DBAllocGroupelmap(int, DBdatatype);

SILO_API extern void     DBFreeMatspecies(DBmatspecies *);
SILO_API extern void     DBFreeMaterial(DBmaterial *);
SILO_API extern void     DBFreeFacelist(DBfacelist *);
SILO_API extern void     DBFreeEdgelist(DBedgelist *);
SILO_API extern void     DBFreeZonelist(DBzonelist *);
SILO_API extern void     DBFreePHZonelist(DBphzonelist *);
SILO_API extern void     DBFreeCSGZonelist(DBcsgzonelist *);
SILO_API extern void     DBResetUcdvar(DBucdvar *);
SILO_API extern void     DBFreeUcdvar(DBucdvar *);
SILO_API extern void     DBResetQuadvar(DBquadvar *);
SILO_API extern void     DBFreeCsgvar(DBcsgvar *);
SILO_API extern void     DBFreeQuadvar(DBquadvar *);
SILO_API extern void     DBFreeUcdmesh(DBucdmesh *);
SILO_API extern void     DBFreeMeshvar(DBmeshvar *);
SILO_API extern void     DBFreePointmesh(DBpointmesh *);
SILO_API extern void     DBFreeQuadmesh(DBquadmesh *);
SILO_API extern void     DBFreeCsgmesh(DBcsgmesh *);
SILO_API extern void     DBFreeDefvars(DBdefvars*);
SILO_API extern void     DBFreeMultimesh(DBmultimesh *);
SILO_API extern void     DBFreeMultimeshadj(DBmultimeshadj *);
SILO_API extern void     DBFreeMultivar(DBmultivar *);
SILO_API extern void     DBFreeMultimat(DBmultimat *);
SILO_API extern void     DBFreeMultimatspecies(DBmultimatspecies *);
SILO_API extern void     DBFreeCompoundarray(DBcompoundarray *);
SILO_API extern void     DBFreeCurve(DBcurve *);
SILO_API extern void     DBFreeNamescheme(DBnamescheme *);

SILO_API extern long     DBSetDataReadMask(long);
SILO_API extern long     DBGetDataReadMask(void);
SILO_API extern int      DBSetAllowOverwrites(int allow);
SILO_API extern int      DBGetAllowOverwrites(void);
SILO_API extern int      DBSetEnableChecksums(int enable);
SILO_API extern int      DBGetEnableChecksums(void);
SILO_API extern void     DBSetCompression(const char *);
SILO_API extern char    *DBGetCompression(void);
SILO_API extern int      DBSetFriendlyHDF5Names(int enable);
SILO_API extern int      DBGetFriendlyHDF5Names(void);
SILO_API extern int      DBGuessHasFriendlyHDF5Names(DBfile *f);
SILO_API extern int      DBSetDeprecateWarnings(int max);
SILO_API extern int      DBGetDeprecateWarnings();
SILO_API extern int     *DBSetUnknownDriverPriorities(const int *);
SILO_API extern int     *DBGetUnknownDriverPriorities();
SILO_API extern int      DBRegisterFileOptionsSet(const DBoptlist *opts);
SILO_API extern int      DBUnregisterFileOptionsSet(int opts_set_id);
SILO_API extern void     DBUnregisterAllFileOptionsSets();
SILO_API extern void    *DBGrabDriver(DBfile *);
SILO_API extern int      DBUngrabDriver(DBfile *, const void *);
SILO_API extern int      DBGetDriverType(const DBfile *);
SILO_API extern int      DBGetDriverTypeFromPath(const char *);
SILO_API extern char    *DBJoinPath(const char *, const char *);
SILO_API extern char    *DBVersion(void);
SILO_API extern int      DBVersionGE(int Maj, int Min, int Pat);
SILO_API extern char    *DBFileVersion(DBfile *dbfile);
SILO_API extern int      DBFileVersionGE(DBfile *dbfile, int Maj, int Min, int Pat);
SILO_API extern void     DBShowErrors(int, DBErrFunc_t);
SILO_API extern char    *DBErrString(void);
SILO_API extern char    *DBErrFunc(void);
SILO_API extern char    *DBErrFuncname(void);
SILO_API extern DBErrFunc_t DBErrfunc(void);
SILO_API extern int      DBErrno(void);
SILO_API extern int      DBErrlvl(void);
SILO_API extern int      DBClose(DBfile *);
SILO_API extern int      DBPause(DBfile *);
SILO_API extern int      DBContinue(DBfile *);
SILO_API extern int      DBInqVarExists(DBfile *, const char *);
SILO_API extern int      DBForceSingle(int);
SILO_API extern int      DBUninstall(DBfile *);
SILO_API extern DBoptlist *DBMakeOptlist(int);
SILO_API extern int      DBClearOptlist(DBoptlist *);
SILO_API extern int      DBFreeOptlist(DBoptlist *);
SILO_API extern int      DBAddOption(DBoptlist *, int, void *);
SILO_API extern void    *DBGetOption(const DBoptlist *, int);
SILO_API extern int      DBClearOption(DBoptlist *, int);
SILO_API extern DBtoc   *DBGetToc(DBfile *);
SILO_API extern int      DBNewToc(DBfile *);
SILO_API extern int      DBSortObjectsByOffset(DBfile *, int nobjs,
                             const char *const *const obj_names, int *ranks);
SILO_API extern int      DBFilters(DBfile *, FILE *);
SILO_API extern int      DBFilterRegistration(const char *, int (*init) (DBfile *, char *),
                                     int (*open) (DBfile *, char *));
SILO_API extern void    *DBGetAtt(DBfile *, const char *, const char *);
SILO_API extern DBobject *DBGetObject(DBfile *, const char *);
SILO_API extern int      DBChangeObject(DBfile *, DBobject *);
SILO_API extern int      DBWriteObject(DBfile *, DBobject *, int);
SILO_API extern void    *DBGetComponent(DBfile *, const char *, const char *);
SILO_API extern int      DBGetComponentType(DBfile *, const char *, const char *);
SILO_API extern int      DBWriteComponent(DBfile *, DBobject *, const char *, const char *, const char *,
                                 const void *, int, long *);
SILO_API extern int      DBWrite(DBfile *, const char *, void *, int *, int, int);
SILO_API extern int      DBWriteSlice(DBfile *, const char *, void *, int, int[], int[],
                             int[], int[], int);
SILO_API extern DBfacelist *DBCalcExternalFacelist(int *, int, int, int *, int *, int,
                                          int *, int);
SILO_API extern DBfacelist *DBCalcExternalFacelist2(int *, int, int, int, int, int *,
                                           int *, int *, int, int *, int);
SILO_API extern int      DBGetDir(DBfile *, char *);
SILO_API extern int      DBSetDir(DBfile *, const char *);
SILO_API extern int      DBSetDirID(DBfile *, int);
SILO_API extern int      DBListDir(DBfile *, char **, int);
SILO_API extern int      DBMkDir(DBfile *, const char *);
SILO_API extern int      DBCpDir(DBfile *dbfile, const char *srcDir,
                             DBfile *dstFile, const char *dstDir);

#define DBMkdir DBMkDir
SILO_API extern int      DBReadAtt(DBfile *, const char *, const char *, void *);
SILO_API extern int      DBRead(DBfile *, const char *, void *);
SILO_API extern int      DBReadVar(DBfile *, const char *, void *);
SILO_API extern int      DBReadVar1(DBfile *, const char *, int, void *);
SILO_API extern int      DBReadVarSlice(DBfile *, const char *, int *, int *, int *, int,
                               void *);
SILO_API extern DBobject *DBMakeObject(const char *, int, int);
SILO_API extern int      DBFreeObject(DBobject *);
SILO_API extern int      DBClearObject(DBobject *);
SILO_API extern int      DBAddVarComponent(DBobject *, const char *, const char *);
SILO_API extern int      DBAddIntComponent(DBobject *, const char *, int);
SILO_API extern int      DBAddFltComponent(DBobject *, const char *, double);
SILO_API extern int      DBAddDblComponent(DBobject *, const char *, double);
SILO_API extern int      DBAddStrComponent(DBobject *, const char *, const char *);
SILO_API extern int      DBGetComponentNames(DBfile *, const char *, char ***, char ***);

SILO_API extern DBcompoundarray *DBGetCompoundarray(DBfile *, const char *);
SILO_API extern DBcurve *DBGetCurve(DBfile *, const char *);
SILO_API extern DBdefvars *DBGetDefvars(DBfile *, const char *);
SILO_API extern DBmaterial *DBGetMaterial(DBfile *, const char *);
SILO_API extern DBmatspecies *DBGetMatspecies(DBfile *, const char *);
SILO_API extern DBmultimesh *DBGetMultimesh(DBfile *, const char *);
SILO_API extern DBmultimeshadj *DBGetMultimeshadj(DBfile *, const char *,
                                                  int, const int *);
SILO_API extern DBmultivar *DBGetMultivar(DBfile *, const char *);
SILO_API extern DBmultimat *DBGetMultimat(DBfile *, const char *);
SILO_API extern DBmultimatspecies *DBGetMultimatspecies(DBfile *, const char *);
SILO_API extern DBpointmesh *DBGetPointmesh(DBfile *, const char *);
SILO_API extern DBmeshvar *DBGetPointvar(DBfile *, const char *);
SILO_API extern DBquadmesh *DBGetQuadmesh(DBfile *, const char *);
SILO_API extern DBquadvar *DBGetQuadvar(DBfile *, const char *);
SILO_API extern int      DBGetQuadvar1(DBfile *, const char *, DB_DTPTR1, int *, int *,
                              DB_DTPTR1, int *, int *, int *);
SILO_API extern int      DBAnnotateUcdmesh(DBucdmesh *);
SILO_API extern DBucdmesh *DBGetUcdmesh(DBfile *, const char *);
SILO_API extern DBucdvar *DBGetUcdvar(DBfile *, const char *);
SILO_API extern DBcsgmesh *DBGetCsgmesh(DBfile *, const char *);
SILO_API extern DBcsgvar *DBGetCsgvar(DBfile *, const char *);
SILO_API extern DBcsgzonelist *DBGetCSGZonelist(DBfile *, const char *);
SILO_API extern DBfacelist *DBGetFacelist(DBfile *, const char *);
SILO_API extern DBzonelist *DBGetZonelist(DBfile *, const char *);
SILO_API extern DBphzonelist *DBGetPHZonelist(DBfile *, const char *);
SILO_API extern void    *DBGetVar(DBfile *, const char *);
SILO_API extern int      DBGetVarByteLength(DBfile *, const char *);
SILO_API extern int      DBGetVarLength(DBfile *, const char *);
SILO_API extern int      DBGetVarDims(DBfile *, const char *, int, int *);
SILO_API extern int      DBGetVarType(DBfile *, const char *);
SILO_API extern int      DBInqFileHasObjects(DBfile *);
SILO_API extern int      DBInqMeshname(DBfile *, const char *, const char *);
SILO_API extern int      DBInqMeshtype(DBfile *, const char *);
SILO_API extern int      DBInqCompoundarray(DBfile *, const char *, char ***,
                                   int **, int *, int *, int *);
SILO_API extern DBObjectType DBInqVarType(DBfile *, const char *);

SILO_API extern int      DBPutCompoundarray(DBfile *, const char *, char **, int *, int,
                                   void *, int, int, DBoptlist *);
SILO_API extern int      DBPutCurve(DBfile *, const char *, void *, void *, int, int,
                           DBoptlist *);
SILO_API extern int      DBPutDefvars (DBfile *, const char *, int, char **, 
                            const int *, char **, DBoptlist **);
SILO_API extern int      DBPutFacelist(DBfile *, const char *, int, int, int *, int, int,
                            int *, int *, int *, int, int *, int *, int);
SILO_API extern int      DBPutMaterial(DBfile *, const char *, const char *, int, int *, int *,
                           int *, int, int *, int *, int *, DB_DTPTR1, int,
                              int, DBoptlist *);
SILO_API extern int      DBPutMatspecies(struct DBfile *, const char *, const char *, int, int *,
                                int *, int *, int, int, DB_DTPTR1, int *, int, int,
                                DBoptlist *);
SILO_API extern int      DBPutMultimesh(DBfile *, const char *, int, char **, int *,
                               DBoptlist *);
SILO_API extern int      DBPutMultimeshadj(DBfile *, const char *, int, const int *,
                               const int *, const int *, const int *, const int *,
                               int **, const int *, int **, DBoptlist *optlist);
SILO_API extern int      DBPutMultivar(DBfile *, const char *, int, char **, int *,
                              DBoptlist *);
SILO_API extern int      DBPutMultimat(DBfile *, const char *, int, char **, DBoptlist *);
SILO_API extern int      DBPutMultimatspecies(DBfile *, const char *, int, char **,
                                     DBoptlist *);
SILO_API extern int      DBPutPointmesh(DBfile *, const char *, int, DB_DTPTR2, int, int,
                               DBoptlist *);
SILO_API extern int      DBPutPointvar(DBfile *, const char *, const char *, int, DB_DTPTR2, int,
                              int, DBoptlist *);
SILO_API extern int      DBPutPointvar1(DBfile *, const char *, const char *, DB_DTPTR1, int, int,
                               DBoptlist *);
SILO_API extern int      DBPutQuadmesh(DBfile *, const char *, char **, DB_DTPTR2, int *, int,
                              int, int, DBoptlist *);
SILO_API extern int      DBPutQuadvar(DBfile *, const char *, const char *, int, char **, DB_DTPTR2,
                             int *, int, DB_DTPTR2, int, int, int, DBoptlist *);
SILO_API extern int      DBPutQuadvar1(DBfile *, const char *, const char *, DB_DTPTR1, int *, int,
                              DB_DTPTR1, int, int, int, DBoptlist *);
SILO_API extern int      DBPutUcdmesh(DBfile *, const char *, int, char **, DB_DTPTR2, int,
                             int, const char *, const char *, int, DBoptlist *);
SILO_API extern int      DBPutUcdsubmesh(DBfile *, const char *, const char *, int,
                             const char *, const char *, DBoptlist *);
SILO_API extern int      DBPutUcdvar(DBfile *, const char *, const char *, int, char **, DB_DTPTR2,
                            int, DB_DTPTR2, int, int, int, DBoptlist *);
SILO_API extern int      DBPutUcdvar1(DBfile *, const char *, const char *, DB_DTPTR1, int, DB_DTPTR1,
                             int, int, int, DBoptlist *);
SILO_API extern int      DBPutZonelist(DBfile *, const char *, int, int, int *, int, int,
                              int *, int *, int);
SILO_API extern int      DBPutZonelist2(DBfile *, const char *, int, int, int *, int, int,
                               int, int, int *, int *, int *, int, DBoptlist*);
SILO_API extern int      DBPutPHZonelist(DBfile *, const char *, int, int *, int, int *, const char *,
                                                           int, int *, int, int *,
                                                           int, int, int, DBoptlist *);
SILO_API extern int      DBPutCsgmesh(DBfile *, const char *, int, int, const int *, const int *,
                                      const void *, int, int, const double *, const char *,
                                      DBoptlist *);
SILO_API extern int      DBPutCSGZonelist(DBfile *, const char *, int, const int *,
                                          const int *, const int *, const void *, int, int,
                                          int, const int *, DBoptlist *);
SILO_API extern int      DBPutCsgvar(DBfile *, const char *, const char *, int, char **,
                                     void **, int, int, int, DBoptlist *);

SILO_API extern void DBFreeMrgtree(DBmrgtree *tree);
SILO_API extern void DBPrintMrgtree(DBmrgtnode *tnode, int walk_order, void *data);
SILO_API extern void DBLinearizeMrgtree(DBmrgtnode *tnode, int walk_order, void *data);
SILO_API extern void DBWalkMrgtree(DBmrgtree *tree, DBmrgwalkcb cb, void *wdata,
                         int traversal_order);

SILO_API extern DBmrgtree *DBMakeMrgtree(int source_mesh_type, int mrgtree_info,
                               int max_root_descendents, DBoptlist *opts);
SILO_API extern int      DBAddRegion(DBmrgtree *tree, const char *region_name,
                             int type_info_bits, int max_descendents,
                             const char *maps_name, int nsegs, int *seg_ids,
                             int *seg_sizes, int *seg_types, DBoptlist *opts);

SILO_API extern int      DBAddRegionArray(DBmrgtree *tree, int nregn,
                             char **regn_names, int type_info_bits,
                             const char *maps_name, int nsegs, int *seg_ids,
                             int *seg_sizes, int *seg_types, DBoptlist *opts);

SILO_API extern int      DBSetCwr(DBmrgtree *tree, const char *path);
SILO_API extern const char *DBGetCwr(DBmrgtree *tree);

SILO_API extern int      DBPutMrgtree(DBfile *dbfile, const char *mrg_tree_name,
                             const char *mesh_name, DBmrgtree *tree, DBoptlist *opts);
SILO_API extern DBmrgtree *DBGetMrgtree(DBfile *dbfile, const char *mrg_tree_name);

SILO_API extern int      DBPutMrgvar(DBfile *dbfile, const char *name,
                             const char *mrgt_name,
			     int ncomps, char **compnames,
			     int nregns, char **reg_pnames,
                             int datatype, void **data, DBoptlist *opts);
SILO_API extern DBmrgvar *DBGetMrgvar(DBfile *dbfile, const char *name);
SILO_API extern void DBFreeMrgvar(DBmrgvar *mrgv);

SILO_API extern int      DBPutGroupelmap(DBfile *dbfile, const char *map_name,
                             int num_segments, int *groupel_types, int *segment_lengths,
                             int *segment_ids, int **segment_data, void **segment_fracs,
                             int fracs_data_type, DBoptlist *opts);
SILO_API extern DBgroupelmap *DBGetGroupelmap(DBfile *dbfile, const char *name);
SILO_API extern void DBFreeGroupelmap(DBgroupelmap *map);

SILO_API extern void *   DBFortranAccessPointer(int value);
SILO_API extern int      DBFortranAllocPointer(void *pointer);
SILO_API extern void     DBFortranRemovePointer(int value);

SILO_API extern int      DBVariableNameValid(const char *s);
SILO_API extern char *safe_strdup(const char *s);

SILO_API extern int      DBFreeCompressionResources(DBfile *dbfile, const char *meshname);

SILO_API extern DBnamescheme *DBMakeNamescheme(const char *fmt, ...);
SILO_API const char *DBGetName(DBnamescheme *ns, int natnum);

SILO_API extern void DBStringArrayToStringList(char **strArray, int n,
                         char **strList, int *m);
SILO_API char ** DBStringListToStringArray(char *strList, int n, int handleSlashSwap,
                     int skipFirstSemicolon);

/*-------------------------------------------------------------------------
 * Public global variables.
 *-------------------------------------------------------------------------
 */
SILO_API extern int     DBDebugAPI;      /*file desc for debug messages, or zero */
SILO_API extern int     db_errno;        /*error number of last error */
SILO_API extern char    db_errfunc[];    /*name of erring function */

#ifndef DB_MAIN
SILO_API extern DBfile *(*DBOpenCB[])(const char *, int, int);
SILO_API extern DBfile *(*DBCreateCB[])(const char *, int, int, int, const char *);
SILO_API extern int     (*DBFSingleCB[])(int);
#endif

#ifdef __cplusplus
}
#endif
#undef NO_FORTRAN_DEFINE
#endif /* !SILO_H */
