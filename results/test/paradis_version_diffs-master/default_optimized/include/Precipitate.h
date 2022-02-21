/*--------------------------------------------------------------------------
 *
 *	Precipitate.h	Define the struct that holds all relevant data for a single
 *		precipitate, either native or ghost
 *
 *		Notice: make depend if new files are added, otherwise always
 *			make clean after changing this file  Wei Cai 04/09/2002
 *
 *------------------------------------------------------------------------*/

#ifndef _Precipitate_h
#define _Precipitate_h

#include "Typedefs.h"
#include "Tag.h"
#include "ParadisThread.h"

/*
 *      Define the various precipitate 'constraints' available.  Note: these
 *      constraints are mutually exclusive.
 */
#define UNCONSTRAINED  0
#define SURFACE_PRECIPITATE   1
#define PINNED_PRECIPITATE    7

/*
 *      Define the bit flags that can be set for the precipitate.  Note: these
 *      flags may be OR'ed together.
 */
#define PRECIPITATE_RESET_FORCES    0x01
#define PRECIPITATE_OUTSIDE_SURF    0x02
#define NO_COLLISIONS        0x04
#define NO_MESH_COARSEN      0x08
#define PRECIPITATE_CHK_DBL_LINK    0x10

/*
 *      Used as bit flags to indicate the type of nodal data items
 *      preserved by calls to PreserveNodalData().
 */
#define PRECIPITATE_POSITION 0x01
#define PRECIPITATE_CURR_VEL 0x02
#define PRECIPITATE_OLD_VEL  0x04

struct _precipitate {
	int	flags;
        
	real8	x, y, z;		/* nodal position */
	real8	fX, fY, fZ;		/* nodal force: units=Pa*b^2) */
	real8	vX, vY, vZ;		/* nodal velocity: units=burgMag/sec */
   real8 forcep;			/*precipitate force parameter*/	
   real8 r;					/*precipitate radius*/	
        
	real8	oldx, oldy, oldz;	/* for strain increment, Moono.Rhee */
	real8	oldvX, oldvY, oldvZ;	/* nodal velocity at previous step */
    real8   currvX, currvY, currvZ; /* nodal velocity at beginning of the */
                                        /* current step.  Only used during    */
                                        /* timestep integration while vX, vY  */
                                        /* and vZ are being determined.       */

	Tag_t	myTag;

/*
 *	nbrTag and burgID are dynamically allocated to size numNbrs, when the
 *	precipitate is allocated.
 */
	int	numNbrs;
	Tag_t	*nbrTag;

/*
 *	Arm information
 */
	real8	*armfx, *armfy, *armfz;	/* arm specific force contribution */
	real8	*burgX, *burgY, *burgZ;	/* burgers vector */        
	real8	*nx, *ny, *nz;		/* glide plane */        

	real8	*sigbLoc;		/* sig.b on arms (numNbr*3) */
	real8	*sigbRem;		/* sig.b on arms (numNbr*3) */

	int	*armCoordIndex;		/* Array of indices (1 per arm) into */
					/* the mirror domain's arrays of     */
					/* coordinates (armX, armY, armZ)    */
					/* of precipitates' neighbors.  This array  */
					/* is only present and useful while  */
					/* task zero is downloading the data */
					/* from the remote domain for        */
					/* generating output                 */

	int	constraint;     /* constraint =  1 : any surface precipitate   */
				/* constraint =  7 : Frank-Read end	*/
				/*                   points, fixed	*/
				/* constraint =  9 : Jog		*/
				/* constraint = 10 : Junction		*/
				/* constraint = 20 : cross slip(default)*/
				/*                   2-precipitates non-deletable*/

	int	cellIdx;	/* cell precipitate is currently sorted into */
	int	cell2Idx;	/* cell2 precipitate is currently sorted into */
	int	cell2QentIdx;	/* Index of this precipitate's entry in the */
				/* home->cell2QentArray.             */

	int	native;		/* 1 = native precipitate, 0 = ghost precipitate */

	Precipitate_t	*next;		/* pointer to the next precipitate in the queue */
				/* (ghost or free)			 */

	Precipitate_t	*nextInCell;	/* used to queue precipitate onto the current */
				/* containing cell		       */

	int	sgnv;		/* +1: if contribute to strain rate,	*/
				/* -1: if moving in opposite direction	*/

#ifdef DEBUG_LOG_MULTI_PRECIPITATE_SPLITS
        int     multiPrecipitateLife;
#endif

#ifdef _FEM
        int     fem_Surface[2];
        real8   fem_Surface_Norm[3];
#endif

#ifdef _OPENMP
        omp_lock_t precipitateLock;
#endif
};

struct _precipitateblock {
	PrecipitateBlock_t	*next;
	Precipitate_t		*precipitates;

};

#endif

