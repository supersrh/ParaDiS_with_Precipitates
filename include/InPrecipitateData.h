/*****************************************************************************
 *
 *   InpData.h  Define the data structures used to hold the input data
 *
 ****************************************************************************/
#ifndef _InPrecipitateData_h
#define _InPrecipitateData_h

#include "Typedefs.h"
#include "Home.h"
#include "Tag.h"

/*
 *	During the initial problem generation, and parallel application
 *	initialization, we want to limit the number of precipitates for which
 *	data will be held in memory at any given point.  When that
 *	number exceeds the below value, the current nodal data will be
 * 	written to disk (problem generation) or distributed to remote
 *	domains (parallel execution).
 */
#define	MAX_PRECIPITATES_PER_BLOCK	50000

struct _inprecipitatedata {
	Param_t	*param;

	Precipitate_t  *precipitate;		/* array of precipitate structs */
	int	precipitateCount;	/* number of precipitates in <precipitate> array */
	int precipitatezero;	/* Beginning index for quick and dirty precipitate appending (AL) */
	real8	*burgX;		/* Array of XYZ burgers vector values */
	real8	*burgY;		/* which are now only used to support */
	real8	*burgZ;		/* old format ctrl files which include*/
				/* burger's vector arrays	      */

	int	nburg;		/* number of burgers vectors in the */
                                /* burgX, burgY and burgZ arrays    */

        void    *decomp;        /* pointer to memory allocated to hold   */
                                /* the domain decomposition read in from */
                                /* the restart file(s)                   */
};
#endif

