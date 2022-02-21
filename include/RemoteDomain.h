/****************************************************************************
 *
 *	RemoteDomain.h	Define the struct that holds various data used in 
 *			communicating with a neighboring domain
 *
 ****************************************************************************/

#ifndef _RemoteDomain_h
#define  _RemoteDomain_h

#include "Typedefs.h"
#include "Tag.h"

struct _remotedomain {
	int	domainIdx;	/* encoded index of this domain */
	int	numExpCells;	/* number of native cells exported to */
				/* this domain                        */
	int	*expCells;	/* list of encoded indices of the */
				/* exported cells                 */

	int	maxTagIndex;	/* error if tag.index exceeds this value */
	Node_t	**nodeKeys;	/* indexed by node's tag.index, points */
				/* to Node_t  
				 *
				                          */
	int	maxprecipitateTagIndex;			                          
	Precipitate_t	**precipitateKeys;
	/* indexed by precipitates's tag.index, points */
				/* to Precipitate_t */

	int	inBufLen;
	int	inPBufLen;
	char	*inBuf;
	char	*inPBuf;
	int	outBufLen;
	int	outPBufLen;
	char	*outBuf;
	char	*outPBuf;
};

#endif
