11a12
> #include "Precipitate.h"
200a202,384
> 
> 
> 
> Precipitate_t *PopFreePrecipitateQ(Home_t *home)
> {
> 	int		i;
> 	PrecipitateBlock_t	*precipitateBlock;
> 	Precipitate_t		*currPrecipitate, *freePrecipitate;
> 
> 	if (home->freePrecipitateQ == 0) {
> 
> 		precipitateBlock = (PrecipitateBlock_t *) malloc(sizeof(PrecipitateBlock_t));
> 		precipitateBlock->next = home->precipitateBlockQ;
> 		home->precipitateBlockQ = precipitateBlock;
> 
> 		precipitateBlock->precipitates = (Precipitate_t *)calloc(1, PRECIPITATE_BLOCK_COUNT *
> 						    sizeof(Precipitate_t));
> 
> 		currPrecipitate = precipitateBlock->precipitates;
> 		home->freePrecipitateQ = currPrecipitate;
> 
> 		for (i = 1; i < PRECIPITATE_BLOCK_COUNT; i++) {
> 			INIT_LOCK(&currPrecipitate->precipitateLock);
> 			currPrecipitate->next = currPrecipitate + 1;
> 			currPrecipitate++;
> 		}
> 
> 		currPrecipitate->next = 0;    /* last free precipitate */
> 		INIT_LOCK(&currPrecipitate->precipitateLock);
> 		home->lastFreePrecipitate = currPrecipitate;
> 	}
> 
> /*
>  *	Dequeue the first free precipitate and return it to the caller
>  */
> 	freePrecipitate = home->freePrecipitateQ;
> 	home->freePrecipitateQ = freePrecipitate->next;
> 
> 	return(freePrecipitate);
> }
> 
> 
> /*--------------------------------------------------------------------------
>  *
>  *	Function:	PushFreePrecipitateQ
>  *	Description:	Push a precipitate onto the top of the free queue
>  *
>  *-------------------------------------------------------------------------*/
> void PushFreePrecipitateQ(Home_t *home, Precipitate_t *precipitate)
> {
> 	precipitate->next = home->freePrecipitateQ;
> 	home->freePrecipitateQ = precipitate;
> 
> 	return;
> }
> 
> 
> /*--------------------------------------------------------------------------
>  *
>  *	Function:	PushNativePrecipitateQ
>  *	Description:	Push a precipitate onto the top of the native queue
>  *
>  *-------------------------------------------------------------------------*/
> void PushNativePrecipitateQ(Home_t *home, Precipitate_t *precipitate)
> {
> 	precipitate->next = home->nativePrecipitateQ;
> 	home->nativePrecipitateQ = precipitate;
> 
> 	return;
> }
> 
> 
> /*--------------------------------------------------------------------------
>  *
>  *	Function:	PushGhostPrecipitateQ
>  *	Description:	Push a precipitate onto the top of the ghost queue
>  *
>  *-------------------------------------------------------------------------*/
> void PushGhostPrecipitateQ(Home_t *home, Precipitate_t *precipitate)
> {
> 	if (home->ghostPrecipitateQ == 0) home->lastGhostPrecipitate = precipitate;
> 	precipitate->next = home->ghostPrecipitateQ;
> 	home->ghostPrecipitateQ = precipitate;
> 
> 	return;
> }
> 
> 
> /*--------------------------------------------------------------------------
>  *
>  *	Function:	RecycleGhostPrecipitates
>  *	Description:	Put all the precipitates in the ghostPrecipitateQ back onto
>  *			the freePrecipitateQ and reset the ghostPrecipitateQ to empty
>  *
>  *-------------------------------------------------------------------------*/
> void RecycleGhostPrecipitates(Home_t *home)
> {
> 	if (home->ghostPrecipitateQ == NULL) return;  /* nothing to do */
> 
> 	if (home->freePrecipitateQ == NULL) {
> 		home->freePrecipitateQ = home->ghostPrecipitateQ;
> 	} else {
> 		home->lastFreePrecipitate->next = home->ghostPrecipitateQ;
> 	}
> 
> 	home->lastFreePrecipitate = home->lastGhostPrecipitate;
> 	home->lastGhostPrecipitate = NULL;
> 	home->ghostPrecipitateQ = NULL;
> 
> 	return;
> }
> 
> 
> /*--------------------------------------------------------------------------
>  *
>  *	Function:	RemovePrecipitateFromCellQ
>  *	Description:	Remove the specified precipitate from the queue for
>  *                      the cell containing the precipitate.  Typically only
>  *                      needed when deleting a precipitate during topological
>  *                      changes.
>  *
>  *-------------------------------------------------------------------------*/
> void RemovePrecipitateFromCellQ(Home_t *home, Precipitate_t *precipitate)
> {
>         Cell_t *cell;
>         Precipitate_t *prevPrecipitate, *nextPrecipitate;
> 
>         if (precipitate == (Precipitate_t *)NULL) return;
>         if (precipitate->cellIdx < 0) return;
> 
>         cell = home->cellKeys[precipitate->cellIdx];
> 
>         if (precipitate == cell->precipitateQ) {
>             cell->precipitateQ = precipitate->nextInCell;
>             cell->precipitateCount--;
>         } else {
>             nextPrecipitate = cell->precipitateQ;
>             while (nextPrecipitate != (Precipitate_t *)NULL) {
>                 if (nextPrecipitate == precipitate) {
>                     prevPrecipitate->nextInCell = precipitate->nextInCell;
>                     cell->precipitateCount--;
>                     return;
>                 }
>                 prevPrecipitate = nextPrecipitate;
>                 nextPrecipitate = nextPrecipitate->nextInCell;
>             }
>         }
> 
>         return;
> }
> 
> 
> /*--------------------------------------------------------------------------
>  *
>  *	Function:	RemovePrecipitateFromCell2Q
>  *	Description:	Remove the specified precipitate from the queue for
>  *                      the cell2 containing the precipitate.  Typically only
>  *                      needed when deleting a precipitate during topological
>  *                      changes.
>  *
>  *-------------------------------------------------------------------------*/
> void RemovePrecipitateFromCell2Q(Home_t *home, Precipitate_t *precipitate)
> {
> 
> 	if (home->cell2QentArray == (C2Qent_t *)NULL) {
> 		return;
> 	}
> 
>         if (precipitate == (Precipitate_t *)NULL) {
> 		return;
> 	}
> 
>         if ((precipitate->cell2Idx < 0) || (precipitate->cell2QentIdx < 0)) {
> 		return;
> 	}
> 
>         home->cell2QentArray[precipitate->cell2QentIdx].precipitate = (Precipitate_t *)NULL;
>         home->cell2QentArray[precipitate->cell2QentIdx].next = -1;
> 
>         return;
> }
> 
> 
