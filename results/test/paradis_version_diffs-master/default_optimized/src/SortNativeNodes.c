13a14
> #include "Precipitate.h"
142a144,292
> 
> 
> 
> 
> 
> void AssignPrecipitateToCell(Home_t *home, Precipitate_t *precipitate)
> {
>         Param_t *param;
>         int iCell, jCell, kCell, cellIdx;
>         real8 probXmin, probYmin, probZmin;
>         real8 cellXsize, cellYsize, cellZsize;
>         Cell_t *cell;
> 
>         if (precipitate == (Precipitate_t *)NULL) return;
> 
> /*
>  *      set the lower limit of the base area (excluding possible periodic cells)
>  *      and the size of each cell
>  */
>         param = home->param;
>         probXmin = param->minSideX;
>         probYmin = param->minSideY;
>         probZmin = param->minSideZ;
> 
>         cellXsize = (param->maxSideX - param->minSideX) / param->nXcells;
>         cellYsize = (param->maxSideY - param->minSideY) / param->nYcells;
>         cellZsize = (param->maxSideZ - param->minSideZ) / param->nZcells;
> 
> /*
>  *      put the precipitate on its proper cell. If the index exceeds this domains
>  *      range of native cells, put in the nearest native cell
>  */
>         iCell = (int)((precipitate->x - probXmin) / cellXsize);
>         if (iCell < param->iCellNatMin) iCell = param->iCellNatMin;
>         if (iCell > param->iCellNatMax) iCell = param->iCellNatMax;
>         iCell++;   /* compensate for periodic cells */
> 
>         jCell = (int)((precipitate->y - probYmin) / cellYsize);
>         if (jCell < param->jCellNatMin) jCell = param->jCellNatMin;
>         if (jCell > param->jCellNatMax) jCell = param->jCellNatMax;
>         jCell++;   /* compensate for periodic cells */
> 
>         kCell = (int)((precipitate->z - probZmin) / cellZsize);
>         if (kCell < param->kCellNatMin) kCell = param->kCellNatMin;
>         if (kCell > param->kCellNatMax) kCell = param->kCellNatMax;
>         kCell++;   /* compensate for periodic cells */
> 
>         cellIdx = EncodeCellIdx (home, iCell, jCell, kCell);
>         cell = home->cellKeys[cellIdx];
>         precipitate->nextInCell = cell->precipitateQ;
>         cell->precipitateQ = precipitate;
>         cell->precipitateCount++;
> 
>         precipitate->cellIdx = cellIdx;
> 
>         return;
> }
> 
> 
> 
> void SortNativePrecipitates (Home_t *home)
> {
>    Param_t *param;
>    int i, iCell, jCell, kCell, cellIdx ;
>    real8 probXmin, probYmin, probZmin ;
>    real8 cellXsize, cellYsize, cellZsize ;
>    Cell_t *cell ;
>    Precipitate_t *precipitate ;
> 
>    
> 
> /* set the lower limit of the base area (excluding possible periodic cells)
>  * and the size of each cell
>  */
> 
>    param = home->param;
>    probXmin = param->minSideX ;
>    probYmin = param->minSideY ;
>    probZmin = param->minSideZ ;
> 
>    cellXsize = (param->maxSideX - param->minSideX) / param->nXcells ;
>    cellYsize = (param->maxSideY - param->minSideY) / param->nYcells ;
>    cellZsize = (param->maxSideZ - param->minSideZ) / param->nZcells ;
> 
> /* Loop through cells and reinitialize their precipitate 
>    queues to empty.  */
> 
>    for (i = 0 ; i < home->cellCount ; i++) {
> 
>       cell = home->cellKeys[home->cellList[i]] ;
>       cell->precipitateQ = 0 ;
>       cell->precipitateCount = 0 ;
>    }
> 
> /* Loop thru active precipitates, putting them in their proper cell. If the
>  * index exceeds this domains range of native cells, put in the nearest 
>  * native cell
>  */
> 
>    for (i = 0 ; i < home->newPrecipitateKeyPtr ; i++) {
> 
>       precipitate = home->precipitateKeys[i] ;
>       if (!precipitate) continue ;
> 
>       iCell = (int)((precipitate->x - probXmin) / cellXsize) ;
>       if (iCell < param->iCellNatMin) iCell = param->iCellNatMin ;
>       if (iCell > param->iCellNatMax) iCell = param->iCellNatMax ;
>       iCell++ ;   /* compensate for periodic cells */
> 
>       jCell = (int)((precipitate->y - probYmin) / cellYsize) ;
>       if (jCell < param->jCellNatMin) jCell = param->jCellNatMin ;
>       if (jCell > param->jCellNatMax) jCell = param->jCellNatMax ;
>       jCell++ ;   /* compensate for periodic cells */
> 
>       kCell = (int)((precipitate->z - probZmin) / cellZsize) ;
>       if (kCell < param->kCellNatMin) kCell = param->kCellNatMin ;
>       if (kCell > param->kCellNatMax) kCell = param->kCellNatMax ;
>       kCell++ ;   /* compensate for periodic cells */
> 
>       cellIdx = EncodeCellIdx (home, iCell, jCell, kCell) ;
>       cell = home->cellKeys[cellIdx] ;
>       precipitate->nextInCell = cell->precipitateQ ;
>       cell->precipitateQ = precipitate ;
>       cell->precipitateCount++ ;
> 
>       precipitate->cellIdx = cellIdx ;
>       precipitate->native = 1 ;
> 
>    }
> 
>      
> 
>         return;
> }
> 
> 
> 
> 
> 
> 
> 
> 
> 
> 
> 
> 
> 
> 
> 
