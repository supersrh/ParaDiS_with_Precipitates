36a37
>         Precipitate_t      *precipitate;
38a40
>         PrecipitateBlock_t *precipitateBlock, *thisPrecipitateBlock;
86a89
>                     FreeNodePNbrs(node);
113a117,152
> 
> 		precipitateBlock = home->precipitateBlockQ;
> 
>         while (precipitateBlock) {
> 
>             precipitate = precipitateBlock->precipitates;
> 
>             for (i = 0; i < PRECIPITATE_BLOCK_COUNT; i++) {
> 
>                 
> 
>                 DESTROY_LOCK(&precipitate->precipitateLock);
> 
>                 precipitate++;
>             }
> 
>             thisPrecipitateBlock = precipitateBlock;
>             precipitateBlock = precipitateBlock->next;
> 
>             memset(thisPrecipitateBlock->precipitates, 0, PRECIPITATE_BLOCK_COUNT * sizeof(Precipitate_t));
>             free(thisPrecipitateBlock->precipitates);
>             free(thisPrecipitateBlock);
>         }
> 
>         home->precipitateBlockQ = 0;
> 
>         if(home->precipitateKeys) {
>             free(home->precipitateKeys);
>             home->precipitateKeys = NULL;
>         }
> 
> 
> 
> 
> 
> 
120a160,165
> 
> 		 if (home->recycledPrecipitateHeapSize > 0) {
>             free(home->recycledPrecipitateHeap);
>             home->recycledPrecipitateHeap = (int *)NULL;
>         }
> 
