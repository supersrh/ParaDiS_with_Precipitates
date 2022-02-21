17a18
>         Precipitate_t *precipitate;
72a74
>             free(remDom->precipitateKeys);
83a86,90
>         
>         free(home->inPRequests);
>         free(home->outPRequests);
>         free(home->inPStatus);
>         free(home->outPStatus);
112a120,143
>         
>         for (i = 0; i < home->newPrecipitateKeyPtr; i++) {
> 
>             if ((precipitate = home->precipitateKeys[i]) == (Precipitate_t *)NULL) {
>                 continue;
>             }
> 
>             precipitate->cellIdx      = -1;
>             precipitate->cell2Idx     = -1;
>             precipitate->cell2QentIdx = -1;
>             precipitate->nextInCell   = (Precipitate_t *)NULL;
>         }
> 
>         precipitate = home->ghostPrecipitateQ;
> 
>         while (precipitate) {
>             precipitate->cellIdx      = -1;
>             precipitate->cell2Idx     = -1;
>             precipitate->cell2QentIdx = -1;
>             precipitate->nextInCell   = (Precipitate_t *)NULL;
>             precipitate = precipitate->next;
>         }
>         
>         
