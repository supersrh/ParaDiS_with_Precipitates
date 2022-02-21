41a42
> #include "Precipitate.h"
173a175
>         ParamList_t *precipitateParamList;
189a192,204
>         
>         Precipitate_t    *nativePrecipitateQ;
>         Precipitate_t    *ghostPrecipitateQ;
>         Precipitate_t     *freePrecipitateQ;
> 
>         Precipitate_t     *lastFreePrecipitate;
>         Precipitate_t     *lastGhostPrecipitate;
> 
>         PrecipitateBlock_t  *precipitateBlockQ;
>         
>         
>         
>         
203a219,223
> 		Precipitate_t    **precipitateKeys;
>         int       newPrecipitateKeyPtr;
>         int       newPrecipitateKeyMax;
> 
> 
206a227,230
>         
>         int       *recycledPrecipitateHeap;
>         int       recycledPrecipitateHeapSize;
>         int       recycledPrecipitateHeapEnts;
286a311,315
>         
>         MPI_Request  *inPRequests; /* used for asynchronous I/O */
>         MPI_Request  *outPRequests;
>         MPI_Status   *inPStatus;
>         MPI_Status   *outPStatus;
348a378,384
>         
>         TagMap_t  *precipitatetagMap;
>         int       precipitatetagMapSize;
>         int       precipitatetagMapEnts;
>         
>         
>         
