13a14
> #include "InPrecipitateData.h"
29a31
>             FreeNodePNbrs(&inData->node[i]);
85a88,161
> 
> void FreeInPrecipitateArray(InPrecipitateData_t *inprecipitateData, int numPrecipitates)
> {
>         int i;
> 
>         if (inprecipitateData->precipitate == (Precipitate_t *)NULL) {
>             return;
>         }
> 
> /*
>  *      For each precipitate in inprecipitateData->precipitate, free the precipitate's arm arrays
>  */
>         //for (i = 0; i < numPrecipitates; i++) {
>             //FreePrecipitateArms(&inprecipitateData->precipitate[i]);
>         //}
> 
>         free(inprecipitateData->precipitate);
> 
>         inprecipitateData->precipitate      = (Precipitate_t *)NULL;
>         inprecipitateData->precipitateCount = 0;
> 
>         return;
> }
> 
> 
> void FreeInitPrecipitateArrays(Home_t *home, InPrecipitateData_t *inprecipitateData)
> {
> /*
>  *      Only the domains actively involved in reading nodal data
>  *      during initialization will have allocated InData_t
>  *      structures, so don't try to free structures that have
>  *      not been allocated.
>  */
>         if (inprecipitateData == (InPrecipitateData_t *)NULL) {
>             return;
>         }
> 
> /*
>  *      Free the inprecipitateData precipitate array and all arrays associated
>  *      with each precipitate's arms.
>  */
>         FreeInPrecipitateArray(inprecipitateData, inprecipitateData->precipitateCount);
> 
> /*
>  *      Free all remaining arrays in the inprecipitateData structure
>  */
>         //if (inprecipitateData->burgX != (real8 *)NULL) {
>             //free(inprecipitateData->burgX);
>             //inprecipitateData->burgX = (real8 *)NULL;
>         //}
>         //if (inprecipitateData->burgY != (real8 *)NULL) {
>             //free(inprecipitateData->burgY);
>             //inprecipitateData->burgY = (real8 *)NULL;
>         //}
>         //if (inprecipitateData->burgZ != (real8 *)NULL) {
>             //free(inprecipitateData->burgZ);
>             //inprecipitateData->burgZ = (real8 *)NULL;
>         //}
> 
> /*
>  *      Memory associated with the domain decomposition is dependent
>  *      on the type of decomposition used, so invoke a generic function
>  *      that will take the appropriate actions based on the decomposition
>  *      type.
>  */
>         //FreeDecomp(home, inprecipitateData->decomp);
>         //inprecipitateData->decomp = (void *)NULL;
> 
>         return;
> }
> 
> 
> 
> 
