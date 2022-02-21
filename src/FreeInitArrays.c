/************************************************************************
 *
 *      Module:      FreeInitArrays
 *      Description: Contains functions needed to release temporary 
 *                   most memory allocated for the inData struct.
 *
 *      Includes functions:
 *          FreeInNodeArray()
 *          FreeInitArrays()
 *
 ************************************************************************/
#include "Home.h"
#include "InData.h"
#include "InPrecipitateData.h"
#include "Decomp.h"


void FreeInNodeArray(InData_t *inData, int numNodes)
{
        int i;

        if (inData->node == (Node_t *)NULL) {
            return;
        }

/*
 *      For each node in inData->node, free the node's arm arrays
 */
        for (i = 0; i < numNodes; i++) {
            FreeNodeArms(&inData->node[i]);
			FreeNodePNbrs(&inData->node[i]);
        }

        free(inData->node);

        inData->node      = (Node_t *)NULL;
        inData->nodeCount = 0;

        return;
}


void FreeInitArrays(Home_t *home, InData_t *inData)
{
/*
 *      Only the domains actively involved in reading nodal data
 *      during initialization will have allocated InData_t
 *      structures, so don't try to free structures that have
 *      not been allocated.
 */
        if (inData == (InData_t *)NULL) {
            return;
        }

/*
 *      Free the inData node array and all arrays associated
 *      with each node's arms.
 */
        FreeInNodeArray(inData, inData->nodeCount);

/*
 *      Free all remaining arrays in the inData structure
 */
        if (inData->burgX != (real8 *)NULL) {
            free(inData->burgX);
            inData->burgX = (real8 *)NULL;
        }
        if (inData->burgY != (real8 *)NULL) {
            free(inData->burgY);
            inData->burgY = (real8 *)NULL;
        }
        if (inData->burgZ != (real8 *)NULL) {
            free(inData->burgZ);
            inData->burgZ = (real8 *)NULL;
        }

/*
 *      Memory associated with the domain decomposition is dependent
 *      on the type of decomposition used, so invoke a generic function
 *      that will take the appropriate actions based on the decomposition
 *      type.
 */
        FreeDecomp(home, inData->decomp);
        inData->decomp = (void *)NULL;

        return;
}
 
 void FreeInPrecipitateArray(InPrecipitateData_t *inprecipitateData, int numPrecipitates)
 {
         int i;
 
         if (inprecipitateData->precipitate == (Precipitate_t *)NULL) {
             return;
         }
 
 /*
  *      For each precipitate in inprecipitateData->precipitate, free the precipitate's arm arrays
  */
         //for (i = 0; i < numPrecipitates; i++) {
             //FreePrecipitateArms(&inprecipitateData->precipitate[i]);
         //}
 
         free(inprecipitateData->precipitate);
 
         inprecipitateData->precipitate      = (Precipitate_t *)NULL;
         inprecipitateData->precipitateCount = 0;
 
         return;
 }
 
 
 void FreeInitPrecipitateArrays(Home_t *home, InPrecipitateData_t *inprecipitateData)
 {
 /*
  *      Only the domains actively involved in reading nodal data
  *      during initialization will have allocated InData_t
  *      structures, so don't try to free structures that have
  *      not been allocated.
  */
         if (inprecipitateData == (InPrecipitateData_t *)NULL) {
             return;
         }
 
 /*
  *      Free the inprecipitateData precipitate array and all arrays associated
  *      with each precipitate's arms.
  */
         FreeInPrecipitateArray(inprecipitateData, inprecipitateData->precipitateCount);
 
 /*
  *      Free all remaining arrays in the inprecipitateData structure
  */
         //if (inprecipitateData->burgX != (real8 *)NULL) {
             //free(inprecipitateData->burgX);
             //inprecipitateData->burgX = (real8 *)NULL;
         //}
         //if (inprecipitateData->burgY != (real8 *)NULL) {
             //free(inprecipitateData->burgY);
             //inprecipitateData->burgY = (real8 *)NULL;
         //}
         //if (inprecipitateData->burgZ != (real8 *)NULL) {
             //free(inprecipitateData->burgZ);
             //inprecipitateData->burgZ = (real8 *)NULL;
         //}
 
 /*
  *      Memory associated with the domain decomposition is dependent
  *      on the type of decomposition used, so invoke a generic function
  *      that will take the appropriate actions based on the decomposition
  *      type.
  */
         //FreeDecomp(home, inprecipitateData->decomp);
         //inprecipitateData->decomp = (void *)NULL;
 
         return;
 }
