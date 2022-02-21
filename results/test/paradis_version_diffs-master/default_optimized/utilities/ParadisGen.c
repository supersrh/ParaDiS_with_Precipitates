86a87
> #include "InPrecipitateData.h"
105c106,110
<         {FTYPE_EDGE,             FNAME_EDGE}
---
>         {FTYPE_EDGE,             FNAME_EDGE},
>         {FTYPE_PRECIPITATES,     FNAME_PRECIPITATES},
>         {FTYPE_FCC2D,     FNAME_FCC2D},
>         {FTYPE_SEDGE,     FNAME_SEDGE},
>         {FTYPE_D2EDGE,     FNAME_D2EDGE}
120c125
<         {OPT_FRLEN,     "frlen",    1, 1},
---
>         {OPT_FRLEN,     "frlen",    4, 1},
132c137
<         {OPT_TYPE,      "type",     1, 1},
---
>         {OPT_TYPE,      "type",     4, 1},
136c141,145
<         {OPT_ZSURF,     "zsurf",    1, 1}
---
>         {OPT_ZSURF,     "zsurf",    1, 1},
>         {OPT_FORCEP,  	"forcep", 	4, 1},
>         {OPT_NPRECIPITATES,   "nprecipitates",  1, 1},
>         {OPT_NODEZERO,  "nodezero", 1, 1}
>         
329d337
< 
335a344,347
>         inArgs->forcep	   = 1.0;
>         inArgs->numPrecipitates  = 1;
>         inArgs->nodezero   = 0;
>         	
412c424,425
<         printf("Maximum Segment Length    %lf\n",inArgs->maxSegLen);
---
>         printf("Maximum Segment Length:    %lf\n",inArgs->maxSegLen);
>         printf("Zeronode:    %d\n",inArgs->nodezero);
447a461
>             case FTYPE_FCC2D:
457a472,482
>                 
>             case FTYPE_PRECIPITATES:
>                 printf("Number of precipitates:         %d\n", inArgs->numPrecipitates);
>                 printf(" Forcep:          %e\n", inArgs->forcep);
>                 printf("Radius:  %e\n", inArgs->radius);
>                 printf("Seed Value:               %d\n", inArgs->seed);
>                 break;    
>                 
>                 
>                 
>                 
500a526
>            
510a537
> 						
524a552
> 				 
574a603,604
>                  
>                     
587a618
>           
590a622,625
>                
>                 
>                
>                     
648a684,705
>                     
>                     
>                   
>                  case OPT_FORCEP:
> 	
> 					inArgs->forcep = atof(argValue);
>                  break;  
> 
>                  case OPT_NPRECIPITATES:
>                
>                     inArgs->numPrecipitates = atoi(argValue);
>                     break;     
>                     
>                     
>                  case OPT_NODEZERO:
>                     inArgs->nodezero = atoi(argValue);
>                     break; 
>                     
>                      
>                     
>                           
>                     
675c732
< void WriteInitialNodeData(Home_t *home, InData_t *inData, int lastBlock)
---
> void WriteInitialNodeData(Home_t *home, InData_t *inData, int lastBlock, int nodezero)
691c748
<             fp = fopen(param->node_data_file, "w");
---
>             fp = fopen(param->node_data_file, "a");
729c786
<                     node->myTag.domainID, node->myTag.index, 
---
>                     node->myTag.domainID, node->myTag.index+nodezero, 
737c794
<                         node->nbrTag[iArm].index,
---
>                         node->nbrTag[iArm].index+nodezero,
758a816,920
> /*---------------------------------------------------------------------------
>  *
>  *      Function:       WriteInitialPrecipitateData
>  *      Description:    Writes the precipitate data for all precipitates contained
>  *                      in the inprecipitateData->precipitate list to the specified data
>  *                      file.  The first call to ths function will
>  *                      open up the output file, write the version
>  *                      and header data and so on before writing precipitate
>  *                      data.  The final call (i.e. lastBlock == 1) will
>  *                      insert the final precipitate count into the data file
>  *                      and clean up.
>  *      Arguments:
>  *
>  *          lastBlock       set to 1 if if the supplied precipitate data
>  *                          is the last block of dfata to be written
>  *                          to the specified file, set to zero in
>  *                          all other cases.
>  *
>  *-------------------------------------------------------------------------*/
> void WriteInitialPrecipitateData(Home_t *home, InPrecipitateData_t *inprecipitateData, int lastBlock, int precipitatezero)
> {
>         int         i, iArm;
>         Precipitate_t      *precipitate;
>         Param_t     *param;
>         static int  totalPrecipitateCount = 0;
>         static FILE *fp = (FILE *)NULL;
> 
> 
>         param = inprecipitateData->param;
> /*
>  *      If this is the first block of data being written to
>  *      this file (i.e. fp == NULL), do some basic initialization
>  */
>         if (fp == (FILE *)NULL) {
> 			  
>             fp = fopen(param->precipitate_data_file, "a");
>           
>             if (!fp) {
>                 Fatal("%s: error %d opening file %s\n",
>                       "WriteInitialNodeData", errno,
>                       param->precipitate_data_file);
>             }
>              
> 
> /*
>  *          Write data file parameters
>  */
> 
>             WriteParam(home->precipitateParamList, -1, fp);
> 
> /*
>  *          Write the domain decomposition into the data file
>  */
>             fprintf(fp, "\n#\n#  END OF DATA FILE PARAMETERS\n#\n\n");
> 
>             fprintf(fp, "domainDecomposition = \n");
>             WriteDecompBounds(home, fp);
> 
>             fprintf(fp, "precipitateData = \n");
>             fprintf(fp, "#  Primary lines: precipitate_tag, x, y, z, "
>                     "forcep, r , status (0=active 1=dead)  \n");
>             
>         }
> 
> 
> 
> /*
>  *      Dump all the precipitate data in this block to the file
>  */
>         totalPrecipitateCount += inprecipitateData->precipitateCount;
> 
>         for (i = 0; i < inprecipitateData->precipitateCount; i++) {
> 			
>             precipitate = &inprecipitateData->precipitate[i];
>             
>             
>            
> 
>             fprintf(fp,
>                     " %d,%d %.4f %.4f %.4f %.4f %.4f %d\n",
>                     precipitate->myTag.domainID, precipitate->myTag.index+precipitatezero, 
>                     precipitate->x, precipitate->y, precipitate->z, precipitate->forcep,
>                     precipitate->r,precipitate->constraint);
> 
> 
>         }
> 
> 
> 
> /*
>  *      If this is the last block of data being written to the
>  *      file, seek back to the node count in the file, overwrite
>  *      it with the final total, and cleanup.
>  */
>         if (lastBlock) {
>             fclose(fp);
>             fp = (FILE *)NULL;
>             printf("\nTotal precipitate count:         %d\n", totalPrecipitateCount);
>         }
> 
>         return;
> }
> 
> 
> 
767a930
>         InPrecipitateData_t        inprecipitateData;
775a939
>         memset(&inprecipitateData, 0, sizeof(InData_t));
777a942
>         inprecipitateData.param = (Param_t *)calloc(1, sizeof(Param_t));
778a944,945
>         home.precipitateParamList = (ParamList_t *)calloc(1, sizeof(ParamList_t));
>          
787c954
< 
---
> 		PrecipitateParamInit(inprecipitateData.param, home.precipitateParamList);
816a984,985
>                 
>                 
867a1037,1038
>  
>   
872c1043
<                                   &totDislocLen, inArgs.type);
---
>                                   &totDislocLen, inArgs.type,inArgs.nodezero);
877c1048,1104
<                             &totDislocLen, inArgs.type);
---
>                             &totDislocLen, inArgs.type,inArgs.nodezero);
>                 break; 
>                 
>         case FTYPE_SEDGE:
>                 CreateSingleEdges(&home, &inData, inArgs.cubeLength,
>                             inArgs.numChains, inArgs.seed,
>                             &totDislocLen, inArgs.type,inArgs.nodezero);
>                 break;        
>                 
>        case FTYPE_D2EDGE:
>                 CreateD2Edges(&home, &inData, inArgs.cubeLength,
>                             inArgs.numChains, inArgs.seed,
>                             &totDislocLen, inArgs.type,inArgs.nodezero);
>                 break;           
>                         
>         case FTYPE_PRECIPITATES:
> 				param = inprecipitateData.param;
> 				maxSide = inArgs.cubeLength / 2;
> 				minSide = -maxSide;
> 
> 				param->minSideX = minSide;
> 				param->minSideY = minSide;
> 				param->minSideZ = minSide;
> 
> 				param->maxSideX = maxSide;
> 				param->maxSideY = maxSide;
> 				param->maxSideZ = maxSide;
> 
> 				param->Lx = param->maxSideX - param->minSideX;
> 				param->Ly = param->maxSideY - param->minSideY;
> 				param->Lz = param->maxSideZ - param->minSideZ;
> 
> 				param->invLx = 1.0 / param->Lx;
> 				param->invLy = 1.0 / param->Ly;
> 				param->invLz = 1.0 / param->Lz;
> 
> 				param->maxSeg = inArgs.maxSegLen;
> 				
> 				param->minCoordinates[X] = minSide;
> 				param->minCoordinates[Y] = minSide;
> 				param->minCoordinates[Z] = minSide;
> 
> 				param->maxCoordinates[X] = maxSide;
> 				param->maxCoordinates[Y] = maxSide;
> 				param->maxCoordinates[Z] = maxSide;
> 				
> 				
> 				
> 				strncpy(param->precipitate_data_file, inArgs.outputFile,
>                 sizeof(param->precipitate_data_file)-1);
>                 
>                 printf("forcep %e \n",inArgs.forcep); 
>                 
> 				CreatePrecipitates(&home, &inprecipitateData, inArgs.cubeLength,
>                  inArgs.numPrecipitates, inArgs.seed, inArgs.radius,inArgs.forcep
>                  ,inArgs.nodezero);			
>                             
883c1110
<  *              provided boundary conditions here.  If no surface coordinates
---
>  *              provided boundary conditions here.  If ,no surface coordinates
917c1144
<                                         &totDislocLen, inArgs.type);
---
>                                         &totDislocLen, inArgs.type,inArgs.nodezero);
924c1151
<                                     inArgs.seed, &totDislocLen, inArgs.type);
---
>                                     inArgs.seed, &totDislocLen, inArgs.type,inArgs.nodezero);
930c1157
<                                &totDislocLen, inArgs.type);
---
>                                &totDislocLen, inArgs.type,inArgs.nodezero);
935c1162
<                                 inArgs.type);
---
>                                 inArgs.type,inArgs.nodezero);
936a1164,1171
>                 
>         case FTYPE_FCC2D:
>                 CreateFCC2DConfig(&home, &inData, inArgs.cubeLength,
>                                 inArgs.numChains, inArgs.seed, &totDislocLen,
>                                 inArgs.type,inArgs.nodezero);
>                 break;       
>                 
>                 
941c1176
<                                      &totDislocLen, inArgs.type);
---
>                                      &totDislocLen, inArgs.type,inArgs.nodezero);
