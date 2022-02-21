27a28
> #include "InPrecipitateData.h"
33a35
> #include "Comm.h"
690a693,700
>             
>          } else if (strcmp(param->mobilityLaw, "BCC_Fe_nl") == 0) {
>             param->materialType = MAT_TYPE_BCC;
>             param->mobilityType = MOB_BCC_FE_NL;
>             param->mobilityFunc = Mobility_BCC_Fe_nl;
>             param->numBurgGroups = 5;   
>             
>             
1109,1111c1119,1121
<         char         *sep, *start;
<         char         tmpDataFile[256], testFile[256];
<         char         *ctrlFile, *dataFile;
---
>         char         *sep, *start,*sep2,*start2;
>         char         tmpDataFile[256], tmpDataFile2[256], testFile[256];
>         char         *ctrlFile, *dataFile,*pdataFile;
1113a1124
>         InPrecipitateData_t *inprecipitateData;
1124a1136,1139
>         
>         home->precipitatetagMap     = (TagMap_t *)NULL;
>         home->precipitatetagMapSize = 0;
>         home->precipitatetagMapEnts = 0;
1132a1148
>         inprecipitateData = (InPrecipitateData_t *) calloc(1, sizeof(InPrecipitateData_t));
1192a1209
>                 strcpy(tmpDataFile2, ctrlFile);
1193a1211
>                 start2 = strrchr(tmpDataFile2, '/');
1196a1215,1220
>                 
>                 if (start2 == (char *)NULL) start2 = tmpDataFile2;
>                 sep2 = strrchr(start2, '.');
>                 if ((sep2 != (char *)NULL) && (sep2 > start2)) *sep2 = 0;
>                 
>                 
1219a1244
>                         strcat(tmpDataFile2, PRECIPITATEDATA_FILE_SUFFIX);
1225c1250
< 
---
> 				pdataFile = tmpDataFile2;
1241a1267
>             home->precipitateParamList = (ParamList_t *)calloc(1,sizeof(ParamList_t));
1248a1275,1276
>             PrecipitateParamInit(param, home->precipitateParamList);
>             
1255c1283,1285
< 
---
>   
>     
> 	
1359a1390
> 			ReadPrecipitateDataFile(home,inprecipitateData,inData,pdataFile);
1360a1392,1395
> #ifdef PARALLEL
> 			MPI_Barrier(MPI_COMM_WORLD);
> #endif
>             
1362a1398,1399
> 
> 
1415a1453
>  * 
1423a1462,1466
>         
>         //Tassa joku vika
>         
>         FreeInitPrecipitateArrays(home, inprecipitateData);
>         free(inprecipitateData);
1453a1497
> 		DistributePrecipitateTagMaps(home);		
1454a1499
>        
1526a1572,1577
>         SortNativePrecipitates(home);
>         
>         
>          //printf("Task %d beginning of Initiliaze flag 1, cycle %d\n", home->myDomain,home->cycle);
>         
> 		CommSendPrecipitateGhosts(home); 
1527a1579,1585
>         
>         
> 		 
>         //printf("Tassa2 %d %d \n",home->myDomain,home->cycle);
>         #ifdef PARALLEL
>         MPI_Barrier(MPI_COMM_WORLD);
> #endif
