180a181,187
>         
>         
>         BindVar(CPList, "stopTime", &param->stopTime, V_DBL, 1, VFLAG_NULL);
>         param->stopTime = 1.0e-06;
>         
>         BindVar(CPList, "relaxTime", &param->relaxTime, V_DBL, 1, VFLAG_NULL);
>         param->relaxTime = 1.0e-06;
288a296,320
>                 
>         BindVar(CPList, "stressAmp", &param->stressAmp, V_DBL, 1,
>                 VFLAG_NULL); 
>                 
>         BindVar(CPList, "stressFreq", &param->stressFreq, V_DBL, 1,
>                 VFLAG_NULL);  
>                 
>         BindVar(CPList, "stresscomponentIndex", &param->stresscomponentIndex, V_INT, 1,
>                 VFLAG_NULL);        
>          
>             param->stresscomponentIndex = 1;  
>                 
>         BindVar(CPList, "avalancheLimit", &param->avalancheLimit, V_DBL, 1,
>                 VFLAG_NULL);          
>                      
>         BindVar(CPList, "stressRate", &param->stressRate, V_DBL, 1,
>                 VFLAG_NULL);     
>                 
>         BindVar(CPList, "vintegWindow", &param->vintegWindow, V_DBL, 1,
>                 VFLAG_NULL); 
>                 
>         BindVar(CPList, "nbrListFreq", &param->nbrListFreq, V_INT, 1,
>                 VFLAG_NULL);            
>         param->nbrListFreq = 10;       
>                             
850a883,884
>         
>         
863a898,944
> 
> 
> 
> void PrecipitateParamInit(Param_t *param, ParamList_t *DPList)
> {
> 
> /*
>  *      Note: Parameters need only be initialized if their
>  *      default values are non-zero.
>  */
>         DPList->paramCnt = 0;
>         DPList->varList = (VarData_t *)NULL;
> 
>         BindVar(DPList, "dataFileVersion", &param->dataFileVersion, V_INT,
>                 1, VFLAG_NULL);
>         param->dataFileVersion = NODEDATA_FILE_VERSION;
> 
>         BindVar(DPList, "numFileSegments", &param->numFileSegments, V_INT,
>                 1, VFLAG_NULL);
>         param->numFileSegments = 1;
> 
>         BindVar(DPList, "minCoordinates", param->minCoordinates, V_DBL, 3,
>                 VFLAG_NULL);
> 
>         BindVar(DPList, "maxCoordinates", param->maxCoordinates, V_DBL, 3,
>                 VFLAG_NULL);
> 
>         
>         
>         BindVar(DPList, "precipitateCount", &param->precipitateCount, V_INT, 1, VFLAG_NULL);
> 
>         BindVar(DPList, "dataDecompType", &param->dataDecompType, V_INT, 1,
>                 VFLAG_NULL);
>         param->dataDecompType = 2;
> 
>         BindVar(DPList, "dataDecompGeometry", param->dataDecompGeometry,
>                 V_INT, 3, VFLAG_NULL);
>         param->dataDecompGeometry[X] = 1;
>         param->dataDecompGeometry[Y] = 1;
>         param->dataDecompGeometry[Z] = 1;
> 
>         return;
> }
> 
> 
> 
> 
