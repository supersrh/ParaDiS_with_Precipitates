21a22
> #define PRECIPITATEDATA_FILE_SUFFIX ".pdata"
206a208,216
>         real8 stressAmp;	/*Stress Amplitude for oscillating stress (AL) */
>         real8 stressFreq; /* Frequency for oscillating stress (AL) */
>         real8 stressRate; /* Stress rate for quasistatic loading (AL) */
>         real8 avalancheLimit;		/*When the average velocity is higher than this, an avalanche is occuring (AL) */
>         real8 stopTime;  /* Time to stop the simulation (AL) */
>         real8 relaxTime;  /* Time to stop the relaxation simulation (AL) */
>         real8 vintegWindow;  /* velocity integration window (AL) */
>         
>         int   stresscomponentIndex; /* index of stress component xx=1, yy=2, zz=3, */
282a293,298
>         
>  /* Neighborlist parameters*/
>          
> 		int	  nbrListFreq;		/*  Tell's you how often neigborlist is updated. Units are in cycles */	
> 
> 
288a305,317
>         
>         real8 distAverage;      /*the "distance" from the average velocity, used in relaxation calculations (AL)  */
>         real8 vMoment2;			/* Moments of velocity (AL)*/
>         real8 vMoment3;
>         real8 vMoment05;
>         real8 distMoment2; 
>         real8 distMoment3;
>         real8 distMoment05;
> 		real8 velstatTime;
>         real8 velstatT;
>         
>         
>         
436a466,468
> 
> 		 char precipitate_data_file[MAX_STRING_LEN];  /*name of precipitate file*/
> 
444a477,479
>                                      
>                                      /*total number of precipitates (AL)*/
>         int precipitateCount;                             
531a567
> extern "C" void PrecipitateParamInit(Param_t *param, ParamList_t *PDPList);
534a571
> void PrecipitateParamInit(Param_t *param, ParamList_t *PDPList); 
