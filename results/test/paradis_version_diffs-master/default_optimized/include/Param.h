21a22
> #define PRECIPITATEDATA_FILE_SUFFIX ".pdata"
193a195,197
>                              /*6 Oscillating stress ramp (AL) */
>                              /*7 velocity controlled quasitatic stress ramp (AL) */
>                              /*8 strainrate controlled quasistatic stress ramp*/
206a211,219
>         real8 stressAmp;	/*Stress Amplitude for oscillating stress (AL) */
>         real8 stressFreq; /* Frequency for oscillating stress (AL) */
>         real8 stressRate; /* Stress rate for quasistatic loading (AL) */
>         real8 avalancheLimit;		/*Threshold for quasistatic stress ramp (AL) */
>         real8 stopTime;  /* Time to stop the simulation (AL) */
>         real8 relaxTime;  /* Relaxation time  (AL) */
>         real8 vintegWindow;  /* Lenght of velocity integration window (AL) */
>         
>         int   stresscomponentIndex; /* index of stress component xx=1, yy=2, zz=3, */
282a296,301
>         
>  /* Neighborlist parameters*/
>          
> 		int	  nbrListFreq;		/*  Tell's you how often neigborlist is updated. Units are in cycles (AL) */	
> 
> 
288a308,317
>         
>         real8 distAverage;      /*the "distance" from the average velocity, used in relaxation calculations (AL)  */
>         real8 totdisLength;		/*Total length of dislocation network (AL)*/
>         real8 intStartStrain;				/*strain from the start of integration window (AL)*/
>         real8 integStrainrate;				/*Integrated strain rate (AL)*/
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
