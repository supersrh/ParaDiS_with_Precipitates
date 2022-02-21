37c37,38
<  *                                   is implemented.   
---
>  *                                   is implemented. 
>  * 						02/01/2012-  Lehtinen Oscillatory stress implemented  
55c56
<  *      Function:     SpinMatrix
---
>  *      function:     SpinMatrix
91a93,94
>         real8   stressFreq,stressRate;
>         real8   stressAmp;
122a126,127
>         stressFreq		= param->stressFreq;
>         stressAmp		= param->stressAmp;
551a557,562
>             
>             
>   /* Oscillatory stressramp (AL)*/         
>             if(param->timeNow>= param->relaxTime){
>                 param->appliedStress[param->stresscomponentIndex] = stressAmp*sin((timeNow-param->relaxTime)*stressFreq);
> 			 }
552a564,576
>    
>              case 7:
>   /* Velocity controlled quasistatic stressramp (AL)*/ 
> 		stressRate=param->stressRate;
> 
> 				if(param->vAverage <= param->avalancheLimit && param->timeNow>= param->relaxTime )	{
> 					param->appliedStress[param->stresscomponentIndex] =param->appliedStress[param->stresscomponentIndex]+stressRate*(param->realdt);
> 					}
>               
>                
>             
>                 break;
>      
553a578,593
> 	  case 8:
> 
> 		/* Strain rate controlled quasistatic stressramp (AL)*/ 
> 		stressRate=param->stressRate;
> 
> 				if(param->integStrainrate <= param->avalancheLimit && param->timeNow>= param->relaxTime )	{
> 					param->appliedStress[param->stresscomponentIndex] =param->appliedStress[param->stresscomponentIndex]+stressRate*(param->realdt);
> 					}
>                 
>                
>             
>                 break;
> 	
> 
> 
> 
