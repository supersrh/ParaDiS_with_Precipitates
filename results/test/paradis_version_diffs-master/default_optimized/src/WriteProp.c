224c224
<                 if ((param->loadType == 1) || (param->loadType == 4)) {
---
>                 if ((param->loadType == 1) || (param->loadType == 4) || (param->loadType == 6) ) {
275a276,283
> 					sigijk = param->appliedStress[0]*al*al     +
>                          param->appliedStress[1]*am*am     +
>                          param->appliedStress[2]*an*an     +
>                          2.0*param->appliedStress[3]*am*an +
>                          2.0*param->appliedStress[4]*an*al +
>                          2.0*param->appliedStress[5]*al*am;
>         
>         
279,280c287,288
<                     fprintf(fp, "%e %e\n", param->timeNow,
<                             fabs(dpstnijk/param->realdt));
---
>                     fprintf(fp, "%e %e %e %e\n", param->timeNow,
>                             dpstnijk/param->realdt,sigijk,param->appliedStress[2]);
281a290,325
>                     //Added stress and removed the absoluteness of strain rate (AL)
>                      //fprintf(fp, "%e %e\n", param->timeNow,
>                             //fabs(dpstnijk/param->realdt));
>                             //New file for all strain rates and applied stresses. Format is stress1,strainrate1,stress2,straintrate2..
>                     //(AL)
>                         snprintf(fileName, sizeof(fileName), "%s/allepsdot",
>                              DIR_PROPERTIES);
>                     fp = fopen(fileName, "a");
>                     fprintf(fp, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n", param->timeNow,
>                             param->delpStrain[0]/param->realdt,param->appliedStress[0],
>                             param->delpStrain[1]/param->realdt,param->appliedStress[1],
>                             param->delpStrain[2]/param->realdt,param->appliedStress[2], 
>                             param->delpStrain[3]/param->realdt,param->appliedStress[3], 
>                             param->delpStrain[4]/param->realdt,param->appliedStress[4],
>                             param->delpStrain[5]/param->realdt,param->appliedStress[5]);
>                             
>                     fclose(fp);     
>                             
>                     
>                     
>                     
>                      snprintf(fileName, sizeof(fileName), "%s/avalanche",
>                              DIR_PROPERTIES);
>                     fp = fopen(fileName, "a");
> 					//last one is the strain rate
> 					fprintf(fp, " %e %e %e %e %e %e \n",
> 					param->timeNow,param->vAverage,param->totpStn[param->stresscomponentIndex],param->appliedStress[param->stresscomponentIndex],param->totdisLength, param->integStrainrate);
>        
>                             
>                     fclose(fp);     
>                     
>                     
>                     
>                     
>                     
>                     
