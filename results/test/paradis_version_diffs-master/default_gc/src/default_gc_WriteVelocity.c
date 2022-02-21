28,29c28,29
<         int      i, newNodeKeyPtr;
<         real8    vx, vy, vz, burgMag;
---
>         int      i, newNodeKeyPtr,nodecount,nodeid;
>         real8    vx, vy, vz,vaver,burgMag;
85c85,87
< 
---
> 		nodecount=0;
> 		vaver=0.0;
> 		nodeid=-1;
95,97c97,104
< 
<             fprintf(fp, "%e %e %e %d  (%d,%d)\n", vx, vy, vz, node->sgnv,
<                     node->myTag.domainID, node->myTag.index);
---
>             
> 			nodecount=nodecount+1;
> 			vaver=vaver+sqrt(vx*vx+vy*vy+vz*vz);
> 			
> 			
> 			nodeid=node->myTag.domainID;
>             fprintf(fp, "%e %e %e %d  (%d,%d, %d)\n", vx, vy, vz, node->sgnv,
>                     node->myTag.domainID, node->myTag.index, node->constraint);
98a106
>         
101c109,118
< 
---
>         
>         //snprintf(fileName, sizeof(fileName), "%s/nopeudet",
>                          //DIR_PROPERTIES);
>                 //fp = fopen(fileName, "a");
>                 //fprintf(fp,"%e %e %d \n", param->timeNow, vaver/nodecount,nodeid);
>                 //fclose(fp);
>         
>         
>         
> //printf("TASSA! %e ",vaver/nodecount);
