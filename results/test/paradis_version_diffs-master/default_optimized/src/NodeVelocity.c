29,31c29,31
<         V_AVERAGE_X,
<         V_AVERAGE_Y,
<         V_AVERAGE_Z,
---
>         V_AVERAGE_X, 
>         V_AVERAGE_Y, 
>         V_AVERAGE_Z, 
33c33,38
<         V_MAXSTATS
---
>         V_DISL,
>         V_MOMENT05,
>         V_MOMENT1,
>         V_MOMENT2,
>         V_MOMENT3, 
>         V_MAXSTATS 
47c52
<         int      i, nodeCount;
---
>         int      i,iarm;
49c54,56
<         real8    v2, vx, vy, vz;
---
>         real8    v2,v3, vx, vy, vz,nodeCount,burgMag;
>         real8	 rt[3],vn[3],v[3],vnbr[3];	
>         real8    vnbrx, vnbry, vnbrz,vproj;
51a59,61
>         real8	 v05,v4,rtNorm,vnNorm,maxSeg;
>         real8	 dx, dy, dz,nbrx, nbry, nbrz,disMag,seglength;
>         real8	 vMoment2,vMoment3,vMoment05;
53,54c63
<         Node_t   *node;
< 
---
>         Node_t   *node,*nbrNode;
58c67,68
< 
---
> 		burgMag = param->burgMag;
> 		maxSeg=2.0*param->maxSeg;
60c70
<  *      Zero out some stuff before starting
---
>  *      Zero out some stuff before starting 
67,70c77,82
< /*
<  *      Loop over all the nodes, accumulating all the necessary data
<  */
<         for (i=0; i<home->newNodeKeyPtr;i++) {
---
> 
> 
> 
> 
> 
>      for (i=0; i<home->newNodeKeyPtr;i++) {
76,78c88,151
<             vx = node->vX;
<             vy = node->vY;
<             vz = node->vZ;
---
> disMag=0.0;
> v2=0.0;
>  for (iarm = 0; iarm < node->numNbrs; iarm++) {
> 
>                 nbrNode = GetNeighborNode(home, node, iarm);
> 
>                 if (nbrNode == (Node_t *)NULL) {
>                     printf("WARNING: Neighbor not found at %s line %d\n",
>                            __FILE__, __LINE__);
>                     continue;
>                 }
> 
>       /*          if (OrderNodes(node, nbrNode) > 0) {
> 
> /*
>  *                  Get the neighbor's coordinates and the line
>  *                  direction vector
>  */
>                     nbrx = nbrNode->x;
>                     nbry = nbrNode->y;
>                     nbrz = nbrNode->z;
>                     
>                     dx = nbrx - node->x;
>                     dy = nbry - node->y;
>                     dz = nbrz - node->z;
>                     
>                     rt[0]=dx;
>                     rt[1]=dy;
>                     rt[2]=dz;
>                     
>                     rtNorm = Normal(rt);
>                     
>                     seglength=rtNorm;
>                     
>                     /*For the case of periodic boundaries. We dont want to count fictional segments of the cube size */
>                    if (seglength>maxSeg){
> 					   seglength=0.0;
> 					   }
>                     
>                     
>                     vx = node->vX*burgMag;
> 					vy = node->vY*burgMag;
> 					vz = node->vZ*burgMag;
> 					
> 					v[0]=vx;
>                     v[1]=vy;
>                     v[2]=vz;
>                     
>                     vproj=DotProduct(v, rt)/(rtNorm * rtNorm);
>                     
>                     vn[0]=v[0]-vproj*rt[0];
>                     vn[1]=v[1]-vproj*rt[1];
>                     vn[2]=v[2]-vproj*rt[2];
> 					
> 					vnNorm=Normal(vn);
> 
>           disMag =disMag+seglength;
>                     
>          
>           
> 			
> 
>             v2 =v2+vnNorm*seglength;                          
> }
80c153,157
<             v2 = vx*vx+vy*vy+vz*vz;
---
>  /*printf("disMag,seglength %e %e   \n",v2,disMag);*/	
> 	
> 			v05 = sqrt(v2);
> 			v3=  sqrt(v2*v2*v2);
> 			v4= 1.0/v05;
84a162
>  * Lisää momentit ja painotussumma, paradissteppiin summailu
89,90c167,176
<             velStatsLocal[V_VAR]       += v2;
<             velStatsLocal[V_NODE_COUNT]++;
---
>             velStatsLocal[V_MOMENT1]   += v2;
>             velStatsLocal[V_MOMENT05]  += v05;
>             velStatsLocal[V_MOMENT3]   += v3;
>             velStatsLocal[V_MOMENT2]   += v4;
>             velStatsLocal[V_DISL]	   +=disMag	; 
>             
>             
>             
>             //~ velStatsLocal[V_NODE_COUNT]++;
>             velStatsLocal[V_NODE_COUNT]+=1.0;
92,109c178,210
< 
<         if (velStatsLocal[V_NODE_COUNT] > 0) {
< 
<             MPI_Allreduce(velStatsLocal, velStatsGlobal, V_MAXSTATS,
<                           MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
< 
<             nodeCount = velStatsGlobal[V_NODE_COUNT];
< 
<             vAveragex = velStatsGlobal[V_AVERAGE_X] / nodeCount;
<             vAveragey = velStatsGlobal[V_AVERAGE_Y] / nodeCount;
<             vAveragez = velStatsGlobal[V_AVERAGE_Z] / nodeCount;
<  
<             vAveragesq = vAveragex*vAveragex +
<                          vAveragey*vAveragey +
<                          vAveragez*vAveragez;
<  
<             vStDev = sqrt(velStatsGlobal[V_VAR] / nodeCount) - vAveragesq;
< 
---
>         
>         
>         
> 
>         
>         
>         
> 		//~ MPI_Barrier(MPI_COMM_WORLD);
>          //~ if (velStatsLocal[V_NODE_COUNT] > 0) {
> 			//~ printf("Domains in  beginning of IF, NODECOUNT if %e %d   \n",velStatsLocal[V_NODE_COUNT],home->myDomain);
> #ifdef PARALLEL
> MPI_Barrier(MPI_COMM_WORLD);
> #endif 
>            /* MPI_Allreduce(velStatsLocal, velStatsGlobal, V_MAXSTATS,
>                           MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);*/
>             //~ printf("Domains in after ALLreduce  %d   \n",home->myDomain);              
> 
> 			//~ printf("Domains in  beginning of IF, NODECOUNT if %e %d   \n",velStatsLocal[V_NODE_COUNT],home->myDomain);
> #ifdef PARALLEL
> MPI_Barrier(MPI_COMM_WORLD);
> #endif 
> 
>             nodeCount = velStatsLocal[V_NODE_COUNT];
> 			
>             vAveragex = velStatsLocal[V_AVERAGE_X] / nodeCount;
>             vAveragey = velStatsLocal[V_AVERAGE_Y] / nodeCount;
>             vAveragez = velStatsLocal[V_AVERAGE_Z] / nodeCount;
> #ifdef PARALLEL
>         MPI_Barrier(MPI_COMM_WORLD);
> #endif  
>             
> 			
>             vStDev = sqrt(velStatsLocal[V_AVERAGE_Z]/ nodeCount-(velStatsLocal[V_VAR] / nodeCount)*(velStatsLocal[V_VAR] / nodeCount) );
111,112c212,219
<             param->vAverage = sqrt(vAveragesq);
<         }
---
>             param->vAverage = velStatsLocal[V_MOMENT1] / velStatsLocal[V_DISL];
>              //printf("velocity %e   \n",velStatsLocal[V_VAR] / nodeCount);
>         //~ }
>         
>       //~ printf("Domains in  velocitystastistics after if %d %e  \n",home->myDomain,velStatsLocal[V_NODE_COUNT]);
>  
>         
>         
115a223
> //~ MPI_Barrier(MPI_COMM_WORLD);
246c354,361
<         GetVelocityStatistics(home);
---
> 
> 
> 
> 
> 
> 
>         /*This was commented out because we calculate  averagevelocities in ParadisStep (AL)
>          * GetVelocityStatistics(home);*/ 
