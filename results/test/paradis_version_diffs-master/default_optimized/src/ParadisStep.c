54a55,67
> typedef enum {
>         V_NODE_COUNT = 0,
>         V_AVERAGE_X, 
>         V_AVERAGE_Y, 
>         V_AVERAGE_Z, 
>         V_VAR,
>         V_DISL,
>         V_MOMENT1,
>         V_MAXSTATS 
> } VStatTypes_t;
> 
> 
> 
193a207,634
> void SetinitialIntegrationpos(Home_t *home)
> {
> 	
> 	int      i;
> 	Param_t  *param;
>     Node_t   *node,*nbrNode;
> 
> 	
> 	 for (i=0; i<home->newNodeKeyPtr;i++) {
> 
>             if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
>                 continue;
>             }
>          node->nodedx=node->x;
> 	     node->nodedy=node->y;
> 	     node->nodedz=node->z;
>             
>             
> 	}
> 	
> 	}
> 
> 
> 
> /* Modified the default version of NodeVelocityCollection(). It now calculates the segment length weighted average velocity of dislocations in a time integration window (AL) */
> 
> 
> void NodeVelocityCollection(Home_t *home)
> {
> 
>         int      i,iarm;
>         real8    velStatsLocal[V_MAXSTATS], velStatsGlobal[V_MAXSTATS];
>         real8    v2,v3, vx, vy, vz,nodeCount,burgMag;
>         real8	 rt[3],vn[3],v[3],vnbr[3];	
>         real8   xCellSize, yCellSize, zCellSize, minCellSize;
>         real8    vnbrx, vnbry, vnbrz,vproj;
>         real8    v2sq, vAveragesq, vStDev;
>         real8    vAverage, vAveragex, vAveragey, vAveragez;
>         real8	 v05,v4,rtNorm,vnNorm,maxSeg;
>         real8	 dx, dy, dz,nbrx, nbry, nbrz,disMag,seglength;
>         real8	 integdx,integdy,integdz,integdist;
>         real8	 vMoment2,vMoment3,totdisLength;
>         real8    minX,minY,minZ,maxX,maxY,maxZ;
>         Param_t  *param;
>         Node_t   *node,*nbrNode;
> 
>     param = home->param;
>     	/* Determening the minimun cell size */
>         xCellSize = param->Lx / param->nXcells;
>         yCellSize = param->Ly / param->nYcells;
>         zCellSize = param->Lz / param->nZcells;
> 
>         minCellSize = MIN(xCellSize, yCellSize);
>         minCellSize = MIN(minCellSize, zCellSize);
>       
>    
> 	burgMag = param->burgMag;
> 	maxSeg=2.0*param->maxSeg;
> 	
> 	minX=param->minCoordinates[0];
> 	minY=param->minCoordinates[1];
> 	minZ=param->minCoordinates[2];
> 	
> 	maxX=param->maxCoordinates[0];
> 	maxY=param->maxCoordinates[1];
> 	maxZ=param->maxCoordinates[2];
> 	
> 	
> /*
>  *      Zero out some stuff before starting 
>  */
>         for (i = 0; i < V_MAXSTATS; i++) {
>             velStatsLocal[i]  = 0.0;
>             velStatsGlobal[i] = 0.0;
>         }
> 
> 
> 
>  
> 
> 
> 
> 
> 
> 
> /*
>  *      Loop over all the nodes, accumulating all the necessary data
>  */
>         for (i=0; i<home->newNodeKeyPtr;i++) {
> 
>             if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
>                 continue;
>             }
> 
> 					disMag=0.0;
> 					v2=0.0;				
> 
> 					integdx=node->x-node->nodedx;
> 					integdy=node->y-node->nodedy;
> 					integdz=node->z-node->nodedz;
> 					
> 					/*PBC*/
> 					
> 					if( (abs(integdx)>minCellSize) || (abs(integdy)>minCellSize) || (abs(integdz)>minCellSize)){
> 					ZImage(param, &integdx, &integdy, &integdz);
> 					
> 						}
> 							
> 					integdist=sqrt(integdx*integdx+integdy*integdy+integdz*integdz);
> 					
> 					
> 					
> 									
> 					vx = (integdx)/(param->timeNow-param->velstatT)*burgMag;
> 					vy = (integdy)/(param->timeNow-param->velstatT)*burgMag;
> 					vz = (integdz)/(param->timeNow-param->velstatT)*burgMag;	
> 					
> 					v[0]=vx;
>                     v[1]=vy;
>                     v[2]=vz;
> 
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
>                    
>                     
>                     /*For the case of periodic boundaries. We dont want to count fictional segments of the cube size */
>                     /*PBC Korjaa ottamalla käyttöön paradisin functio*/	
>                     
>                     if( (abs(dx)>minCellSize) || (abs(dy)>minCellSize) || (abs(dz)>minCellSize)){
> 					
> 					ZImage(param, &dx, &dy, &dz);
> 						}
>                     
>                     
> 		          
>                     
>                     /*Direction vector*/
>                     rt[0]=dx;
>                     rt[1]=dy;
>                     rt[2]=dz;
>                     
>                     rtNorm = Normal(rt);
>                     
>                     seglength=(rtNorm*burgMag);
>                      /* projecting paraller velocity components away*/
>                     vproj=DotProduct(v, rt)/(rtNorm * rtNorm);
>                     
>                     vn[0]=v[0]-vproj*rt[0];
>                     vn[1]=v[1]-vproj*rt[1];
>                     vn[2]=v[2]-vproj*rt[2];
>                    				
> 					vnNorm=Normal(vn);
> 					
> 					disMag =disMag+seglength;
> 					v2 =v2+vnNorm*seglength;                          
> }
> 
> 	     node->nodedx=node->x;
> 	     node->nodedy=node->y;
> 	     node->nodedz=node->z;
> 
> 	
> 	
> 			v05 = sqrt(v2);
> 			v3=  sqrt(v2*v2*v2);
> 			v4= 1.0/v05;
> 
> /* V_moments etc are old variables that have no real meaning anymore*/
>             velStatsLocal[V_AVERAGE_X] += vx;
>             velStatsLocal[V_AVERAGE_Y] += vy;
>             velStatsLocal[V_AVERAGE_Z] += vz;
>             velStatsLocal[V_MOMENT1]   += v2/2.0;         
>             velStatsLocal[V_DISL]	   +=disMag/2.0	; 
>             
>             
>             
>             //~ velStatsLocal[V_NODE_COUNT]++;
>             velStatsLocal[V_NODE_COUNT]+=1.0;
>             
>         }
>         
>         
>         
> 
> 
> 
>         
>       
> #ifdef PARALLEL
> 			MPI_Barrier(MPI_COMM_WORLD);
>             MPI_Allreduce(velStatsLocal, velStatsGlobal, V_MAXSTATS,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);                       
>                            
>             
> 
>             nodeCount = velStatsGlobal[V_NODE_COUNT];
> 			
>             vAveragex = velStatsGlobal[V_AVERAGE_X] / nodeCount;
>             vAveragey = velStatsGlobal[V_AVERAGE_Y] / nodeCount;
>             vAveragez = velStatsGlobal[V_AVERAGE_Z] / nodeCount;
> 
>         MPI_Barrier(MPI_COMM_WORLD);
> 
> 
>             vStDev = sqrt(velStatsGlobal[V_AVERAGE_Z]/ nodeCount-(velStatsGlobal[V_VAR] / nodeCount)*(velStatsGlobal[V_VAR] / nodeCount) );
>             param->vStDev   = vStDev;
>             param->vAverage = velStatsGlobal[V_MOMENT1];
>             param->totdisLength = velStatsGlobal[V_DISL];
>            
>             
>        
>              return;
>  #endif
>  
>  
> 			nodeCount = velStatsLocal[V_NODE_COUNT];
>             param->vAverage = velStatsLocal[V_MOMENT1] / velStatsLocal[V_DISL];
>             param->totdisLength = velStatsLocal[V_DISL];
>             param->velstatTime=param->velstatTime+param->vintegWindow;
>             param->velstatT=param->timeNow;
>  
>  
>   return;
>  
>         
> }
> 
> 
> 
> 
> 
> 
> 
> 
> 
> 
> 
> 
> 
> 
> /*Function for assigning precipitates to nodes neighborlists (AL) */
> 
> void AssignPrecipitatestoNodes(Home_t *home)
> {
> 
>         int      i,ci,iarm,ip,NiP,nop,j,precipitatecount,cellnbrcount;
>         int		 *precipitateindex,*cellindex;
>         real8   xCellSize, yCellSize, zCellSize, minCellSize;
>         real8    R,xp,yp,zp,x0,y0,z0,x1,y1,z1,x2,y2,z2,matka;
>         real8    dx,dy,dz;
>         real8	 r0[3],rp[3];	
>         Param_t  *param;
>         Node_t   *node,*nbrNode;
>         Precipitate_t *precipitate,*precipitatetemp,*precipitatetool;
> 		Cell_t	 *homecell,*nbrcell,*tempcell;
> 		Tag_t		 *tags;
> 
>     param = home->param;
> 	precipitatecount=param->precipitateCount;
> 	
> 	/* Determening the minimun cell size */
>         xCellSize = param->Lx / param->nXcells;
>         yCellSize = param->Ly / param->nYcells;
>         zCellSize = param->Lz / param->nZcells;
> 
>         minCellSize = MIN(xCellSize, yCellSize);
>         minCellSize = MIN(minCellSize, zCellSize);
> 	
> 	
> 	
> 	
> 	
> 	
> 	precipitateindex=(int*) calloc (precipitatecount,sizeof(int));
> 	tags = (Tag_t *)calloc(precipitatecount,sizeof(Tag_t) );
> /*
>  *      Loop over all the nodes, accumulating all the necessary data
>  */
> 	for (i=0; i<home->newNodeKeyPtr;i++) {
> 
> 		if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
> 			continue;
>         }
>         
>         
>         //home cell
> 		homecell=home->cellKeys[node->cellIdx];
> 		cellnbrcount=homecell->nbrCount;
> 		//cellnbrcount=home->nativeCellCount;
> 		
> 		cellindex=(int*) calloc (cellnbrcount+1,sizeof(int));
> 		
> 		cellindex[0]= node->cellIdx;
> 		//List of the neighbors of the current cell
> 		
> 		for(ci=1;ci<cellnbrcount+1;ci++){
> 			cellindex[ci]=homecell->nbrList[ci-1];
> 			
> 			//cellindex[ci]=home->cellList[ci-1];
> 		}
> 		nop=0;
> 		for(ci=0;ci<cellnbrcount+1;ci++){
> 			
> 			/*Sifting for periodicity */
> 			tempcell=home->cellKeys[cellindex[ci]];
> 			if(tempcell->baseIdx>0){
> 				tempcell=home->cellKeys[tempcell->baseIdx];
> 				
> 			}	
> 			//Lisaa periodisuus naapurilistaukseen ja voiman laskentaan!
> 			//printf(" cell %d %d \n",cellindex[ci],tempcell->baseIdx);
> 			
> 			precipitate=tempcell->precipitateQ;
> 
> 			//printf("Lippu %d\n",cellindex[ci]);
> 			for ( ; precipitate != (Precipitate_t *)NULL; precipitate = precipitate->nextInCell) {
> 			//printf("node  ( %d %d )ghost %d  preciptate %d home %d  \n ",node->myTag.domainID,node->myTag.index,node->native,precipitate->native,home->myDomain);
> 			//printf("node  ( %d %d )precipitate (%d %d)ghost %d domain %d \n ",node->myTag.domainID,node->myTag.index,precipitate->myTag.domainID,precipitate->myTag.index,precipitate->native,home->myDomain);
> 			
> 				NiP=0;
> 			//radius of the neighbor search is the half of the minimun cell size or ten time the radius of the precipitate
> 				
> 				R=5.0*precipitate->r;
> 				if((0.5*minCellSize)<R){
> 					R=0.5*minCellSize;
> 					}
> 				
> 				
> 				
> 				
> 				
> 				
> 				r0[0]=node->x;
> 				r0[1]=node->y;
> 				r0[2]=node->z;
> 				
> 				rp[0]=precipitate->x;
> 				rp[1]=precipitate->y;
> 				rp[2]=precipitate->z;
> 				
> 				
> 				
> 				
> 				
> 				dx=rp[0]-r0[0];
> 				dy=rp[1]-r0[1];
> 				dz=rp[2]-r0[2];
> 				matka=sqrt(dx*dx+dy*dy+dz*dz);
> 				
> 				if((abs(dx)>minCellSize) || (abs(dy)>minCellSize) || (abs(dz)>minCellSize)){
> 				PBCPOSITION(param,	r0[0],r0[1],r0[2],&rp[0],&rp[1],&rp[2]);
> 				}
> 				
> 				
> 				
> 				
> 				
> 				NiP=NodeinPrecipitateforNbrList(R, rp[0],rp[1],rp[2],r0[0],r0[1],r0[2]);
> 				
> 				if(NiP>0){
> 					//printf("node,nbrnode,precipitate  ( %d %d ) (%d %d) (%d %d) %e %e %e %e %e %e\n ",node->myTag.domainID,node->myTag.index,nbrNode->myTag.domainID,nbrNode->myTag.index,precipitate->myTag.domainID,precipitate->myTag.index,x1,y1,z1,xp,yp,zp);
> 					precipitateindex[nop]=precipitate->myTag.index;
> 					tags[nop]=precipitate->myTag;
> 					//printf("node,precipitate,(%d,%d) (%d,%d)\n",node->myTag.domainID,node->myTag.index,precipitate->myTag.domainID,precipitate->myTag.index);
> 					
> 					nop=nop+1;
> 					
> 				}			
> 		}
> 		
> 			
> 	}		
> 					
> 			
> 			
> 			AllocNodePrecipitates(node, nop);
> 			node->numPNbrs=nop;
> 			
> 			for (j =0; j < nop; j++) {
> 				node->PnbrTag[j]=tags[j];
> 				
> 			}
> 				
> 		
> 		
> 		
> 		
> 		
> 	free(cellindex);	
> 	}
>     
> free(tags);
> free(precipitateindex);
>         
>        
>       
>         return;
> }
> 
> 
> 
> 
> 
> 
> 
> 
> 
197c638
<         int        i;
---
>         int        i,prcount;
199a641
>         real8		deltastrain;
201a644
>         Precipitate_t *precipitate;
205d647
<         param = home->param;
206a649,653
>         param = home->param;
>        
>     
> 	// printf("Task %d beginning of paradistep flag 1, cycle %d\n", home->myDomain,home->cycle);
>    // MPI_Barrier(MPI_COMM_WORLD);   	  
218a666
>  	
221c669
< 
---
> 			
227,228c675,687
<             Migrate(home);
<             RecycleGhostNodes(home);
---
>    
> #ifdef PARALLEL
> 	   MigratePrecipitate(home); 
> 	   Migrate(home);		
>        
> #endif       
> 
> #ifdef PARALLEL
>         MPI_Barrier(MPI_COMM_WORLD);
> #endif  
> 			RecycleGhostPrecipitates(home);
>             RecycleGhostNodes(home);          
>             SortNativePrecipitates(home);
229a689,692
>             CommSendPrecipitateGhosts(home);
> #ifdef PARALLEL
>         MPI_Barrier(MPI_COMM_WORLD);
> #endif
231c694
< 
---
> 			
235a699
>        
245a710,715
>  
>  
>  
>  
>  
>  
246a717,742
>    
>         
>         if (firstTime) {
> 			
> 			AssignPrecipitatestoNodes(home);
> 			 
> 			//SortNativeNodes(home);
> 						 
> 			//CommSendGhosts(home);
> 			
> 			}else{
>         
>         
> 			if(home->cycle % (param->nbrListFreq) == 0){
> 			//printf("Assigning precipitates\n");	
> 			AssignPrecipitatestoNodes(home);
> 			//SortNativeNodes(home);
> 			//CommSendGhosts(home);
> 			}
> 		}
>         
>         
>           
>       
>    
>          
259a756
>             SetinitialIntegrationpos(home);
261a759
> 			  
262a761
>         
263a763
>            
265c765,783
<         }
---
>             
>             
>          /* Calculate the averagevelocities of dislocations (AL) */
>         if(param->timeNow>=param->velstatTime){
> 			NodeVelocityCollection(home);
>     
> 	
> 
>     }
>             
>      }
>         
>           
> 		
>         
>    
>         
>         
>      
284a803,806
> 
> 
> //printf("Task %d after Timestepintegrator, cycle %d\n", home->myDomain,home->cycle);
> 
305c827,837
< 
---
>         //We use still the old variabels in param, these should be changed as soon as possible
>     if (!firstTime) {     
>      if(param->timeNow>=param->velstatTime){
> 	    /*strain rate during the integration window (AL)*/
>        param->integStrainrate=(param->totpStn[param->stresscomponentIndex]-param->intStartStrain)/(param->timeNow-param->velstatT);  
> 	  	 
>       param->velstatTime=param->velstatTime+param->vintegWindow;
>       param->velstatT=param->timeNow;
>       param->intStartStrain=param->totpStn[param->stresscomponentIndex];  
>     }
> }
403c935
< 
---
> 		
420c952,953
< 
---
>        
> 	  
425c958,971
<         Migrate(home);
---
>  
> #ifdef PARALLEL 
> 		MigratePrecipitate(home);		
> 		Migrate(home);   		          
> #endif        
> 
> 
>     // printf("aftermigrate domain %d  cycle %d  \n",home->myDomain,home->cycle);
> #ifdef PARALLEL
>         MPI_Barrier(MPI_COMM_WORLD);
> #endif
> 
> 
> 
429a976
> 		RecycleGhostPrecipitates(home);
430a978
>         
431a980
>  
434a984
> 		SortNativePrecipitates(home);
435a986,988
>         
>  
>  
436a990
>  
439a994,997
> 		CommSendPrecipitateGhosts(home);
> #ifdef PARALLEL
>         MPI_Barrier(MPI_COMM_WORLD);
> #endif
441c999,1007
< 
---
>       
>         
> 		
> 		
> 		 
> 		
> 		
> 		
> 		
451a1018
> 
460c1027
< 
---
>   
