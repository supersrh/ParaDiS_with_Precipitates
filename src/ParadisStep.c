/*-------------------------------------------------------------------------
 *
 *      Function:     ParadisStep
 *      Description:  This function controls everything needed for a
 *                    single step of a ParaDiS simulation including
 *                    force calculations, ghost cell communications,
 *                    node migration, dynamic load balance, output
 *                    generation, etc.
 *
 *-----------------------------------------------------------------------*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "Home.h"
#include "Util.h"
#include "DisplayC.h"
#include "Comm.h"
#include "Mobility.h"
#include "Decomp.h"
#include "QueueOps.h"

#if defined _FEM | defined _FEMIMGSTRESS
#include "FEM.h"
#endif

#ifdef PARALLEL
#include "mpi.h"
#endif

/*
 *      By default, there are no runtime checks to see if all of
 *      the dislocations have annihilated themselves.  To enable
 *      a check with an abort if it happens, simply define the
 *      DEBUG_CHECK_FOR_ZERO_SEG value below to 1 rather than zero.
 */
#define DEBUG_CHECK_FOR_ZERO_SEG 0


/*
 *      For debugging only.  If DEBUG_STEP is not defined, all
 *      calls to Synchronize() will be replaced with an empty
 *      block of code, but if it is defined, the calls will
 *      be replaced with a call to syncronize the code and log
 *      a message.
 */
#ifdef DEBUG_STEP
#define Synchronize(a,b) _Synchronize((a),(b))
#else
#define Synchronize(a,b) {}
#endif
 typedef enum {
         V_NODE_COUNT = 0,
         V_AVERAGE_X, 
         V_AVERAGE_Y, 
         V_AVERAGE_Z, 
         V_VAR,
         V_DISL,
         V_MOMENT1,
         V_MAXSTATS 
 } VStatTypes_t;
/*
 *      Explicitly synchronize parallel execution via an MPI
 *      barrier and print a message when all tasks have reached
 *      the barrier.  For debug only.
 */
void _Synchronize(Home_t *home, char *msg)
{

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        if (home->myDomain == 0) {
            printf(" *** %s: All tasks synchronized\n", msg);
            fflush(NULL);
        }
}


/*-------------------------------------------------------------------------
 *
 *      Function:    ApplyDeltaStress
 *      Description: Increment the force/vel for all native nodes based
 *                   on the change in applied stress (deltaStress)
 *                   during a timestep.  Since the applied stress plays
 *                   no impact on the segment/segment interactions this
 *                   is much less expensive than recomputing all the n^2
 *                   segment interactions.
 *
 *------------------------------------------------------------------------*/
static void ApplyDeltaStress(Home_t *home, real8 deltaStress[3][3])
{
        int     i, j, nbrArm, nbrIsLocal;
        real8   x1, y1, z1, x2, y2, z2;
        real8   bx1, by1, bz1, dx, dy, dz;
        real8   f1[3], f2[3];
        Node_t  *node, *nbr;
        Param_t *param;
        
        param = home->param;


/*
 *      Loop over all native nodes
 */        
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            node = home->nodeKeys[i];
            if (!node) continue;
        
            x1=node->x;
            y1=node->y;
            z1=node->z;
        
/*
 *          For each node, recalculate forces for all the node's arms
 *          that are either owned by this node, or terminate non-locally.
 */
            for (j = 0; j < node->numNbrs; j++) {

                nbr = GetNeighborNode(home, node, j);

                if (nbr == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                nbrIsLocal = (nbr->myTag.domainID == home->myDomain);

                if (nbrIsLocal) {
                    if (OrderNodes(node, nbr) >= 0) continue;
                    nbrArm = GetArmID(home, nbr, node);
                }
        
                bx1 = node->burgX[j];
                by1 = node->burgY[j];
                bz1 = node->burgZ[j];
        
                dx=nbr->x-x1;
                dy=nbr->y-y1;
                dz=nbr->z-z1;
        
                ZImage(param, &dx, &dy, &dz) ;
        
                x2=x1+dx;
                y2=y1+dy;
                z2=z1+dz;
        
                ExtPKForce(deltaStress, bx1, by1, bz1, x1, y1, z1,
                        x2, y2, z2, f1, f2);
        
                AddtoNodeForce(node,f1);
                AddtoArmForce(node, j, f1);
        
                if (nbrIsLocal) {
                    AddtoNodeForce(nbr, f2);
                    AddtoArmForce(nbr, nbrArm, f2);
                }
            }

            (void)EvaluateMobility(home, node);
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    ReevaluateForces
 *      Description: Look for any local nodes whose force/velocity data
 *                   have been flagged as obsolete.  For all such nodes
 *                   recompute the forces and velocities.
 *
 *------------------------------------------------------------------------*/
void ReevaluateForces(Home_t *home)
{
        int     i;
        Node_t  *node;
        Param_t *param;

        param = home->param;

        for (i = 0; i < home->newNodeKeyPtr; i++) {
            node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) continue;
            if (node->flags & NODE_RESET_FORCES) {

                SetOneNodeForce(home, node);
                EvaluateMobility(home, node);

                node->flags &= (~NODE_RESET_FORCES);
            }
        }

        return;
}
 void SetinitialIntegrationpos(Home_t *home)
 {
 	
 	int      i;
 	Param_t  *param;
     Node_t   *node,*nbrNode;
 
 	
 	 for (i=0; i<home->newNodeKeyPtr;i++) {
 
             if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                 continue;
             }
          node->nodedx=node->x;
 	     node->nodedy=node->y;
 	     node->nodedz=node->z;
             
             
 	}
 	
 	}
 
 
 
 /* Modified the default version of NodeVelocityCollection(). It now calculates the segment length weighted average velocity of dislocations in a time integration window (AL) */
 
 
 void NodeVelocityCollection(Home_t *home)
 {
 
         int      i,iarm;
         real8    velStatsLocal[V_MAXSTATS], velStatsGlobal[V_MAXSTATS];
         real8    v2,v3, vx, vy, vz,nodeCount,burgMag;
         real8	 rt[3],vn[3],v[3],vnbr[3];	
         real8   xCellSize, yCellSize, zCellSize, minCellSize;
         real8    vnbrx, vnbry, vnbrz,vproj;
         real8    v2sq, vAveragesq, vStDev;
         real8    vAverage, vAveragex, vAveragey, vAveragez;
         real8	 v05,v4,rtNorm,vnNorm,maxSeg;
         real8	 dx, dy, dz,nbrx, nbry, nbrz,disMag,seglength;
         real8	 integdx,integdy,integdz,integdist;
         real8	 vMoment2,vMoment3,totdisLength;
         real8    minX,minY,minZ,maxX,maxY,maxZ;
         Param_t  *param;
         Node_t   *node,*nbrNode;
 
     param = home->param;
     	/* Determening the minimun cell size */
         xCellSize = param->Lx / param->nXcells;
         yCellSize = param->Ly / param->nYcells;
         zCellSize = param->Lz / param->nZcells;
 
         minCellSize = MIN(xCellSize, yCellSize);
         minCellSize = MIN(minCellSize, zCellSize);
       
    
 	burgMag = param->burgMag;
 	maxSeg=2.0*param->maxSeg;
 	
 	minX=param->minCoordinates[0];
 	minY=param->minCoordinates[1];
 	minZ=param->minCoordinates[2];
 	
 	maxX=param->maxCoordinates[0];
 	maxY=param->maxCoordinates[1];
 	maxZ=param->maxCoordinates[2];
 	
 	
 /*
  *      Zero out some stuff before starting 
  */
         for (i = 0; i < V_MAXSTATS; i++) {
             velStatsLocal[i]  = 0.0;
             velStatsGlobal[i] = 0.0;
         }
 
 
 
  
 
 
 
 
 
 
 /*
  *      Loop over all the nodes, accumulating all the necessary data
  */
         for (i=0; i<home->newNodeKeyPtr;i++) {
 
             if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                 continue;
             }
 
 					disMag=0.0;
 					v2=0.0;				
 
 					integdx=node->x-node->nodedx;
 					integdy=node->y-node->nodedy;
 					integdz=node->z-node->nodedz;
 					
 					/*PBC*/
 					
 					if( (abs(integdx)>minCellSize) || (abs(integdy)>minCellSize) || (abs(integdz)>minCellSize)){
 					ZImage(param, &integdx, &integdy, &integdz);
 					
 						}
 							
 					integdist=sqrt(integdx*integdx+integdy*integdy+integdz*integdz);
 					
 					
 					
 									
 					vx = (integdx)/(param->timeNow-param->velstatT)*burgMag;
 					vy = (integdy)/(param->timeNow-param->velstatT)*burgMag;
 					vz = (integdz)/(param->timeNow-param->velstatT)*burgMag;	
 					
 					v[0]=vx;
                     v[1]=vy;
                     v[2]=vz;
 
  for (iarm = 0; iarm < node->numNbrs; iarm++) {
 
                 nbrNode = GetNeighborNode(home, node, iarm);
 
                 if (nbrNode == (Node_t *)NULL) {
                     printf("WARNING: Neighbor not found at %s line %d\n",
                            __FILE__, __LINE__);
                     continue;
                 }
 
       /*          if (OrderNodes(node, nbrNode) > 0) {
 
 /*
  *                  Get the neighbor's coordinates and the line
  *                  direction vector
  */
                     nbrx = nbrNode->x;
                     nbry = nbrNode->y;
                     nbrz = nbrNode->z;
                     
                     dx = nbrx - node->x;
                     dy = nbry - node->y;
                     dz = nbrz - node->z;
                     
                    
                     
                     /*For the case of periodic boundaries. We dont want to count fictional segments of the cube size */
                     /*PBC Korjaa ottamalla käyttöön paradisin functio*/	
                     
                     if( (abs(dx)>minCellSize) || (abs(dy)>minCellSize) || (abs(dz)>minCellSize)){
 					
 					ZImage(param, &dx, &dy, &dz);
 						}
                     
                     
 		          
                     
                     /*Direction vector*/
                     rt[0]=dx;
                     rt[1]=dy;
                     rt[2]=dz;
                     
                     rtNorm = Normal(rt);
                     
                     seglength=(rtNorm*burgMag);
                      /* projecting paraller velocity components away*/
                     vproj=DotProduct(v, rt)/(rtNorm * rtNorm);
                     
                     vn[0]=v[0]-vproj*rt[0];
                     vn[1]=v[1]-vproj*rt[1];
                     vn[2]=v[2]-vproj*rt[2];
                    				
 					vnNorm=Normal(vn);
 					
 					disMag =disMag+seglength;
 					v2 =v2+vnNorm*seglength;                          
 }
 
 	     node->nodedx=node->x;
 	     node->nodedy=node->y;
 	     node->nodedz=node->z;
 
 	
 	
 			v05 = sqrt(v2);
 			v3=  sqrt(v2*v2*v2);
 			v4= 1.0/v05;
 
 /* V_moments etc are old variables that have no real meaning anymore*/
             velStatsLocal[V_AVERAGE_X] += vx;
             velStatsLocal[V_AVERAGE_Y] += vy;
             velStatsLocal[V_AVERAGE_Z] += vz;
             velStatsLocal[V_MOMENT1]   += v2/2.0;         
             velStatsLocal[V_DISL]	   +=disMag/2.0	; 
             
             
             
             //~ velStatsLocal[V_NODE_COUNT]++;
             velStatsLocal[V_NODE_COUNT]+=1.0;
             
         }
         
         
         
 
 
 
         
       
 #ifdef PARALLEL
 			MPI_Barrier(MPI_COMM_WORLD);
             MPI_Allreduce(velStatsLocal, velStatsGlobal, V_MAXSTATS,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);                       
                            
             
 
             nodeCount = velStatsGlobal[V_NODE_COUNT];
 			
             vAveragex = velStatsGlobal[V_AVERAGE_X] / nodeCount;
             vAveragey = velStatsGlobal[V_AVERAGE_Y] / nodeCount;
             vAveragez = velStatsGlobal[V_AVERAGE_Z] / nodeCount;
 
         MPI_Barrier(MPI_COMM_WORLD);
 
 
             vStDev = sqrt(velStatsGlobal[V_AVERAGE_Z]/ nodeCount-(velStatsGlobal[V_VAR] / nodeCount)*(velStatsGlobal[V_VAR] / nodeCount) );
             param->vStDev   = vStDev;
             param->vAverage = velStatsGlobal[V_MOMENT1];
             param->totdisLength = velStatsGlobal[V_DISL];
            
             
        
              return;
  #endif
  
  
 			nodeCount = velStatsLocal[V_NODE_COUNT];
             param->vAverage = velStatsLocal[V_MOMENT1] / velStatsLocal[V_DISL];
             param->totdisLength = velStatsLocal[V_DISL];
             param->velstatTime=param->velstatTime+param->vintegWindow;
             param->velstatT=param->timeNow;
  
  
   return;
  
         
 }
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 /*Function for assigning precipitates to nodes neighborlists (AL) */
 
 void AssignPrecipitatestoNodes(Home_t *home)
 {
 
         int      i,ci,iarm,ip,NiP,nop,j,precipitatecount,cellnbrcount;
         int		 *precipitateindex,*cellindex;
         real8   xCellSize, yCellSize, zCellSize, minCellSize;
         real8    R,xp,yp,zp,x0,y0,z0,x1,y1,z1,x2,y2,z2,matka;
         real8    dx,dy,dz;
         real8	 r0[3],rp[3];	
         Param_t  *param;
         Node_t   *node,*nbrNode;
         Precipitate_t *precipitate,*precipitatetemp,*precipitatetool;
 		Cell_t	 *homecell,*nbrcell,*tempcell;
 		Tag_t		 *tags;
 
     param = home->param;
 	precipitatecount=param->precipitateCount;
 	
 	/* Determening the minimun cell size */
         xCellSize = param->Lx / param->nXcells;
         yCellSize = param->Ly / param->nYcells;
         zCellSize = param->Lz / param->nZcells;
 
         minCellSize = MIN(xCellSize, yCellSize);
         minCellSize = MIN(minCellSize, zCellSize);
 	
 	
 	
 	
 	
 	
 	precipitateindex=(int*) calloc (precipitatecount,sizeof(int));
 	tags = (Tag_t *)calloc(precipitatecount,sizeof(Tag_t) );
 /*
  *      Loop over all the nodes, accumulating all the necessary data
  */
 	for (i=0; i<home->newNodeKeyPtr;i++) {
 
 		if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
 			continue;
         }
         
         
         //home cell
 		homecell=home->cellKeys[node->cellIdx];
 		cellnbrcount=homecell->nbrCount;
 		//cellnbrcount=home->nativeCellCount;
 		
 		cellindex=(int*) calloc (cellnbrcount+1,sizeof(int));
 		
 		cellindex[0]= node->cellIdx;
 		//List of the neighbors of the current cell
 		
 		for(ci=1;ci<cellnbrcount+1;ci++){
 			cellindex[ci]=homecell->nbrList[ci-1];
 			
 			//cellindex[ci]=home->cellList[ci-1];
 		}
 		nop=0;
 		for(ci=0;ci<cellnbrcount+1;ci++){
 			
 			/*Sifting for periodicity */
 			tempcell=home->cellKeys[cellindex[ci]];
 			if(tempcell->baseIdx>0){
 				tempcell=home->cellKeys[tempcell->baseIdx];
 				
 			}	
 			//Lisaa periodisuus naapurilistaukseen ja voiman laskentaan!
 			//printf(" cell %d %d \n",cellindex[ci],tempcell->baseIdx);
 			
 			precipitate=tempcell->precipitateQ;
 
 			//printf("Lippu %d\n",cellindex[ci]);
 			for ( ; precipitate != (Precipitate_t *)NULL; precipitate = precipitate->nextInCell) {
 			//printf("node  ( %d %d )ghost %d  preciptate %d home %d  \n ",node->myTag.domainID,node->myTag.index,node->native,precipitate->native,home->myDomain);
 			//printf("node  ( %d %d )precipitate (%d %d)ghost %d domain %d \n ",node->myTag.domainID,node->myTag.index,precipitate->myTag.domainID,precipitate->myTag.index,precipitate->native,home->myDomain);
 			
 				NiP=0;
 			//radius of the neighbor search is the half of the minimun cell size or ten time the radius of the precipitate
 				
 				R=5.0*precipitate->r;
 				if((0.5*minCellSize)<R){
 					R=0.5*minCellSize;
 					}
 				
 				
 				
 				
 				
 				
 				r0[0]=node->x;
 				r0[1]=node->y;
 				r0[2]=node->z;
 				
 				rp[0]=precipitate->x;
 				rp[1]=precipitate->y;
 				rp[2]=precipitate->z;
 				
 				
 				
 				
 				
 				dx=rp[0]-r0[0];
 				dy=rp[1]-r0[1];
 				dz=rp[2]-r0[2];
 				matka=sqrt(dx*dx+dy*dy+dz*dz);
 				
 				if((abs(dx)>minCellSize) || (abs(dy)>minCellSize) || (abs(dz)>minCellSize)){
 				PBCPOSITION(param,	r0[0],r0[1],r0[2],&rp[0],&rp[1],&rp[2]);
 				}
 				
 				
 				
 				
 				
 				NiP=NodeinPrecipitateforNbrList(R, rp[0],rp[1],rp[2],r0[0],r0[1],r0[2]);
 				
 				if(NiP>0){
 					//printf("node,nbrnode,precipitate  ( %d %d ) (%d %d) (%d %d) %e %e %e %e %e %e\n ",node->myTag.domainID,node->myTag.index,nbrNode->myTag.domainID,nbrNode->myTag.index,precipitate->myTag.domainID,precipitate->myTag.index,x1,y1,z1,xp,yp,zp);
 					precipitateindex[nop]=precipitate->myTag.index;
 					tags[nop]=precipitate->myTag;
 					//printf("node,precipitate,(%d,%d) (%d,%d)\n",node->myTag.domainID,node->myTag.index,precipitate->myTag.domainID,precipitate->myTag.index);
 					
 					nop=nop+1;
 					
 				}			
 		}
 		
 			
 	}		
 					
 			
 			
 			AllocNodePrecipitates(node, nop);
 			node->numPNbrs=nop;
 			
 			for (j =0; j < nop; j++) {
 				node->PnbrTag[j]=tags[j];
 				
 			}
 				
 		
 		
 		
 		
 		
 	free(cellindex);	
 	}
     
 free(tags);
 free(precipitateindex);
         
        
       
         return;
 }

void ParadisStep(Home_t *home)
{
        int        i,prcount;
        int        doAll = 1;
        real8      deltaStress[3][3];
		real8		deltastrain;
        Param_t    *param;
        Node_t     *node;
		Precipitate_t *precipitate;
        static int firstTime = 1;
   

        param = home->param;
        
     
 	// printf("Task %d beginning of paradistep flag 1, cycle %d\n", home->myDomain,home->cycle);
    // MPI_Barrier(MPI_COMM_WORLD);   	
/*
 *      If this step is an initial load-balance-only step just
 *      perform the minimal work needed to estimate per-process
 *      load, shift boundaries, and migrate nodes among processors.
 */
        if (param->numDLBCycles > 0) {
/*
 *          Note:  When param->numDLBCycles > 0, NodeForce() assumes
 *          the cycle is a DLB-only cycle and only counts the number
 *          of force calcs that would be done without actually calculating
 *          any forces.
 */
            NodeForce(home, FULL);
            Rebalance(home, DLB_USE_FORCECALC_COUNT);

/*
 *          Any time the boundaries are changed, we need to migrate
 *          nodes to their new owning domains and go through all the
 *          ghost node communication stuff.
 */
    
 #ifdef PARALLEL
 	   MigratePrecipitate(home); 
 	   Migrate(home);		
        
 #endif       
 
 #ifdef PARALLEL
         MPI_Barrier(MPI_COMM_WORLD);
 #endif  
 			RecycleGhostPrecipitates(home);
             RecycleGhostNodes(home);          
             SortNativePrecipitates(home);
            SortNativeNodes(home);
			CommSendPrecipitateGhosts(home);
 #ifdef PARALLEL
         MPI_Barrier(MPI_COMM_WORLD);
 #endif
            CommSendGhosts(home);

            home->cycleForceCalcCount = 0;

            return;
        }

#ifndef NO_XWINDOW
        while (WinIsPaused()) {
            sleep(1);
        }
#endif

/*
 *      Calculate the net charge tensor for each cell (includes global comm)
 */
        CellCharge(home);
		if (firstTime) {
			
			AssignPrecipitatestoNodes(home);
			 
			//SortNativeNodes(home);
						 
			//CommSendGhosts(home);
			
			}else{
        
        
			if(home->cycle % (param->nbrListFreq) == 0){
			//printf("Assigning precipitates\n");	
			AssignPrecipitatestoNodes(home);
			//SortNativeNodes(home);
			//CommSendGhosts(home);
			}
		}

/*
 *      Calculate new force and velocity data for all nodes or a selected
 *      subset and distribute the new data out to neighboring domains.
 *      The first cycle we'll recalculate all forces and velocities
 *      We do this to get an initial estimate of forces on the first cycle,
 *      After that, we only need to recompute values for nodes that were
 *      involved in topological changes the previous step.
 */
        if (firstTime) {
            NodeForce(home, FULL);
            CalcNodeVelocities(home, 1, doAll);
            CommSendVelocity(home);
			SetinitialIntegrationpos(home);
            firstTime = 0;
        } else {
            NodeForce(home, PARTIAL);
            CalcNodeVelocities(home, 0, doAll);
            CommSendVelocity(home);
            
            
         /* Calculate the averagevelocities of dislocations (AL) */
        if(param->timeNow>=param->velstatTime){
			NodeVelocityCollection(home);
    
	

    }
            
     }

/*
 *      Invoke the selected time step integration method.  The
 *      selected method will calculate the time step as well as
 *      move nodes to their correct locations and communicate 
 *      the new nodal force/velocity data to neighboring domains.
 */
        if (strcmp(param->timestepIntegrator, "forward-euler") == 0) {
            ForwardEulerIntegrator(home);
        } else if (strcmp(param->timestepIntegrator, "trapezoid") == 0) {
            TrapezoidIntegrator(home);
        } else {
/*
 *          Used to be specified as 'backard-euler', so if integration
 *          method is unrecognized, use trapezoid as default
 */
            TrapezoidIntegrator(home);
        }
//printf("Task %d after Timestepintegrator, cycle %d\n", home->myDomain,home->cycle);
/*
 *      In some simulations, it is necessary to recalculate and distribute
 *      the glide plane infromation for segments after the nodes have been
 *      repositioned.  Do so now if needed.
 */
        ResetGlidePlanes(home);

/*
 *      Increment the per-burgers vector density gain/loss with
 *      changes for this cycle.  This must be done immediately
 *      after timestep integration!
 *
 *      Note: This is currently only applicable to BCC simulations.
 */
        GetDensityDelta(home);

/*
 *      Calculate the new plastic strain.
 */
        DeltaPlasticStrain(home);
         //We use still the old variabels in param, these should be changed as soon as possible
     if (!firstTime) {     
      if(param->timeNow>=param->velstatTime){
 	    /*strain rate during the integration window (AL)*/
        param->integStrainrate=(param->totpStn[param->stresscomponentIndex]-param->intStartStrain)/(param->timeNow-param->velstatT);  
 	  	 
       param->velstatTime=param->velstatTime+param->vintegWindow;
       param->velstatT=param->timeNow;
       param->intStartStrain=param->totpStn[param->stresscomponentIndex];  
     }
 }
/*
 *      The call to GenerateOutput will update the time and cycle counters,
 *      determine if any output needs to be generated at this stage, and
 *      call the appropriate I/O functions if necessary.
 */
        GenerateOutput(home, STAGE_CYCLE);

/*
 *      Before doing topological changes, set flags indicating any
 *      nodes exempt from topological changes.  These flags are used
 *      in both splitting multi-arm nodes and collisions, so this
 *      function should be invoked before either of those items are done.
 */
        InitTopologyExemptions(home);

/*
 *      Now do all the topological changes from segment interactions
 *      (collisions, multinode splitting)...  Clear the list of local
 *      operations that will be sent to the remote domains for processsing,
 *      then split any multi-arm nodes that need splitting, cross slip
 *      nodes (as needed/allowed), handle all local collisions, then
 *      send remote nodes the list of ops needed to keep their data in sync.
 */
        ClearOpList(home);
        SortNodesForCollision(home);

        SplitMultiNodes(home);

#ifdef _FEM
        SplitSurfaceNodes(home);
#endif

/*
 *      Call a generic cross-slip dispatch function that will call
 *      (if necessary) the cross-slip function appropriate to the
 *      type of material in use.
 */
        CrossSlip(home);

/*
 *      Search for dislocation segments in close proximity to each other
 *      and if necessary handle any collision between them.
 */
        HandleCollisions(home);
 
#ifdef PARALLEL
#ifdef SYNC_TIMERS
        TimerStart(home, POST_COLLISION_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, POST_COLLISION_BARRIER);
#endif
#endif

        TimerStart(home, COL_SEND_REMESH);
        CommSendRemesh(home);
        TimerStop(home, COL_SEND_REMESH);

        TimerStart(home, COL_FIX_REMESH);
        FixRemesh(home);
        TimerStop(home, COL_FIX_REMESH);

#ifdef _FEM
	AdjustNodePosition(home, 1); 
#endif 

/*
 *      Under certain circumstances, parallel topological changes can
 *      create double links between nodes; links which can not be detected
 *      until after FixRemesh() is called... so, a quick check has to be
 *      done to clean up these potential double-links here, or they will
 *      cause problems later on.  Should only have to check nodes local
 *      to this domain.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
            (void)RemoveDoubleLinks(home, node, 0);
            node->flags &= ~NODE_CHK_DBL_LINK;
        }

/*
 *      If memory debugging is enabled, run a consistency check on all
 *      allocated memory blocks looking for corruption in any of the
 *      block headers or trailers.
 */
#ifdef DEBUG_MEM
        ParadisMemCheck();
#endif
       
/*
 *      Invoke mesh coarsen/refine
 */
        Remesh(home);

/*
 *      Define load curve and calculate change in applied stress this cycle
 */
        LoadCurve(home, deltaStress);

/*
 *      This is only needed when we do force calcs for only
 *      a subset of the nodes at the beginning of the timestep.  It will
 *      adjust the nodal forces based on the current delta stress and
 *      recalculate the nodal velocities so we have more accurate values
 *      when we enter the timestep integrator at the beginning of the next
 *      cycle.
 */
        ApplyDeltaStress(home, deltaStress);

/*
 *      If necessary, use the current load data to generate a new
 *      domain decomposition to rebalance the workload among the
 *      processors.
 */
        Rebalance(home, DLB_USE_WALLCLK_TIME);

/*
 *      Send any nodes that have moved beyond the domain's
 *      boundaries to the domain the node now belongs to.
 */
 #ifdef PARALLEL 
 		MigratePrecipitate(home);		
 		Migrate(home);   		          
 #endif        
 
 
     // printf("aftermigrate domain %d  cycle %d  \n",home->myDomain,home->cycle);
 #ifdef PARALLEL
         MPI_Barrier(MPI_COMM_WORLD);
 #endif

/*
 *      Recycle all the ghost nodes: move them back to the free Queue
 */
		RecycleGhostPrecipitates(home);
        RecycleGhostNodes(home);

/*
 *      Sort the native nodes into their proper subcells.  
 */
		SortNativePrecipitates(home);
        SortNativeNodes(home);

/*
 *      Communicate ghost cells to/from neighbor domains
 */
 		CommSendPrecipitateGhosts(home);
 #ifdef PARALLEL
         MPI_Barrier(MPI_COMM_WORLD);
 #endif
        CommSendGhosts(home);

#ifdef NAN_CHECK
/*
 *      For debug only:  Abort if any of the nodes have position or
 *      velocity values that are NaNs or infinites.  Be sure to do this
 *      before we write the restart files so we don't get bad data
 *      in the restart.
 */
        CheckForNANS(home);
#endif

/*
 *      If memory debugging is enabled, run a consistency check on all
 *      allocated memory blocks looking for corruption in any of the
 *      block headers or trailers.
 */
#ifdef DEBUG_MEM
        ParadisMemCheck();
#endif

        CheckMemUsage(home, "ParadisStep-complete");

/*
 *      Zero out the count of force calculations done this cycle
 *      so the load-balancing in the next step is based on accurate
 *      values.
 */
        home->cycleForceCalcCount = 0;

/*
 *      For Debugging and testing only...
 */
#if DEBUG_CHECK_FOR_ZERO_SEG
        CheckForEmptySimulation(home);
#endif

        return;
}
