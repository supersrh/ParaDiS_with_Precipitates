/**************************************************************************
 *
 *      Module:       NodeVelocity.c
 *      Description:  Contains functions to control setting nodal
 *                    velocities, generating velocity statistics,
 *                    and setting the timestep.
 *
 ***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Home.h"
#include "Mobility.h"
#include "Util.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

/*
 *     Define the types of velocity statistics to be accumulated.
 *     NOTE: V_MAXSTATS should always be the final item in the
 *     enumerated list below.
 */
typedef enum {
        V_NODE_COUNT = 0,
        V_AVERAGE_X, 
        V_AVERAGE_Y, 
        V_AVERAGE_Z, 
        V_VAR,
        V_DISL,
        V_MOMENT05,
        V_MOMENT1,
        V_MOMENT2,
        V_MOMENT3, 
        V_MAXSTATS 
} VStatTypes_t;

/*-------------------------------------------------------------------------
 *
 *      Function:    GetVelocityStatistics
 *      Description: If gathering of velocity statistics is enabled,
 *                   gather the statistics (defined by the VStatTypes_t
 *                   above)
 *
 *------------------------------------------------------------------------*/
void GetVelocityStatistics(Home_t *home)
{
#ifdef VEL_STATISTICS
        int      i,iarm;
        real8    velStatsLocal[V_MAXSTATS], velStatsGlobal[V_MAXSTATS];
        real8    v2,v3, vx, vy, vz,nodeCount,burgMag;
       real8	 rt[3],vn[3],v[3],vnbr[3];	
       real8    vnbrx, vnbry, vnbrz,vproj;
        real8    v2sq, vAveragesq, vStDev;
        real8    vAverage, vAveragex, vAveragey, vAveragez;
		real8	 v05,v4,rtNorm,vnNorm,maxSeg;
        real8	 dx, dy, dz,nbrx, nbry, nbrz,disMag,seglength;
        real8	 vMoment2,vMoment3,vMoment05;
        Param_t  *param;
        Node_t   *node,*nbrNode;



        param = home->param;
 		burgMag = param->burgMag;
 		maxSeg=2.0*param->maxSeg;
/*
 *      Zero out some stuff before starting 
 */
        for (i = 0; i < V_MAXSTATS; i++) {
            velStatsLocal[i]  = 0.0;
            velStatsGlobal[i] = 0.0;
        }

      for (i=0; i<home->newNodeKeyPtr;i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            > disMag=0.0;
			v2=0.0;
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
								
								rt[0]=dx;
								rt[1]=dy;
								rt[2]=dz;
								
								rtNorm = Normal(rt);
								
								seglength=rtNorm;
								
								/*For the case of periodic boundaries. We dont want to count fictional segments of the cube size */
								if (seglength>maxSeg){
								seglength=0.0;
								}
								
								
								vx = node->vX*burgMag;
								vy = node->vY*burgMag;
								vz = node->vZ*burgMag;
								
								v[0]=vx;
								v[1]=vy;
								v[2]=vz;
								
								vproj=DotProduct(v, rt)/(rtNorm * rtNorm);
								
								vn[0]=v[0]-vproj*rt[0];
								vn[1]=v[1]-vproj*rt[1];
								vn[2]=v[2]-vproj*rt[2];
								
								vnNorm=Normal(vn);
			
					disMag =disMag+seglength;
								
					
					
						
			
						v2 =v2+vnNorm*seglength;                          
			}

            /*printf("disMag,seglength %e %e   \n",v2,disMag);*/	
	
			v05 = sqrt(v2);
			v3=  sqrt(v2*v2*v2);
			v4= 1.0/v05;

/*
 *        If we're gathering statistics, accumulate the necessary data.
 * Lisää momentit ja painotussumma, paradissteppiin summailu
 */
            velStatsLocal[V_AVERAGE_X] += vx;
            velStatsLocal[V_AVERAGE_Y] += vy;
            velStatsLocal[V_AVERAGE_Z] += vz;
            velStatsLocal[V_MOMENT1]   += v2;
            velStatsLocal[V_MOMENT05]  += v05;
            velStatsLocal[V_MOMENT3]   += v3;
            velStatsLocal[V_MOMENT2]   += v4;
            velStatsLocal[V_DISL]	   +=disMag	; 
            
            
            
            //~ velStatsLocal[V_NODE_COUNT]++;
            velStatsLocal[V_NODE_COUNT]+=1.0;
        }
         
         
         
 
         
         
         
 		//~ MPI_Barrier(MPI_COMM_WORLD);
          //~ if (velStatsLocal[V_NODE_COUNT] > 0) {
 			//~ printf("Domains in  beginning of IF, NODECOUNT if %e %d   \n",velStatsLocal[V_NODE_COUNT],home->myDomain);
 #ifdef PARALLEL
 MPI_Barrier(MPI_COMM_WORLD);
 #endif 
            /* MPI_Allreduce(velStatsLocal, velStatsGlobal, V_MAXSTATS,
                           MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);*/
             //~ printf("Domains in after ALLreduce  %d   \n",home->myDomain);              
 
 			//~ printf("Domains in  beginning of IF, NODECOUNT if %e %d   \n",velStatsLocal[V_NODE_COUNT],home->myDomain);
 #ifdef PARALLEL
 MPI_Barrier(MPI_COMM_WORLD);
 #endif 
 
             nodeCount = velStatsLocal[V_NODE_COUNT];
 			
             vAveragex = velStatsLocal[V_AVERAGE_X] / nodeCount;
             vAveragey = velStatsLocal[V_AVERAGE_Y] / nodeCount;
             vAveragez = velStatsLocal[V_AVERAGE_Z] / nodeCount;
 #ifdef PARALLEL
         MPI_Barrier(MPI_COMM_WORLD);
 #endif  
             
 			
             vStDev = sqrt(velStatsLocal[V_AVERAGE_Z]/ nodeCount-(velStatsLocal[V_VAR] / nodeCount)*(velStatsLocal[V_VAR] / nodeCount) );

            param->vStDev   = vStDev;
            param->vAverage = velStatsLocal[V_MOMENT1] / velStatsLocal[V_DISL];
             //printf("velocity %e   \n",velStatsLocal[V_VAR] / nodeCount);
        
 //~ MPI_Barrier(MPI_COMM_WORLD);
#endif  /* ifdef VEL_STATISTICS */

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    CalcNodeVelocities
 *      Description: Driver function that will invoke the appropriate
 *                   mobility function to update the velocity of every
 *                   native node, then apply the velocity cutoff if
 *                   applicable.
 *
 *      Arguments:
 *          zeroOnErr  Flag indicating if nodal velocity should be
 *                     zeroed for any node for which the mobility
 *                     function was unable to calculate a velocity.
 *          doAll      Flag indicating if ALL nodes are to have
 *                     veolcity recalculated, or just those nodes
 *                     that have had their forces updated.
 *
 *      Returns:  0 if velocity was successfully calculated for all
 *                  nodes
 *                1 if the mobility functions were unable to converge
 *                  on a velocity for one or more nodes.
 *
 *------------------------------------------------------------------------*/
int CalcNodeVelocities(Home_t *home, int zeroOnErr, int doAll)
{
        int     domainMobError;
        Node_t  *node;
        Param_t *param;

        TimerStart(home, CALC_VELOCITY);

        param = home->param;
        domainMobError = 0;


//#pragma omp parallel reduction(+ : domainMobError) shared(doAll, zeroOnErr)
        {
            int    i, threadMobError = 0, nodeMobError;
            int    threadID, threadIterStart, threadIterEnd;
            Node_t *threadNode;

            GetThreadIterationIndices(home->newNodeKeyPtr, &threadID,
                                      &threadIterStart, &threadIterEnd);


            for (i = threadIterStart; i < threadIterEnd; i++) {
/*
 *              If we encountered a mobility error on a previous node,
 *              we'll probably be cutting the timestep, so don't
 *              bother restting the velocity of any subsequent nodes.
 *          
 *              Note: continue the loop rather than breaking out, though,
 *              because a 'break' from the loop would prevent the compiler
 *              from threading the loop.
 */
                if (threadMobError != 0) {
                    continue;
                }

                if ((threadNode = home->nodeKeys[i]) == (Node_t *)NULL) {
                    continue;
                }

/*
 *              If we do not need to recalculate velocity on ALL nodes,
 *              skip this node unless the forces have just been updated.
 */
                if ((doAll == 0) &&
                    ((threadNode->flags & NODE_RESET_FORCES) == 0)) {
                    continue;
                }

/*
 *              We set a pointer to the appropriate mobility function
 *              during initialization, so just invoke the function now.
 *              Need two error flags.  One to indicate if the current node
 *              had a mobility error, the other (which is returned to the
 *              caller) indicates if there were *any* nodes with mobility
 *              errors in the domain.
 */
//		printf("In Novelocity, Mobility is called  \n");
		
                nodeMobError = param->mobilityFunc(home, threadNode);
//		printf("Novelocity is back \n");

                threadMobError |= nodeMobError;

/*
 *              If we had problems calculating the mobility for
 *              this node, do any special handling.
 */
                if (nodeMobError) {
#ifdef DEBUG_TIMESTEP
                    printf("Mobility error on node (%d,%d)\n",
                           threadNode->myTag.domainID, threadNode->myTag.index);
                    PrintNode(threadNode);
#endif
                    if (zeroOnErr) {
                        threadNode->vX = 0.0;
                        threadNode->vY = 0.0;
                        threadNode->vZ = 0.0;
                    }
                }

/*
 *              We also used the 'reset forces' flag to determine if we needed
 *              to recalculate velocity, but now it can be reset.
 */
                threadNode->flags &= ~NODE_RESET_FORCES;
            }

            domainMobError = domainMobError + threadMobError;

        }  /* end omp parallel section */


/*
 *      We need to zero out the 'reset forces' flag for all 
 *      ghost nodes; local nodes taken care of above.
 *
 *      Note: No threading here. Since the ghost node queue is currently
 *      a linked list, we don't have a good way of threading a loop
 *      over its nodes.
 */
        node = home->ghostNodeQ;

        while (node) {
            node->flags &= ~NODE_RESET_FORCES;
            node = node->next;
        }

		/*This was commented out because we calculate  averagevelocities in ParadisStep (AL)
         * GetVelocityStatistics(home);*/ 

        TimerStop(home, CALC_VELOCITY);
#if PARALLEL
#ifdef SYNC_TIMERS
        TimerStart(home, CALC_VELOCITY_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, CALC_VELOCITY_BARRIER);
#endif
#endif

        return(domainMobError);
}
