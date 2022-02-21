54c54
<         int            armCount, totNodeCount, valCount;
---
>         int            armCount,nbrprecipitateCount,totNodeCount, valCount;
56c56
<         int            iNbr, outBufLen;
---
>         int            iNbr, iPNbr,outBufLen;
79c79
< 
---
> 			nbrprecipitateCount= 0;
90a91,92
>                     nbrprecipitateCount +=  node->numPNbrs ;
>                      
101a104
> 					   FLTS_PER_GHOST_PNRB	* nbrprecipitateCount+
105a109
>            
142c146,147
< 
---
> 					int PnbrCount;
> 					
144c149
< 
---
> 					PnbrCount = node->numPNbrs;	
149a155
>                     outBuf[bufIndex++] = node->numPNbrs;
181a188,210
>                     /* Precipitate neighbors (AL)*/
>                      for (iPNbr = 0; iPNbr < PnbrCount; iPNbr++) {
> 
> /*
>  *                  IMPORTANT!  DO NOT REMOVE the following NULL pointer
>  *                              check!  On the BG/P system, this loop was
>  *                              having issues (core dumps) accessing the
>  *                              nbrTag array.  It appears there may be a
>  *                              compiler bug causing the problem.  Adding
>  *                              the 'if' test with the call to Fatal()
>  *                              before accessing nbrTag seems to get us
>  *                              past the problem.  Sigh.
>  */
>                         if (node->PnbrTag == (Tag_t *)NULL) {
>                             Fatal("nbrTag NULL in CommSendGhosts");
>                         }
> 
>                         outBuf[bufIndex++] = node->PnbrTag[iPNbr].domainID;
>                         outBuf[bufIndex++] = node->PnbrTag[iPNbr].index;
>  
>                         
>                     }
>                                   
185a215,218
>                     
>                     outBuf[bufIndex++] = node->nodedx;
>                     outBuf[bufIndex++] = node->nodedy;
>                     outBuf[bufIndex++] = node->nodedz;
222c255
<         
---
>         //printf("bufindex pack %d \n",bufIndex);
233a267,429
> /*-------------------------------------------------------------------------
>  *
>  *      Function:    CommPackPrecipitateGhosts
>  *
>  *      Description: For each neighbor domain, pack into a send buffer
>  *                   all the entities in any cell that borders the neighbor
>  *
>  *------------------------------------------------------------------------*/
> static void CommPackPrecipitateGhosts(Home_t *home) 
> {
> #ifdef PARALLEL
>         int            armCount,nbrprecipitateCount,totPrecipitateCount, valCount;
>         int            idst, domainIdx, bufIndex, iCell, cellIdx;
>         int            iNbr, iPNbr,outPBufLen;
>         int            precipitatesPacked;
>         real8          *outPBuf;
>         Precipitate_t         *precipitate;
>         Cell_t         *cell;
>         RemoteDomain_t *remDom;
> 
> /*
>  *      Loop through all domains neighboring the current domain
>  *      and package up for those remote domains, the precipitate data
>  *      for all local precipitates in cells bordering those domains.
>  */
>         for (idst = 0; idst < home->remoteDomainCount; idst++) {
> 
>             domainIdx = home->remoteDomains[idst];
>             remDom = home->remoteDomainKeys[domainIdx];
> 
> /*
>  *          Get a count of the total precipitates, arms, and inclusions in
>  *          all cells being exported to this remote domain.
>  */
>             armCount = 0;
>             totPrecipitateCount = 0;
> 			nbrprecipitateCount= 0;
>             for (iCell = 0 ; iCell < remDom->numExpCells ; iCell++) {
> 
>                 cellIdx = remDom->expCells[iCell];
>                 cell = home->cellKeys[cellIdx];
> 
>                 totPrecipitateCount += cell->precipitateCount;
> 
>                 precipitate = cell->precipitateQ ;
> 
>                 
> 
>             }
> 
> /*
>  *          Calculate the buffer size needed to hold the data for
>  *          this remote domain and allocate it.
>  */
>             valCount = FLTS_PER_GHOST_PRECIPITATE * totPrecipitateCount +
>                        FLTS_PER_GHOST_CELL * remDom->numExpCells +
>                        EXTRA_GHOST_FLTS;
> 
>             outPBufLen = (valCount * sizeof(real8));
>            
>             outPBuf = (real8 *)malloc(outPBufLen);
> 
>             bufIndex = 0;
>         
> /*
>  *          Send the remote domain the number of precipitates it will receive, 
>  *          the maximum possible tag.index it will receive, and the
>  *          number of cells in the message
>  */
>             bufIndex++;    /* reserve space for precipitate count */
>             bufIndex++;    /* reserve space */
> 
>             outPBuf[bufIndex++] = home->newPrecipitateKeyPtr;
>             outPBuf[bufIndex++] = remDom->numExpCells;
>         
>             for (iCell = 0; iCell < remDom->numExpCells; iCell++) {
>         
>                 cellIdx = remDom->expCells[iCell];
>                 cell = home->cellKeys[cellIdx];
>         
> /*
>  *              Supply the cell Index and the number of precipitates to follow
>  *              for this cell
>  */
>  						
>                 outPBuf[bufIndex++] = cellIdx;
>                // printf("pack %d %d %d \n",cellIdx,bufIndex,home->myDomain);
>                 outPBuf[bufIndex++] = cell->precipitateCount;
>                 outPBuf[bufIndex++] = 0;
>                 
> /*
>  *              loop through all the precipitates on the cell's precipitate queue
>  */
>                 precipitatesPacked = 0;
>                 precipitate = cell->precipitateQ;
> 
>                 while (precipitate) {
>                   
> 					
>                    
> /*
>  *                  precipitate's domainID is known implicitly, don't need to send it
>  */
>                     outPBuf[bufIndex++] = precipitate->myTag.index;
>                     outPBuf[bufIndex++] = precipitate->numNbrs;                               
> 
>                     outPBuf[bufIndex++] = precipitate->x;
>                     outPBuf[bufIndex++] = precipitate->y;
>                     outPBuf[bufIndex++] = precipitate->z;
>                     
>                     
>         
>                     outPBuf[bufIndex++] = precipitate->fX;
>                     outPBuf[bufIndex++] = precipitate->fY;
>                     outPBuf[bufIndex++] = precipitate->fZ;
>         
>                     outPBuf[bufIndex++] = precipitate->vX;
>                     outPBuf[bufIndex++] = precipitate->vY;
>                     outPBuf[bufIndex++] = precipitate->vZ;
>         
> 					outPBuf[bufIndex++] = precipitate->forcep;
>                     outPBuf[bufIndex++] = precipitate->r;
>                            
>         
>                     outPBuf[bufIndex++] = precipitate->oldvX;
>                     outPBuf[bufIndex++] = precipitate->oldvY;
>                     outPBuf[bufIndex++] = precipitate->oldvZ;
>         
>                     outPBuf[bufIndex++] = (real8)precipitate->constraint;
>                     outPBuf[bufIndex++] = (real8)precipitate->flags;
> 
> #ifdef _FEM
>                     outPBuf[bufIndex++] = (real8)precipitate->fem_Surface[0];
>                     outPBuf[bufIndex++] = (real8)precipitate->fem_Surface[1];
>                     outPBuf[bufIndex++] = precipitate->fem_Surface_Norm[0];
>                     outPBuf[bufIndex++] = precipitate->fem_Surface_Norm[1];
>                     outPBuf[bufIndex++] = precipitate->fem_Surface_Norm[2];
> #endif
>                     precipitate = precipitate->nextInCell;
>                     precipitatesPacked++;
>         
>                 }  /* end while (precipitate) */
> 
>                 if (precipitatesPacked != cell->precipitateCount) {
>                     Fatal("CommPackGhosts: dom %d, remDom %d, cell %d, "
>                           "cell->precipitateCount %d doesn't match cell->precipitateQ, "
>                           "precipitatesPacked value %d", home->myDomain, domainIdx,
>                           cellIdx, cell->precipitateCount, precipitatesPacked);
>                 }
> 
>             }   /* end for (iCell = 0 ...) */
> 			//printf("pack precipitates bufindex  %d %d %d \n",bufIndex,valCount,home->myDomain);	
>             outPBuf[0] = totPrecipitateCount;
>             outPBuf[1] = 0;  /* reserved */
>         
>             remDom->outPBuf = (char *)outPBuf;
>             remDom->outPBufLen = outPBufLen;
>         
>         } /* end for (idst = 0; ...) */
>         
> #endif
>         return;
> }
254c450
<         int            iCell, cellIdx, cellNodeCount, iNbr, numNbrs;
---
>         int            iCell, cellIdx, cellNodeCount, iNbr,iPNbr, numNbrs,numPNbrs;
263a460
>  
267c464,467
<         
---
>             
>          //if(home->myDomain == 56  ){
>            //             printf("unpackghosts myDomain %d remotedomain %d \n",home->myDomain,domainIdx);
>              //           }
276c476
<         
---
>          
305a506
>                
331c532,534
< 
---
>                     numPNbrs               = inBuf[bufIndex++];
> 		    
> 	//printf("unpackghosts myDomain %d tag ( %d %d )  \n",home->myDomain,node->myTag.domainID,node->myTag.index);	
333c536,537
< 
---
> 					AllocNodePrecipitates(node,numPNbrs );
> 					
335c539,542
< 
---
> 						
> 				
> 						
> 						
337a545,548
>                         
>                         
>                         
>                         
350a562,576
>                     
>                     /* Precipitate neighbors (AL)*/
>                      for (iPNbr = 0; iPNbr < numPNbrs; iPNbr++) {
>                         
>                         node->PnbrTag[iPNbr].domainID=inBuf[bufIndex++];
>                         node->PnbrTag[iPNbr].index=inBuf[bufIndex++];
>  
>                         
>                     }
>                     
>                     
>                     
>                     
>                     
>                     
354a581,584
>                     
>                     node->nodedx = inBuf[bufIndex++];
>                     node->nodedy = inBuf[bufIndex++];
>                     node->nodedz = inBuf[bufIndex++];
395a626,737
> 		//printf("bufindex unpack %d \n",bufIndex);
>         }  /* end for (isrc = 0; ...) */
>         
> #endif
>         return;
> }
> 
> 
> /*-------------------------------------------------------------------------
>  *
>  *  Function    : CommUnpackPrecipitateGhosts
>  *  Description : For each remote domain, unpack the comm packet which was
>  *                just received into precipitates, and queue the precipitates on the
>  *                ghost precipitate queue
>  *
>  *  Updates:   09/06/01 - add invoice stuff, to support velocity comm - t.p.
>  *             09/14/01 - changed name from CommUnpackPrecipitates, and moved into
>  *                        CommSendGhosts.c file. Converted to arbitrary
>  *                        number of arms for each precipitate - t.p.
>  *
>  *-------------------------------------------------------------------------*/
> static void CommUnpackPrecipitateGhosts(Home_t *home) 
> {
> #ifdef PARALLEL
>         int            i, isrc, domainIdx, numPrecipitates,temp, newSize, reserved;
>         int            maxprecipitateTagIndex, numExpCells, iprecipitate, bufIndex;
>         int            iCell, cellIdx, cellPrecipitateCount, iNbr,iPNbr, numNbrs,numPNbrs;
>         real8          *inPBuf;
>         Precipitate_t         *precipitate;
>         Cell_t         *cell;
>         RemoteDomain_t *remDom;
>         
> /*
>  *      Loop through all the remote domains neighboring the current domain
>  *      and unpack the ghost precipitate data from that remote domain.
>  */
>   
>     
> 
>         for (isrc = 0; isrc < home->remoteDomainCount; isrc++) {
>             domainIdx = home->remoteDomains[isrc];
>             remDom = home->remoteDomainKeys [domainIdx];
>         
> /*
>  *          Free the old remDom->precipitateKeys table, if any
>  */
>             if (remDom->maxprecipitateTagIndex) {
>                 remDom->maxprecipitateTagIndex = 0;
>                 free(remDom->precipitateKeys);
>                 remDom->precipitateKeys = (Precipitate_t **)NULL;
>             }
>                  
> /*
>  *          First unpack the number of items in the buffer,
>  *          the size of the maximum tag index, and the number of cells
>  *          being sent.
>  */
>             inPBuf = (real8 *)remDom->inPBuf;
> 
>             bufIndex = 0;
> 
>             numPrecipitates            = inPBuf[bufIndex++];
>             reserved            = inPBuf[bufIndex++];
>             maxprecipitateTagIndex         = inPBuf[bufIndex++];
>             numExpCells         = inPBuf[bufIndex++];
>         
>             remDom-> maxprecipitateTagIndex =  maxprecipitateTagIndex;
> 
> /*
>  *          Allocate and initialize the precipitateKeys table for this
>  *          remote domain
>  */
>             newSize = maxprecipitateTagIndex * sizeof(Precipitate_t *);
>             remDom->precipitateKeys = (Precipitate_t **)calloc(1, newSize);
>         
> /*
>  *          Loop through the cells exported from this remote domain
>  */
>             for (iCell = 0; iCell < numExpCells; iCell++) {
>         
>                 cellIdx       = inPBuf[bufIndex++];
>                 //printf("UNpack precipitates %d %d %d \n",cellIdx,bufIndex,home->myDomain);
>                
>                 
>                 cellPrecipitateCount = inPBuf[bufIndex++];
>                 reserved      = inPBuf[bufIndex++];
> 
>                 cell = home->cellKeys[cellIdx];
> 				
>                 if (!cell) {
>                     int cX, cY, cZ;
>                     DecodeCellIdx(home, cellIdx, &cX, &cY, &cZ);
>                     Fatal("%s: received an unknown cell %d, (%d,%d,%d)",
>                           "CommUnpackPrecipitateGhosts", cellIdx, cX, cY, cZ);
>                 }
> 
>         
> /*
>  *              Loop through cell precipitates. For each precipitate obtain a free precipitate
>  *              structure, unpack the data into the structure, add it
>  *              to the ghost precipitate queue and the cell's precipitate queue.
>  */
>                 for (iprecipitate = 0; iprecipitate < cellPrecipitateCount; iprecipitate++) {
>         
>                     precipitate = PopFreePrecipitateQ(home);
> 
>                     precipitate->myTag.domainID   = domainIdx; /* known implicitly */
>                     precipitate->myTag.index      = inPBuf[bufIndex++];
>                     precipitate->numNbrs=inPBuf[bufIndex++];
>                      if ( home->myDomain==2  ) {
> 				//printf("Unpackprecipitate   %d cycle %d  tag (%d, %d) \n",home->myDomain,home->cycle,precipitate->myTag.domainID,precipitate->myTag.index ); 
>                     }         
396a739,791
>                     precipitate->x = inPBuf[bufIndex++];
>                     precipitate->y = inPBuf[bufIndex++];
>                     precipitate->z = inPBuf[bufIndex++];
>                     
>                    
>         
>                     precipitate->fX = inPBuf[bufIndex++];
>                     precipitate->fY = inPBuf[bufIndex++];
>                     precipitate->fZ = inPBuf[bufIndex++];
>         
>                     precipitate->vX = inPBuf[bufIndex++];
>                     precipitate->vY = inPBuf[bufIndex++];
>                     precipitate->vZ = inPBuf[bufIndex++];
> 					
> 					precipitate->forcep = inPBuf[bufIndex++];
>                     precipitate->r = inPBuf[bufIndex++];
>                    
>         
>         
>         
>                     precipitate->oldvX = inPBuf[bufIndex++];
>                     precipitate->oldvY = inPBuf[bufIndex++];
>                     precipitate->oldvZ = inPBuf[bufIndex++];
>         
>                     precipitate->constraint = (int)inPBuf[bufIndex++];
>                     precipitate->flags = (int)inPBuf[bufIndex++];
> #ifdef _FEM
>                     precipitate->fem_Surface[0] = (int)inPBuf[bufIndex++];
>                     precipitate->fem_Surface[1] = (int)inPBuf[bufIndex++];
>                     precipitate->fem_Surface_Norm[0] = inPBuf[bufIndex++];
>                     precipitate->fem_Surface_Norm[1] = inPBuf[bufIndex++];
>                     precipitate->fem_Surface_Norm[2] = inPBuf[bufIndex++];
> #endif
>                     precipitate->native = 0;
>         
> /*
>  *                  Register the precipitate in the remote domain's precipitateKeys array,
>  *                  move precipitate onto ghost queue, and add precipitate to the cell's
>  *                  precipitate queue.
>  */
>                     remDom->precipitateKeys[precipitate->myTag.index] = precipitate;
>                     PushGhostPrecipitateQ(home, precipitate);
>         
>                     precipitate->nextInCell = cell->precipitateQ;
>                     cell->precipitateQ = precipitate;
>                     cell->precipitateCount++;
>         
>                     precipitate->cellIdx = cellIdx;
> 
>                 }  /* end for (iprecipitate = 0; ...)  */
> 
>             }  /* end for (iCell = 0; ...) */
> //printf("unpack precipitates bufindex  %d %d \n",bufIndex,home->myDomain);
422c817
< 
---
>  
440c835
<                 Fatal("Missing rmeote domain struct!");
---
>                 Fatal("Missing remote domain struct!");
460c855,856
<         
---
> 			//printf("outbuf %p inBuf %p  \n",remDom->outBuf,remDom->inBuf);
> 
487c883,885
<         
---
>  
>   
>       
504a903,905
>         
>          
>         
573a975
> 	MPI_Barrier(MPI_COMM_WORLD);
589a992
> 
591a995,1206
> 
> 
> /*------------------------------------------------------------------------
>  *
>  *      Function:    CommSendPrecipitateGhosts
>  *      Description: Driver function to send nodal data for local precipitates
>  *                   to neighboring domains that export the precipitates as
>  *                   as ghosts, and to receive similar data for precipitates
>  *                   this domain maintains as ghosts.
>  *
>  *-----------------------------------------------------------------------*/
> void CommSendPrecipitateGhosts(Home_t *home) 
> {
> #ifdef PARALLEL
>         int            i, isrc, domainIdx, idst, idom, valCount;
>         int            localBuffers = 0;
>         int            remDomID, totRemDomCount;
>         RemoteDomain_t *remDom;
>         
>         //TimerStart(home, COMM_SEND_GHOSTS);
> 
> /*
>  *      All ghost precipitates (including secondary ghosts) have been
>  *      recycled, and the precipitateKeys arrays for the primary remote
>  *      domains will be reallocated as necessary when the primary
>  *      ghosts are received, but we have to free the precipitateKeys arrays
>  *      for the secondary remote domains here or we risk leaving
>  *      pointers to structures that have already been freed.
>  */
>         totRemDomCount = home->remoteDomainCount +
>                          home->secondaryRemoteDomainCount;
> 
>         for (i = home->remoteDomainCount; i < totRemDomCount; i++) {
> 			
>             remDomID = home->remoteDomains[i];
>             remDom = home->remoteDomainKeys[remDomID];
> 			
>             if (remDom == (RemoteDomain_t *)NULL) {
>                 Fatal("Missing remote domain struct!");
>             }
> 
>             if (remDom->maxprecipitateTagIndex) {
>                 free(remDom->precipitateKeys);
>                 remDom->precipitateKeys = (Precipitate_t **)NULL;
>                 //kommentoitu pois 1612016 (AL)
>                 //free(remDom);
>                 //home->remoteDomainKeys[remDomID] = (RemoteDomain_t *)NULL;
>             }
> //tama pitää olla pois jotta koodi toimii,sotkee remotedomainit muuten
>             //home->secondaryRemoteDomainCount--;
>         }
>         
>   
> 
> 
> 
> /*
>  *      Pre-issue receives of message lengths from each neighbor
>  */
>         for (isrc = 0; isrc < home->remoteDomainCount; isrc++) {
>         
>             domainIdx = home->remoteDomains[isrc];
>             remDom = home->remoteDomainKeys[domainIdx];
>            
> 
>             MPI_Irecv(&remDom->inPBufLen, 1, MPI_INT, domainIdx, MSG_GHOST_PR_LEN ,
>                       MPI_COMM_WORLD, &home->inPRequests[isrc]);
>         }
>         
>      
>         
>       
>         
> /*
>  *      Package up precipitate data for neighboring domains and send
>  *      out the buffer lengths
>  */
>         CommPackPrecipitateGhosts(home);
>         
>         for (idst = 0; idst < home->remoteDomainCount; idst++) {
>         
>             domainIdx = home->remoteDomains[idst];
>             remDom = home->remoteDomainKeys[domainIdx];
>         
>             MPI_Isend(&remDom->outPBufLen, 1, MPI_INT, domainIdx, MSG_GHOST_PR_LEN,
>                       MPI_COMM_WORLD, &home->outPRequests[idst]);
> 
>             localBuffers += remDom->outPBufLen;
>         }
> 
> /*
>  *      Wait for the length sends/receives to complete
>  */
>         MPI_Waitall(home->remoteDomainCount, home->outPRequests,home->outPStatus);
>         MPI_Waitall(home->remoteDomainCount, home->inPRequests,home->inPStatus);
> 
>     
>       
> /*
>  *      Allocate appropriately sized buffers for the incoming messages
>  *      and pre-issue receives for buffers from all neighboring domains
>  */
>         for (isrc = 0; isrc < home->remoteDomainCount; isrc++) {
>         
>             domainIdx = home->remoteDomains[isrc];
>             remDom = home->remoteDomainKeys[domainIdx];
>            
>             valCount = remDom->inPBufLen / sizeof(real8);
>             
>             remDom->inPBuf = (char *)malloc(remDom->inPBufLen);
>             MPI_Irecv(remDom->inPBuf, valCount, MPI_DOUBLE, domainIdx, 
>                       MSG_GHOST_PR, MPI_COMM_WORLD, &home->inPRequests[isrc]);
> 
>             localBuffers += remDom->inPBufLen;
>         }
>         
>         
>          
>          
> /*
>  *      Send local data out to all the neighboring domains
>  */
>         for (idst = 0; idst < home->remoteDomainCount; idst++) {
>         
>             domainIdx = home->remoteDomains[idst];
>             remDom = home->remoteDomainKeys[domainIdx];
>         
>             valCount = remDom->outPBufLen / sizeof(real8);
>             MPI_Isend(remDom->outPBuf, valCount, MPI_DOUBLE, domainIdx, 
>                       MSG_GHOST_PR, MPI_COMM_WORLD, &home->outPRequests[idst]);
>         }
>         
>         
> /*
>  *      Wait for all data buffer sends/receives to complete and unpack
>  *      the data
>  */
>         MPI_Waitall(home->remoteDomainCount, home->outPRequests,home->outPStatus);
>         MPI_Waitall(home->remoteDomainCount, home->inPRequests,home->inPStatus);
>          
>         CommUnpackPrecipitateGhosts(home);
>         MPI_Barrier(MPI_COMM_WORLD);
> /*
>  *      Just some debug code for printing the maximum buffer space
>  *      used by any domain during the ghost precipitate communications.
>  */
> #if 0
> {
>         int globalBuffers = 0;
> 
>         MPI_Allreduce(&localBuffers, &globalBuffers, 1, MPI_INT, MPI_MAX,
>                       MPI_COMM_WORLD);
> 
>         if (globalBuffers == localBuffers) {
>             printf("  Task %d: Ghost comm total buffers = %dKb\n",
>                    home->myDomain, globalBuffers / 1000);
>         }
> }
> #endif
> 
> /*
>  *      Release the message buffers...
>  */
>         for (idom = 0; idom < home->remoteDomainCount; idom++) {
>         
>             domainIdx = home->remoteDomains[idom];
>             remDom = home->remoteDomainKeys[domainIdx];
>          
>             free(remDom->inPBuf);
>             free(remDom->outPBuf);
> 
>             remDom->inPBuf = (char *)NULL;
>             remDom->outPBuf = (char *)NULL;
> 
>             remDom->inPBufLen = 0;
>             remDom->outPBufLen = 0;
>         }
> 
> /*
>  *      Turns out, in order to be sure we do all the necessary
>  *      direct segment-to-segment force interactions, each domain
>  *      needs a layer of secondary ghosts which are nodes outside
>  *      any native cell or cell immediately neighboring a native
>  *      one, but connected to a primary ghost (i.e. a ghost within
>  *      a native or immediately adjoining cell).
>  *
>  *      So, we have to do one more communication to get those
>  *      secondary ghosts.
>  */
>         
> 
>         //TimerStop(home, COMM_SEND_GHOSTS);
> 
> #endif  /* if PARALLEL */
> 
> /*
>  *      Measure dead time after ghost precipitate communications.
>  */
> #if PARALLEL
> #ifdef SYNC_TIMERS
>         TimerStart(home, GHOST_COMM_BARRIER);
>         MPI_Barrier(MPI_COMM_WORLD);
>         TimerStop(home, GHOST_COMM_BARRIER);
> #endif
> #endif
> 
>         return;
> }
> 
> 
> 
> 
