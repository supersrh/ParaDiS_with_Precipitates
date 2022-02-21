30c30,32
< #define MSG_MIG_NODES  3001
---
> #define MSG_MIGP_LENGTH 3001
> #define MSG_MIG_NODES  3002
> #define MSG_MIG_PRECIPITATES  3003
38a41
> #define VALS_PER_MIG_PRECIPITATE  14
40c43,44
< #define VALS_PER_MIG_NODE  14
---
> #define VALS_PER_MIG_NODE  18 /*originally 14 added 4 (AL)*/
> #define VALS_PER_MIG_PRECIPITATE  14
43a48
> #define VALS_PER_MIG_PNRB       2
163a169,245
> static void BuildPrecipitateMigLists(Home_t *home, int *migCommList, int **migList,
>                           int *migListCnt, int **sendBufDest, int *numSendBufs)
> {
>         int          i, j;
>         int          precipitateCount, numDomains, thisDomain, destDom;
>         int          allocatedValues, sendBufCnt, listLen;
>         int          *sendBufList;
>         Precipitate_t       *precipitate;
> 
>         precipitateCount = home->newPrecipitateKeyPtr;
>         numDomains = home->numDomains;
>         thisDomain = home->myDomain;
> 
>         sendBufCnt = 0;
>         allocatedValues = 0;
>         sendBufList = (int *)NULL;
> 
>         for (i = 0; i < precipitateCount; i++) {
> 
>             if ((precipitate = home->precipitateKeys[i]) == (Precipitate_t *)NULL) {
>                 continue;
>             }
> 
> /*
>  *          Search the domain decomposition for the domain which
>  *          should own this precipitate.  If the current domain should retain
>  *          ownership, move on to the next precipitate.
>  */
>             destDom = FindCoordDomain(home, 1, &precipitate->x, &precipitate->y, &precipitate->z);
> 
>             if (destDom == thisDomain) {
>                 continue;
>             }
> /*
>  *          For any precipitate that needs migrating, add the precipitate to the list
>  *          of precipitates being sent to that particular remote domain, increment
>  *          the count of precipitates beng sent there.  Also set the flag for that
>  *          remote domain in the list of domains to which this domain
>  *          will need to migrate precipitates.
>  */
>             migCommList[destDom] = 1;
>             migList[destDom] = (int *)realloc(migList[destDom],
>                                               (migListCnt[destDom]+1) *
>                                               sizeof(int));
>             migList[destDom][migListCnt[destDom]] = i;
>             migListCnt[destDom] += 1;
> 
> /*
>  *          If the remote domain is not already on the list of domains
>  *          to which this domain will be migrating stuff, add the remote
>  *          domain ID to the list... so we don't have to loop through
>  *          the entire list of remote domains later looking for the
>  *          ones to which we're migrating stuff.
>  */
>             for (j = 0; j < sendBufCnt; j++) {
>                 if (destDom == sendBufList[j]) {
>                     break;
>                 }
>             }
> 
>             if (j == sendBufCnt) {
>                 sendBufCnt += 1;
>                 if (sendBufCnt > allocatedValues) {
>                     allocatedValues += 25;
>                     sendBufList = (int *)realloc(sendBufList, allocatedValues *
>                                                  sizeof(int));
>                 }
>                 sendBufList[sendBufCnt-1] = destDom;
>             }
>             
>         }
> 
>         *sendBufDest = sendBufList;
>         *numSendBufs = sendBufCnt;
> 
>         return;
> }
193c275
<         int          migNodeCount, armCount, numVals;
---
>         int          migNodeCount, armCount,precipitateCount, numVals;
197c279,280
< 
---
> 		
> 		precipitateCount=0;
206a290
> 
212a297
>             precipitateCount += home->nodeKeys[migList[i]]->numPNbrs;
216a302
>                   precipitateCount*VALS_PER_MIG_PNRB +
233a320
>             buf[bufIndex++] = (real8)node->numPNbrs;
239a327,330
>             
>             buf[bufIndex++] = node->nodedx;
>             buf[bufIndex++] = node->nodedy;
>             buf[bufIndex++] = node->nodedz;
273a365,377
>             
>             for (j = 0; j < node->numPNbrs; j++) {
> 
>                 buf[bufIndex++] = (real8)node->PnbrTag[j].domainID;
>                 buf[bufIndex++] = (real8)node->PnbrTag[j].index;
> 
>                 
>             }
>             
>             
>             
>             
>             
277a382
> 			//FreeNodePNbrs(node);
291a397,498
> static void PackPrecipitateMigrators(Home_t *home, int remDomID, int *migList,
>                           int *migCounts, real8 **sendBuf, int *sendBufLen)
> {
>         int          i, j;
>         int          migPrecipitateCount, armCount, numVals;
>         int          bufIndex;
>         real8        *buf;
>         Precipitate_t       *precipitate;
> 
>         armCount = 0;
>         bufIndex = 0;
> 
> /*
>  *      migCounts[remDomID] contains the number of precipitates being
>  *      migrated to the remote domain <remDomID.
>  */
>         migPrecipitateCount = migCounts[remDomID];
> 
> /*
>  *      Determine how large a buffer is required to hold all the
>  *      specified entities and allocate an appropriately sized buffer.
>  */
>         for (i = 0; i < migPrecipitateCount; i++) {
>             armCount += home->precipitateKeys[migList[i]]->numNbrs;
>         }
> 
>         numVals = migPrecipitateCount * VALS_PER_MIG_PRECIPITATE +
>                  +
>                   VALS_PER_MIG_EXTRA;
> 
>         buf = (real8 *)malloc(numVals * sizeof(real8));
> 
>         buf[bufIndex++] = (real8)migPrecipitateCount;
>         buf[bufIndex++] = 0.0;  /* reserved */
> 
> /*
>  *      Loop through all the precipitates to be sent and pack the necessary
>  *      nodal data into the buffer.
>  */
>         for (i = 0; i < migPrecipitateCount; i++) {
>             precipitate = home->precipitateKeys[migList[i]];
> 
>             buf[bufIndex++] = (real8)precipitate->myTag.index;
>             buf[bufIndex++] = (real8)precipitate->forcep;
>             buf[bufIndex++] = (real8)precipitate->r;
>             buf[bufIndex++] = (real8)precipitate->sgnv;
>             buf[bufIndex++] = (real8)precipitate->flags;
> 
>             buf[bufIndex++] = precipitate->x;
>             buf[bufIndex++] = precipitate->y;
>             buf[bufIndex++] = precipitate->z;
> 
>             buf[bufIndex++] = precipitate->vX;
>             buf[bufIndex++] = precipitate->vY;
>             buf[bufIndex++] = precipitate->vZ;
> 
>             buf[bufIndex++] = precipitate->oldvX;
>             buf[bufIndex++] = precipitate->oldvY;
>             buf[bufIndex++] = precipitate->oldvZ;
> 
> #ifdef _FEM
>             buf[bufIndex++] = precipitate->fem_Surface[0];
>             buf[bufIndex++] = precipitate->fem_Surface[1];
> 
>             buf[bufIndex++] = precipitate->fem_Surface_Norm[0];
>             buf[bufIndex++] = precipitate->fem_Surface_Norm[1];
>             buf[bufIndex++] = precipitate->fem_Surface_Norm[2];
> #endif
>             //for (j = 0; j < precipitate->numNbrs; j++) {
> 
>                 //buf[bufIndex++] = (real8)precipitate->nbrTag[j].domainID;
>                 //buf[bufIndex++] = (real8)precipitate->nbrTag[j].index;
> 
>                 //buf[bufIndex++] = precipitate->burgX[j];
>                 //buf[bufIndex++] = precipitate->burgY[j];
>                 //buf[bufIndex++] = precipitate->burgZ[j];
> 
>                 //buf[bufIndex++] = precipitate->nx[j];
>                 //buf[bufIndex++] = precipitate->ny[j];
>                 //buf[bufIndex++] = precipitate->nz[j];
> 
>                 //buf[bufIndex++] = precipitate->armfx[j];
>                 //buf[bufIndex++] = precipitate->armfy[j];
>                 //buf[bufIndex++] = precipitate->armfz[j];
>             //}
> 
> /*
>  *          Free the precipitate and recyle the tag
>  */
>             FreePrecipitate(home, migList[i]);
> 
>         }  /* for (i = 0; i < migCount; ...) */
> 
> /*
>  *      Return pointer to buffer and the buffer size to the caller.
>  *      Caller is responsible for freeing the buffer.
>  */
>         *sendBuf = buf;
>         *sendBufLen = numVals * sizeof(real8);
> 
>         return;
> }
307c514
<         int          i, j, bufIndex, nodeIndex, nodeCount, numNbrs;
---
>         int          i, j, bufIndex, nodeIndex, nodeCount, numNbrs,numPNbrs;
343a551,552
>             numPNbrs           = (int)buf[bufIndex++];
>            
349a559,562
>             
>             node->nodedx = buf[bufIndex++];
>             node->nodedy = buf[bufIndex++];
>             node->nodedz = buf[bufIndex++];
371a585
>             AllocNodePrecipitates(node, numPNbrs);
389a604,616
>             
>              for (j = 0; j < numPNbrs; j++) {
> 
>                 node->PnbrTag[j].domainID = (int)buf[bufIndex++];
>                 node->PnbrTag[j].index = (int)buf[bufIndex++];
> 
>                
>             }
>             
>             
>             
>             
>             
394a622,631
> static void UnpackPrecipitateMigrators(Home_t *home, real8 *buf, int remDomID)
> {
>         int          i, j, bufIndex, precipitateIndex, precipitateCount, numNbrs;
>         int          thisDomain;
>         int          newSize, newIndex;
>         Tag_t        oldTag;
>         Precipitate_t       *precipitate;
> 
>         thisDomain = home->myDomain;
>         bufIndex = 0;
395a633,715
> /*
>  *      Get the precipitate count from the first element of the buffer.
>  *      The second element is reserved,so skip it.
>  */
>         precipitateCount = (int)buf[bufIndex++];
>         bufIndex++;
> 
>         for (i = 0; i < precipitateCount; i++) {
> /*
>  *          Add a new precipitate to the list of local precipitates and populate the
>  *          precipitate structure with data from the remote domain.  Also
>  *          add a mapping between the original precipitate tag from the
>  *          remote domain and the new precipitate tag in the local domain.
>  */
>  
>             precipitateIndex = GetFreePrecipitateTag(home);
>    
>             precipitate = PopFreePrecipitateQ(home);
>            
>             home->precipitateKeys[precipitateIndex] = precipitate;
> 		
> 
> 			
>             precipitate->myTag.domainID = thisDomain;
>             precipitate->myTag.index    = precipitateIndex;
> 
>             oldTag.domainID = remDomID;
>             oldTag.index    = (int)buf[bufIndex++];
> 			//printf("Unpackprecipitates,(%d %d),(%d %d)\n ",thisDomain,precipitateIndex,oldTag.domainID,oldTag.index);
>             AddPrecipitateTagMapping(home, &oldTag, &precipitate->myTag);
> 			
>             precipitate->forcep  = buf[bufIndex++];
>             precipitate->r           =  buf[bufIndex++];
>             precipitate->sgnv        = (int)buf[bufIndex++];
>             precipitate->flags       = (int)buf[bufIndex++];
> 
>             precipitate->x = buf[bufIndex++];
>             precipitate->y = buf[bufIndex++];
>             precipitate->z = buf[bufIndex++];
> 
>             precipitate->vX = buf[bufIndex++];
>             precipitate->vY = buf[bufIndex++];
>             precipitate->vZ = buf[bufIndex++];
> 
>             precipitate->oldvX = buf[bufIndex++];
>             precipitate->oldvY = buf[bufIndex++];
>             precipitate->oldvZ = buf[bufIndex++];
> 
> #ifdef _FEM
>             precipitate->fem_Surface[0] = buf[bufIndex++];
>             precipitate->fem_Surface[1] = buf[bufIndex++];
> 
>             precipitate->fem_Surface_Norm[0] = buf[bufIndex++];
>             precipitate->fem_Surface_Norm[1] = buf[bufIndex++];
>             precipitate->fem_Surface_Norm[2] = buf[bufIndex++];
> #endif
> 
> /*
>  *          Set all the segment specific nodal data values.
>  */
>             //AllocPrecipitateArms(precipitate, numNbrs);
> 
>             //for (j = 0; j < numNbrs; j++) {
> 
>                 //precipitate->nbrTag[j].domainID = (int)buf[bufIndex++];
>                 //precipitate->nbrTag[j].index = (int)buf[bufIndex++];
> 
>                 //precipitate->burgX[j] = buf[bufIndex++];
>                 //precipitate->burgY[j] = buf[bufIndex++];
>                 //precipitate->burgZ[j] = buf[bufIndex++];
> 
>                 //precipitate->nx[j] = buf[bufIndex++];
>                 //precipitate->ny[j] = buf[bufIndex++];
>                 //precipitate->nz[j] = buf[bufIndex++];
> 
>                 //precipitate->armfx[j] = buf[bufIndex++];
>                 //precipitate->armfy[j] = buf[bufIndex++];
>                 //precipitate->armfz[j] = buf[bufIndex++];
>             //}
>         }
> 
>         return;
> }
572a893,1053
> 
> static void CommSendPrecipitateMigrators(Home_t *home, int *migCommList, int **migList,
>                               int *migListCnt, int *sendBufDest,
>                               int numSendBufs)
> {
>         int         i, numValues, numRecvBufs, recvIndex;
>         int         numDomains, thisDomain, domID;
>         int         *glblMigCommList;
>         int         *sendBufLen, *recvBufLen;
>         char        **sendBuf, **recvBuf;
>         MPI_Request *sendRequest, *recvRequest;
>         MPI_Status  *sendStatus, *recvStatus;
>         
> 
>         numDomains = home->numDomains;
>         thisDomain = home->myDomain;
> 
> /*
>  *      First we need to determine how many remote domains will be
>  *      migrating precipitates to this domain.  Each domain has already set
>  *      the flag in the migCommList for all domains to which it will
>  *      migrate precipitates.  When the all-reduce is done, each domain
>  *      will know how many other domains will be migrating precipitates
>  *      to it...
>  */
>         glblMigCommList = (int *)calloc(1, numDomains * sizeof(int));
> 
>         MPI_Allreduce(migCommList, glblMigCommList, numDomains, MPI_INT,
>                       MPI_SUM, MPI_COMM_WORLD);
> 
>         numRecvBufs = glblMigCommList[thisDomain];
> 
>         free(glblMigCommList);
> 
> /*
>  *      Pack buffers for each remote domain to which this domain
>  *      will be migrating precipitates
>  */
>         if (numSendBufs > 0) {
>             sendBuf = (char **)calloc(1, numSendBufs * sizeof(real8 *));
>             sendBufLen = (int *)calloc(1, numSendBufs * sizeof(int));
>             sendRequest = (MPI_Request *)malloc(numSendBufs *
>                                                 sizeof(MPI_Request));
>             sendStatus = (MPI_Status *)malloc(numSendBufs *
>                                               sizeof(MPI_Status));
>         }
> 
>         for (i = 0; i < numSendBufs; i++) {
>             domID = sendBufDest[i];
>             PackPrecipitateMigrators(home, domID, migList[domID], migListCnt,
>                           (real8 **)&sendBuf[i], &sendBufLen[i]);
>         }
> 
> /*
>  *      Allocate arrays for handling incoming migrated precipitates
>  */
>         if (numRecvBufs > 0) {
>             recvBuf = (char **)calloc(1, numRecvBufs * sizeof(real8 *));
>             recvBufLen = (int *)calloc(1, numRecvBufs * sizeof(int));
>             recvRequest = (MPI_Request *)malloc(numRecvBufs *
>                                                 sizeof(MPI_Request));
>             recvStatus = (MPI_Status *)malloc(numRecvBufs *
>                                               sizeof(MPI_Status));
>         }
> 
> /*
>  *      Pre-issue receives of buffer lengths from all domains that
>  *      will be be migrating precipitates to this domain.  Lengths are
>  *      specified in units of bytes.
>  */
>         for (i = 0; i < numRecvBufs; i++) {
>             MPI_Irecv(&recvBufLen[i], 1, MPI_INT, MPI_ANY_SOURCE,
>                       MSG_MIGP_LENGTH, MPI_COMM_WORLD, &recvRequest[i]);
>         }
> 
> /*
>  *      Have the current domain send out the sizes of buffers it will
>  *      be transmitting.
>  */
>         for (i = 0; i < numSendBufs; i++) {
>             MPI_Isend(&sendBufLen[i], 1, MPI_INT, sendBufDest[i],
>                       MSG_MIGP_LENGTH, MPI_COMM_WORLD, &sendRequest[i]);
>         }
> 
> /*
>  *      Wait for all length send/receives to complete
>  */
>         MPI_Waitall(numSendBufs, sendRequest, sendStatus);
>         MPI_Waitall(numRecvBufs, recvRequest, recvStatus);
> 
> /*
>  *      Allocate receive buffers of the appropriate sizes and post the
>  *      receives associated with those buffers.
>  */
>         for (i = 0; i < numRecvBufs; i++) {
>             recvBuf[i] = (char *)malloc(recvBufLen[i]);
>             numValues = recvBufLen[i] / sizeof(real8);
>             MPI_Irecv(recvBuf[i], numValues, MPI_DOUBLE,
>                       recvStatus[i].MPI_SOURCE, MSG_MIG_PRECIPITATES,
>                       MPI_COMM_WORLD, &recvRequest[i]);
>         }
> 
> /*
>  *      Send all the migrating precipitates to the appropriate remote
>  *      domains.
>  */
>         for (i = 0; i < numSendBufs; i++) {
>             numValues = sendBufLen[i] / sizeof(real8);
>             MPI_Isend(sendBuf[i], numValues, MPI_DOUBLE, sendBufDest[i],
>                       MSG_MIG_PRECIPITATES, MPI_COMM_WORLD, &sendRequest[i]);
>         }
> 
> /*
>  *      Process the incoming buffers as soon as they arrive.  Status
>  *      will be placed in first recvStatus struct each time.
>  */
>  
> 
>         for (i = 0; i < numRecvBufs; i++) {
>             MPI_Waitany(numRecvBufs, recvRequest, &recvIndex,
>                         &recvStatus[0]);
>             UnpackPrecipitateMigrators(home, (real8 *)recvBuf[recvIndex],
>                             recvStatus[0].MPI_SOURCE);
>             free(recvBuf[recvIndex]);
>         }
> 
> 
> 
> /*
>  *      Wait for all buffer sends to complete.
>  */
>         MPI_Waitall(numSendBufs, sendRequest, sendStatus);
> 
> /*
>  *      Free all the temporary buffers before returning to the caller.
>  */
>         if (numSendBufs > 0) {
>             for (i = 0; i < numSendBufs; i++) {
>                 free(sendBuf[i]);
>             }
>             free(sendBuf);
>             free(sendBufLen);
>             free(sendRequest);
>             free(sendStatus);
>         }
> 
>         if (numRecvBufs > 0) {
>             free(recvBuf);
>             free(recvBufLen);
>             free(recvRequest);
>             free(recvStatus);
>         }
> 
>         return;
> }
> 
> 
> 
> 
> 
> 
643c1124
< 
---
> //MPI_Barrier(MPI_COMM_WORLD);
646a1128,1130
>  
>  
>  
654a1139
>          
662a1148,1248
>  */
>         TimerStart(home, MIGRATION_BARRIER);
>         MPI_Barrier(MPI_COMM_WORLD);
>         TimerStop(home, MIGRATION_BARRIER);
> #endif
> #endif
> 
>         return;
> }
> 
> 
> 
> 
> 
> 
> void MigratePrecipitate(Home_t *home)
> {
>         int    i, numDomains, numSendBufs;
>         int    *migCommList, **migList, *migListCnt;
>         int    *sendBufDest;
> 
>         TimerStart(home, MIGRATION);
> 		
> #ifdef PARALLEL
> /*
>  *      Precipitates will only migrate if we're running in parallel
>  */
> 
>         numDomains = home->numDomains;
> 
>         numSendBufs = 0;
>         sendBufDest = (int *)NULL;
> 
>         migCommList = (int *)calloc(1, numDomains * sizeof(int));
>         migList = (int **)calloc(1, numDomains * sizeof(int *));
> 
>         migListCnt = (int *)calloc(1, numDomains * sizeof(int));
> 
> /*
>  *      Look through all local precipitates and determine which precipitates need
>  *      to be migrated and the domains to which those precipitates must
>  *      be migrated.  For each remote domain to which this domain
>  *      will migrate one or more precipitates, build a list of the precipitates
>  *      to be sent to that domain.
>  */
>         BuildPrecipitateMigLists(home, migCommList, migList, migListCnt,
>                       &sendBufDest, &numSendBufs);
> 
> #if 0
>         for (i = 0; i < numDomains; i++) {
>             if (migListCnt[i]) {
>                 printf("Task %d migrate %d precipitates to %d\n", home->myDomain,
>                        migListCnt[i], i);
>             }
>         }
> #endif
> 
> /*
>  *      Send out all precipitates (if any) that need to be migrated
>  *      to remote domains and receive any precipitates migrating from
>  *      other domains.
>  */
>         CommSendPrecipitateMigrators(home, migCommList, migList, migListCnt,
>                           sendBufDest, numSendBufs);
>                           
>                           
> 
> /*
>  *      All migrated precipitates have been retagged.  Each domain now needs
>  *      to communicate to its neighboring domains the mappings between
>  *      the old and new tags for all precipitates it received during the
>  *      migration.  Once that is done, the local domains go through
>  *      all their own precipitates reconciling the precipitate tag changes.
>  */
>         DistributePrecipitateTagMaps(home);
>         MPI_Barrier(MPI_COMM_WORLD);
>       
>         
> 
> /*
>  *      Free up all temporary arrays before returning to the caller.
>  */
>         for (i = 0; i < numSendBufs; i++) {
>             free(migList[sendBufDest[i]]);
>         }
> 
>         free(migList);
>         free(migListCnt);
>         free(migCommList);
>         free(sendBufDest);
>         
>         
>          
> #endif  /* ifdef PARALLEL */
> 
>         TimerStop(home, MIGRATION);
> 
> #if PARALLEL
> #ifdef SYNC_TIMERS
> /*
>  *      Measure dead time after precipitate migration
