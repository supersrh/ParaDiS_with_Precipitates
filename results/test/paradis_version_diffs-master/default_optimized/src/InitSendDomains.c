20a21
> #include "InPrecipitateData.h"
23a25
> #include "Precipitate.h"
56c58
<         int    i, armCount, armID, numNbrs;
---
>         int    i, armCount,nbrprecipitateCount,armID,pnbrID, numNbrs,numPNbrs;
69a72,79
> 		
> 		nbrprecipitateCount=0;
> 		 for (i = 0; i < nodeCount; i++) {
>             nbrprecipitateCount += inData->node[nodeList[i]].numPNbrs;
>         }
> 
> 
> 
72c82
<                    (armCount * INIT_VALS_PER_ARM);
---
>                    (armCount * INIT_VALS_PER_ARM)+(nbrprecipitateCount*INIT_VALS_PER_PNRB);
86a97
>         
95a107
>             numPNbrs = inData->node[nIndex].numPNbrs;
96a109
>             dataBuf[bufIndex++] = (real8)numPNbrs;
99c112
< 
---
> 				
119a133,150
>             
>             
>             for (pnbrID = 0; pnbrID < numPNbrs; pnbrID++) {
> 				
> /*
>  *              Tag of the node at the terminating end of the arm.
>  */
>                 dataBuf[bufIndex++] =
>                         (real8)inData->node[nIndex].PnbrTag[pnbrID].domainID;
>                 dataBuf[bufIndex++] =
>                         (real8)inData->node[nIndex].PnbrTag[pnbrID].index;
> 
>             }
>             
>             
>             
>             
>             
125a157,160
>             
>             dataBuf[bufIndex++] = inData->node[nIndex].nodedx; 
>             dataBuf[bufIndex++] = inData->node[nIndex].nodedy; 
>             dataBuf[bufIndex++] = inData->node[nIndex].nodedz;
145c180
<  *      Function:     AddTagMapping
---
>  *      Function:     AddTagMapping for nodes
173a209,242
> void AddPrecipitateTagMapping(Home_t *home, Tag_t *oldTag, Tag_t *newTag)
> {
>         int      newSize;
>         TagMap_t *mapping;
> 
>         if (home->precipitatetagMapEnts >= home->precipitatetagMapSize) {
>               home->precipitatetagMapSize += NEW_PRECIPITATEKEY_INC;
>               newSize = home->precipitatetagMapSize * sizeof(TagMap_t);
>               home->precipitatetagMap = (TagMap_t *)realloc(home->precipitatetagMap, newSize);
>         }
> 
>         mapping = &home->precipitatetagMap[home->precipitatetagMapEnts];
> 
>         mapping->oldTag.domainID = oldTag->domainID;
>         mapping->oldTag.index    = oldTag->index;
> 
>         mapping->newTag.domainID = newTag->domainID;
>         mapping->newTag.index    = newTag->index;
> 
>         home->precipitatetagMapEnts++;
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
> 
> 
196a266,298
> /*---------------------------------------------------------------------------
>  *
>  *      Function:     ExtendPrecipitateKeys
>  *      Description:  
>  *
>  *-------------------------------------------------------------------------*/
> static void ExtendPrecipitateKeys(Home_t *home, int newLength)
> {
>         int i, oldLength;
> 
>         oldLength = home->newPrecipitateKeyMax;
>         home->newPrecipitateKeyMax = newLength;
> 
>         home->precipitateKeys = (Precipitate_t **)realloc(home->precipitateKeys,
>                                             newLength * sizeof(Precipitate_t *));
> 
>         for (i = oldLength; i < newLength; i++) {
>             home->precipitateKeys[i] = (Precipitate_t *)NULL;
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
> 
> 
> 
> 
> 
288a391,484
> static void GetNewPrecipitateTag(Home_t *home, Tag_t *oldTag, Tag_t *newTag,
>                       int *nextAvailableTag)
> {
>         int     nextTag, thisDomain;
> 
>         nextTag    = *nextAvailableTag;
>         thisDomain = home->myDomain;
> 
> /*
>  *      If the old tag belonged to a different domain, we must
>  *      give the precipitate a new tag.
>  */
>         if (oldTag->domainID != thisDomain) {
> 
>             for ( ; ; nextTag++) {
> 
> /*
>  *              Extend the precipitatekeys array if necessary.
>  */
>                 if (nextTag >= home->newPrecipitateKeyMax) {
>                     ExtendPrecipitateKeys(home, home->newPrecipitateKeyMax + NEW_PRECIPITATEKEY_INC);
>                 }
> 
>                 if (home->precipitateKeys[nextTag] == (Precipitate_t *)NULL) {
>                     newTag->domainID = thisDomain;
>                     newTag->index = nextTag++;
>                     AddPrecipitateTagMapping(home, oldTag, newTag);
>                     break;
>                 }
>             }
>         } else {
> /*
>  *          The old tag belonged to this domain...  check if that tag
>  *          is available for this run.
>  *
>  *          Extend the precipitatekeys array if necessary.
>  */
>             if (oldTag->index >= home->newPrecipitateKeyMax) {
>                 ExtendPrecipitateKeys(home, oldTag->index + NEW_PRECIPITATEKEY_INC);
>             }
> 
> /*
>  *          If the old tag is still available, use it and return
>  *          to the caller.
>  */
>             if (home->precipitateKeys[oldTag->index] == (Precipitate_t *)NULL) {
>                 newTag->domainID = oldTag->domainID;
>                 newTag->index    = oldTag->index;
>                 if (newTag->index >= home->newPrecipitateKeyPtr) {
>                     home->newPrecipitateKeyPtr = newTag->index + 1;
>                 }
>                 return;
>             }
> 
> /*
>  *          The old tag is no longer available, so just use
>  *          the next available tag.
>  */
>             for ( ; ; nextTag++) {
> 
> /*
>  *              Extend precipitate keys array if necessary
>  */
>                 if (nextTag >= home->newPrecipitateKeyMax) {
>                     ExtendPrecipitateKeys(home, home->newPrecipitateKeyMax + NEW_PRECIPITATEKEY_INC);
>                 }
> 
>                 if (home->precipitateKeys[nextTag] == (Precipitate_t *)NULL) {
>                     newTag->domainID = thisDomain;
>                     newTag->index = nextTag++;
>                     AddPrecipitateTagMapping(home, oldTag, newTag);
>                     break;
>                 }
>             }
>         }
> 
>         if (newTag->index >= home->newPrecipitateKeyPtr) {
>             home->newPrecipitateKeyPtr = newTag->index + 1;
>         }
> 
>         *nextAvailableTag = nextTag;
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
> 
307c503
<         int        i, armID, numNbrs, nodesInBuf, bufIndex;
---
>         int        i, armID,pnbrID , numNbrs,numPNbrs, nodesInBuf, bufIndex;
338a535
>             
341c538,539
< 
---
>             numPNbrs = (int)buf[bufIndex++];
> 			
343c541,542
< 
---
> 			AllocNodePrecipitates(node,numPNbrs);
> 			
347a547,549
>                 
>               
>                 
356a559,567
>             
>              for (pnbrID = 0; pnbrID < numPNbrs; pnbrID++) {
> 				 
> 				node->PnbrTag[pnbrID].domainID = (int)buf[bufIndex++];
>                 node->PnbrTag[pnbrID].index    = (int)buf[bufIndex++];
> 				//printf("Unpackinitial %d %d \n ",node->PnbrTag[armID].domainID, node->PnbrTag[armID].index );
> 		}
>             
>             
362a574,577
>             
> 			node->nodedx = buf[bufIndex++]; 
>             node->nodedy = buf[bufIndex++]; 
>             node->nodedz = buf[bufIndex++]; 
600a816,1398
> 
> /*---------------------------------------------------------------------------
>  *
>  *      Function:     PackInitialPrecipitateData
>  *      Description:  Pack all data from the current block of precipitates that is
>  *                    relevant to the specified domain into a buffer for
>  *                    trasnmission to a remote domain.  Data includes generic
>  *                    parameters, and precipitate info (precipitate tags of its neighbor
>  *                    precipitates, burgers vector for precipitate, glide plane orientation,
>  *                    coordinates, constraint type).
>  *      Arguments:
>  *          inprecipitateData      structure with initial nodal data
>  *          precipitateList    Integer array of indices into the inprecipitateData precipitate
>  *                      list for precipitates to be packed into the buffer
>  *          precipitateCount   Number of precipitate indices in the <precipitateList> array.
>  *          buf         pointer to location in which to return to the
>  *                      caller a pointer to the buffer of data to
>  *                      be sent to the remote domain
>  *          bufEnts     pointer to the location in which to return
>  *                      to the caller the number of values packed into
>  *                      the buffer
>  *
>  *-------------------------------------------------------------------------*/
> static void PackInitialPrecipitateData(InPrecipitateData_t *inprecipitateData, int *precipitateList,
>                                 int precipitateCount, real8 **buf, int *bufEnts)
> {
>         int    i, armCount, armID, numNbrs;
>         int    bufIndex, maxIndex, nIndex;
>         real8  *dataBuf,forcep,r;
>    
> /*
>  *      To calculate the buffer size we first need a count of all
>  *      the arms associated with all precipitates to be sent to the remote
>  *      domain.  
>  */
>         armCount = 0;
> 
>         for (i = 0; i < precipitateCount; i++) {
>             armCount += inprecipitateData->precipitate[precipitateList[i]].numNbrs;
>         }
> 
>         maxIndex = INIT_PARAMS + (precipitateCount * INIT_VALS_PER_PRECIPITATE) +
>                    (armCount * INIT_VALS_PER_ARM);
> 
>         dataBuf = (real8 *) malloc(maxIndex * sizeof(real8));
> 
> /*
>  *      Pack the precipitate counts as the first items in the buffer
>  *      then pack all the precipitates data.
>  *
>  *      Note: The remote domain will retag the precipitates it receives
>  *      as necessary and reconcile the old/new tags with neighboring
>  *      domains later.
>  */
>         bufIndex = 0;
> 
>         dataBuf[bufIndex++] = (real8)precipitateCount;
>    
>         for (i = 0; i < precipitateCount; i++) {
> 
>             nIndex = precipitateList[i];   /* index into the inprecipitateData precipitate list */
> 
>             dataBuf[bufIndex++] = (real8)inprecipitateData->precipitate[nIndex].myTag.domainID;
>             dataBuf[bufIndex++] = (real8)inprecipitateData->precipitate[nIndex].myTag.index;
> 
> 
> 
>             forcep = inprecipitateData->precipitate[nIndex].forcep;
>             dataBuf[bufIndex++] = forcep;
>             
>             r = inprecipitateData->precipitate[nIndex].r;
>             dataBuf[bufIndex++] = r;
> 
>             //for (armID = 0; armID < numNbrs; armID++) {
> 
> 
>              //Tag of the precipitate at the terminating end of the arm.
>  
>                 //dataBuf[bufIndex++] =
>                         //(real8)inprecipitateData->precipitate[nIndex].nbrTag[armID].domainID;
>                 //dataBuf[bufIndex++] =
>                         //(real8)inprecipitateData->precipitate[nIndex].nbrTag[armID].index;
>              
> 				
>              //bx, by, bz
>  
>                 //dataBuf[bufIndex++]=inprecipitateData->precipitate[nIndex].burgX[armID];
>                 //dataBuf[bufIndex++]=inprecipitateData->precipitate[nIndex].burgY[armID];
>                 //dataBuf[bufIndex++]=inprecipitateData->precipitate[nIndex].burgZ[armID];
> 
>             //nx, ny, nz
>  
>                 //dataBuf[bufIndex++]=inprecipitateData->precipitate[nIndex].nx[armID];
>                 //dataBuf[bufIndex++]=inprecipitateData->precipitate[nIndex].ny[armID];
>                 //dataBuf[bufIndex++]=inprecipitateData->precipitate[nIndex].nz[armID];
>             //}
> 
>             dataBuf[bufIndex++] = (real8)inprecipitateData->precipitate[nIndex].constraint;
> 
>             dataBuf[bufIndex++] = inprecipitateData->precipitate[nIndex].x; 
>             dataBuf[bufIndex++] = inprecipitateData->precipitate[nIndex].y; 
>             dataBuf[bufIndex++] = inprecipitateData->precipitate[nIndex].z; 
>         }
>    
>         if (bufIndex != maxIndex) {
>             Fatal("%s: buffer size %d did not match expected size %d",
>                   "PackInitialPrecipitateData", bufIndex, maxIndex);
>         }
> 
> /*
>  *      Return the buffer and its length to the caller
>  */
>         *buf     = dataBuf;
>         *bufEnts = bufIndex;
> 
>         return;
> }
> 
> 
> /*---------------------------------------------------------------------------
>  *
>  *      Function:     AddTagMapping
>  *      Description:  
>  *
>  *-------------------------------------------------------------------------*/
> //void AddTagMapping(Home_t *home, Tag_t *oldTag, Tag_t *newTag)
> //{
>         //int      newSize;
>         //TagMap_t *mapping;
> 
>         //if (home->tagMapEnts >= home->tagMapSize) {
>               //home->tagMapSize += NEW_NODEKEY_INC;
>               //newSize = home->tagMapSize * sizeof(TagMap_t);
>               //home->tagMap = (TagMap_t *)realloc(home->tagMap, newSize);
>         //}
> 
>         //mapping = &home->tagMap[home->tagMapEnts];
> 
>         //mapping->oldTag.domainID = oldTag->domainID;
>         //mapping->oldTag.index    = oldTag->index;
> 
>         //mapping->newTag.domainID = newTag->domainID;
>         //mapping->newTag.index    = newTag->index;
> 
>         //home->tagMapEnts++;
> 
>         //return;
> //}
> 
> 
> 
> 
> 
> /*---------------------------------------------------------------------------
>  *
>  *      Function:     GetNewTag
>  *      Description:  
>  *
>  *-------------------------------------------------------------------------*/
> //static void GetNewTag(Home_t *home, Tag_t *oldTag, Tag_t *newTag,
>                       //int *nextAvailableTag)
> //{
>         //int     nextTag, thisDomain;
> 
>         //nextTag    = *nextAvailableTag;
>         //thisDomain = home->myDomain;
> 
> /*
>  *      If the old tag belonged to a different domain, we must
>  *      give the precipitate a new tag.
>  */
>         //if (oldTag->domainID != thisDomain) {
> 
>             //for ( ; ; nextTag++) {
> 
> /*
>  *              Extend the precipitatekeys array if necessary.
>  */
>                 //if (nextTag >= home->newPrecipitateKeyMax) {
>                     //ExtendPrecipitateKeys(home, home->newPrecipitateKeyMax + NEW_NODEKEY_INC);
>                 //}
> 
>                 //if (home->precipitateKeys[nextTag] == (Precipitate_t *)NULL) {
>                     //newTag->domainID = thisDomain;
>                     //newTag->index = nextTag++;
>                     //AddTagMapping(home, oldTag, newTag);
>                     //break;
>                 //}
>             //}
>         //} else {
> /*
>  *          The old tag belonged to this domain...  check if that tag
>  *          is available for this run.
>  *
>  *          Extend the precipitatekeys array if necessary.
>  */
>             //if (oldTag->index >= home->newPrecipitateKeyMax) {
>                 //ExtendPrecipitateKeys(home, oldTag->index + NEW_NODEKEY_INC);
>             //}
> 
> /*
>  *          If the old tag is still available, use it and return
>  *          to the caller.
>  */
>             //if (home->precipitateKeys[oldTag->index] == (Precipitate_t *)NULL) {
>                 //newTag->domainID = oldTag->domainID;
>                 //newTag->index    = oldTag->index;
>                 //if (newTag->index >= home->newPrecipitateKeyPtr) {
>                     //home->newPrecipitateKeyPtr = newTag->index + 1;
>                 //}
>                 //return;
>             //}
> 
> /*
>  *          The old tag is no longer available, so just use
>  *          the next available tag.
>  */
>             //for ( ; ; nextTag++) {
> 
> /*
>  *              Extend precipitate keys array if necessary
>  */
>                 //if (nextTag >= home->newPrecipitateKeyMax) {
>                     //ExtendPrecipitateKeys(home, home->newPrecipitateKeyMax + NEW_NODEKEY_INC);
>                 //}
> 
>                 //if (home->precipitateKeys[nextTag] == (Precipitate_t *)NULL) {
>                     //newTag->domainID = thisDomain;
>                     //newTag->index = nextTag++;
>                     //AddTagMapping(home, oldTag, newTag);
>                     //break;
>                 //}
>             //}
>         //}
> 
>         //if (newTag->index >= home->newPrecipitateKeyPtr) {
>             //home->newPrecipitateKeyPtr = newTag->index + 1;
>         //}
> 
>         //*nextAvailableTag = nextTag;
> 
>         //return;
> //}
> 
> 
> /*---------------------------------------------------------------------------
>  *
>  *	Function:	UnpackInitialPrecipitateData
>  *	Description:	Pull the necessary precipitate data from the given
>  *			buffers sent from processor zero.  The nodal
>  *			data may be transmitted in multiple buffers, so
>  *			this function may be called multiple times
>  *			before all this domain's precipitates are received.
>  *			The final buffer will contain the full domain
>  *			decomposition for the problem space.
>  *
>  *	Arguments:
>  *
>  *-------------------------------------------------------------------------*/
> static void UnpackInitialPrecipitateData(Home_t *home, real8 *buf,
>                                   int *nextAvailableTag)
> {
>         int        i, armID, numNbrs, precipitatesInBuf, bufIndex;
>         Tag_t      oldTag, newTag;
>         Precipitate_t     *precipitate;
> 		real8		forcep;
> /*
>  *      Pull the precipitate count out of the buffer then loop through
>  *      all the nodal precipates provided.
>  */
>         bufIndex   = 0;
> 
>         precipitatesInBuf = (int)buf[bufIndex++];
> 
>         for (i = 0; i < precipitatesInBuf; i++) {
> 
>             precipitate = PopFreePrecipitateQ(home);
> 
>             oldTag.domainID = (int)buf[bufIndex++];
>             oldTag.index    = (int)buf[bufIndex++];
> 
> /*
>  *          Determine what tag will be used to identify this precipitate.
>  *          If the user requested that precipitate tags be preserved, the
>  *          new tag returned will be the same as the old tag
>  *          if at all possible.  If the new tag is in fact different
>  *          from the old tag, the GetNewTag() function will add
>  *          a mapping from the old to the new tag for later use
>  *          when reconciling all precipitates that have been retagged.
>  */
>             GetNewPrecipitateTag(home, &oldTag, &newTag, nextAvailableTag);
> 
>             precipitate->myTag.domainID = newTag.domainID;
>             precipitate->myTag.index    = newTag.index;
>             
> 			//printf("Precipitates domid %d tagindex %d \n",precipitate->myTag.domainID, precipitate->myTag.index );
> 			
>             precipitate->forcep = buf[bufIndex++];
>             precipitate->r = buf[bufIndex++];
>             
> 
>            // AllocPrecipitateArms(precipitate, numNbrs);
> 
>             //for (armID = 0; armID < numNbrs; armID++) {
> 
>                 //precipitate->nbrTag[armID].domainID = (int)buf[bufIndex++];
>                 //precipitate->nbrTag[armID].index    = (int)buf[bufIndex++];
> 
>                 //precipitate->burgX[armID] = buf[bufIndex++];
>                 //precipitate->burgY[armID] = buf[bufIndex++];
>                 //precipitate->burgZ[armID] = buf[bufIndex++];
> 
>                 //precipitate->nx[armID] = buf[bufIndex++];
>                 //precipitate->ny[armID] = buf[bufIndex++];
>                 //precipitate->nz[armID] = buf[bufIndex++];
>             //}
> 
>             precipitate->constraint = (int)buf[bufIndex++];
> 
>             precipitate->x = buf[bufIndex++]; 
>             precipitate->y = buf[bufIndex++]; 
>             precipitate->z = buf[bufIndex++]; 
> 
> /*
>  *          Register the precipitate in the precipitateKeys array and push the precipitate
>  *          onto the queue of precipitates native to this domain.
>  */
> 
> 		
> 			  	
>             home->precipitateKeys[precipitate->myTag.index] = precipitate;
>          
>             PushNativePrecipitateQ(home, precipitate);
> 			
>         }  /* for (i = 0; i < precipitatesInBuf; ... ) */
> 
>         return;
> }
> 
> 
> /*---------------------------------------------------------------------------
>  *
>  *      Function:     SendInitialPrecipitateData
>  *      Description:  Assigns all precipitates currently on the inprecipitateData precipitate list
>  *                    to remote domains based on the precipitate's coordinates and
>  *                    the current decomposition data, packages the data
>  *                    up and ships it off to the appropriate domains.
>  *                    Since the initial precipitate data may be read in blocks,
>  *                    this function may be called multiple times before
>  *                    all nodal data has been distributed to the remote
>  *                    domains.
>  *
>  *      Arguments:
>  *          inprecipitateData
>  *          msgCount
>  *          precipitateLists
>  *          listCounts
>  *
>  *-------------------------------------------------------------------------*/
> void SendInitialPrecipitateData(Home_t *home, InPrecipitateData_t *inprecipitateData, int *msgCount,
>                          int **precipitateLists, int *listCounts,
>                          int *nextAvailableTag)
> {
>         int         i, domID, numDomains, thisDomain;
>         int         sendIndex, selfSendIndex;
>         int         numSendBufs, numRecvBufs;
>         int         *sendBufEnts, *sendBufDest, *recvBufEnts;
>         real8       **sendBuf, **recvBuf;
> #ifdef PARALLEL
>         int         recvIndex;
> #endif
> 	
>         numDomains = home->numDomains;
>         thisDomain = home->myDomain;
> 
> /*
>  *      Normally, the following arrays are allocated in 
>  *      InitRemoteDomains(), but that has not yet been
>  *      called when we enter this, so for now, allocate
>  *      the arrays temporarily
>  */
> #ifdef PARALLEL
>         
>         home->inPRequests  = (MPI_Request *) malloc(home->numDomains *
>                                                   sizeof(MPI_Request));
>         home->outPRequests = (MPI_Request *) malloc(home->numDomains *
>                                                   sizeof(MPI_Request));
>         home->inPStatus    = (MPI_Status *) malloc(home->numDomains *
>                                                   sizeof(MPI_Status));
>         home->outPStatus   = (MPI_Status *) malloc(home->numDomains *
>                                                   sizeof(MPI_Status));                                          
>                                                   
>                                                   
>                                                   
>                                                   
>                                                   
>                                                   
>                                                   
>                                                   
> #endif
> 
> /*
>  *      Determine the number of remote domains to  which the
>  *      current domain must send data and allocate arrays
>  *      for the send buffer pointers and the send buffer
>  *      sizes.
>  */
>         numSendBufs   = 0;
>         selfSendIndex = -1;
> 
>         if (listCounts != (int *)NULL) {
>             for (domID = 0; domID < numDomains; domID++) {
>                 if (listCounts[domID] != 0) {
>                     if (domID == thisDomain) {
>                         selfSendIndex = numSendBufs;
>                     }
>                     numSendBufs++;
>                 }
>             }
> 
>             if (numSendBufs > 0) {
>                 sendBuf     = (real8 **)calloc(1, numSendBufs * sizeof(real8 *));
>                 sendBufEnts = (int *)calloc(1, numSendBufs * sizeof(int));
>                 sendBufDest = (int *)calloc(1, numSendBufs * sizeof(int));
>             }
> 
> /*
>  *          Pack up precipitate data for each remote domain to which this
>  *          domain will be sending data
>  */
>             sendIndex = 0;
> 
>             for (domID = 0; domID < numDomains; domID++) {
>                 if (listCounts[domID] != 0) {
>     
>                     PackInitialPrecipitateData(inprecipitateData, precipitateLists[domID],
>                                         listCounts[domID], &sendBuf[sendIndex],
>                                         &sendBufEnts[sendIndex]);
> 
>                     sendBufDest[sendIndex] = domID;
>                     sendIndex++;
> 
>                     if (sendIndex >= numSendBufs) {
>                         break;
>                     }
>                 }
>             }
>         }  /* if (listCount != NULL) */
> 
> /*
>  *      Allocate array in which to receive incoming buffer lengths.
>  *      Do not include any buffer this domain has for itself; any
>  *      such buffer will be handled separately.
>  */
>         numRecvBufs = msgCount[thisDomain] - (selfSendIndex >= 0);
> 
>         if (numRecvBufs > 0) {
>             recvBuf = (real8 **)calloc(1, numRecvBufs * sizeof(real8 *));
>             recvBufEnts = (int *)calloc(1, numRecvBufs * sizeof(int));
>         }
>       
> /*
>  *      Pre-issue receives of buffer lengths from any domain that will
>  *      be sending data to the current domain.
>  */
> #ifdef PARALLEL
>         for (i = 0; i < numRecvBufs; i++) {
>             MPI_Irecv(&recvBufEnts[i], 1, MPI_INT, MPI_ANY_SOURCE,
>                       MSG_INIT_PLENS, MPI_COMM_WORLD, &home->inPRequests[i]);
>         }
> 
> /*
>  *      Have the current domain send out the sizes of the buffers it
>  *      will be transmitting.  Note: If the buffer is for itself, the
>  *      data is not actually transmitted
>  */
>         for (i = 0; i < numSendBufs; i++) {
>             if (i == selfSendIndex) {
>                 home->outPRequests[i] = MPI_REQUEST_NULL;
>             } else {
>                 MPI_Isend(&sendBufEnts[i], 1, MPI_INT, sendBufDest[i],
>                           MSG_INIT_PLENS, MPI_COMM_WORLD,
>                           &home->outPRequests[i]);
>             }
>         }
> #endif
> 
> /*
>  *      If this domain packed up a buffer for itself, unpack that buffer
>  *      now.
>  */
>         if (selfSendIndex >= 0) {
>             UnpackInitialPrecipitateData(home, sendBuf[selfSendIndex],
>                                   nextAvailableTag);
>         }
> 
> /*
>  *      Wait for all length send/receives to complete
>  */
>  
>  
> #ifdef PARALLEL
>         MPI_Waitall(numSendBufs, home->outPRequests, home->outPStatus);
>         MPI_Waitall(numRecvBufs, home->inPRequests, home->inPStatus);
> 
> /*
>  *      Allocate the receive buffers and post the receives
>  *      associated with those buffers.
>  */
>         for (i = 0; i < numRecvBufs; i++) {
>             recvBuf[i] = (real8 *)malloc(recvBufEnts[i] * sizeof(real8));
>             MPI_Irecv(recvBuf[i], recvBufEnts[i], MPI_DOUBLE,
>                       home->inPStatus[i].MPI_SOURCE, MSG_INIT_PRECIPITATES,
>                       MPI_COMM_WORLD, &home->inPRequests[i]);
>         }
> 
> /*
>  *      Send out all the packed buffers to the appropriate
>  *      remote domains.  -- Unless the destination is the current domain
>  *      in which case the buffer has already been processed.
>  */
>         for (i = 0; i < numSendBufs; i++) {
>             
>             if (i == selfSendIndex) {
>                 home->outPRequests[i] = MPI_REQUEST_NULL;
>             } else {
>                 MPI_Isend(sendBuf[i], sendBufEnts[i], MPI_DOUBLE,
>                           sendBufDest[i], MSG_INIT_PRECIPITATES, MPI_COMM_WORLD,
>                           &home->outPRequests[i]);
>             }
>         }
> 
> /*
>  *      Process the incoming buffers as soon as they arrive.
>  */
>         for (i = 0; i < numRecvBufs; i++) {
>             MPI_Waitany(numRecvBufs, home->inPRequests, &recvIndex,
>                         home->inPStatus);
>             UnpackInitialPrecipitateData(home, recvBuf[recvIndex], nextAvailableTag);
>             free(recvBuf[recvIndex]);
>         }
> 
> /*
>  *      Wait for all buffer sends to complete.
>  */
>         MPI_Waitall(numSendBufs, home->outPRequests, home->outPStatus);
> 
>         free(home->inPRequests);
>         free(home->outPRequests);
>         free(home->inPStatus);
>         free(home->outPStatus);
> 
>         home->inPRequests  = (MPI_Request *)NULL;
>         home->outPRequests = (MPI_Request *)NULL;
>         home->inPStatus    = (MPI_Status *)NULL;
>         home->outPStatus   = (MPI_Status *)NULL;
> #endif
>    
>         for (i = 0; i < numSendBufs; i++) {
>             free(sendBuf[i]);
>         }
> 
>         if (numSendBufs > 0) {
>             free(sendBuf);
>             free(sendBufEnts);
>             free(sendBufDest);
>         }
> 
>         if (numRecvBufs > 0) {
>             free(recvBuf);
>             free(recvBufEnts);
>         }
>         
> 	return;
> }
> 
> 
> 
> 
> 
> 
> 
> 
