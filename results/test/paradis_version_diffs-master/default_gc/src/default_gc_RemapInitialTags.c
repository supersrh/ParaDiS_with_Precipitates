100a101,160
> /*-------------------------------------------------------------------------
>  *
>  *      Function:     PackPrecipitateTagMap
>  *      Description:  Packs a buffer with old/new tag mappings for
>  *                    all local nodes which were retagged during
>  *                    the initial node distribution.
>  *      Arguments:
>  *          buf       Pointer to location in which to return to the
>  *                    caller a pointer to the buffer created by
>  *                    this function.
>  *          bufSize   Pointer to location in which to return to the
>  *                    caller the size (in bytes) of the buffer
>  *                    returned in <buf>.
>  *
>  *------------------------------------------------------------------------*/
> static void PackPrecipitateTagMap(Home_t *home, int **buf, int *bufSize)
> {
>         int      i, mapBufEnts, mapBufSize, bufIndex;
>         int      *mapBuf;
>         TagMap_t *mapping;
> 
>         mapBufEnts = 1 + (home->precipitatetagMapEnts * 4);
>         mapBufSize = mapBufEnts * sizeof(int);
>         mapBuf = (int *)malloc(mapBufSize);
> 
>         bufIndex = 0;
> 
>         mapBuf[bufIndex++] = home->precipitatetagMapEnts;
> 
>         for (i = 0; i < home->precipitatetagMapEnts; i++) {
>  
>             mapping = &home->precipitatetagMap[i];
> 
>             mapBuf[bufIndex++] = mapping->oldTag.domainID;
>             mapBuf[bufIndex++] = mapping->oldTag.index;
>             mapBuf[bufIndex++] = mapping->newTag.domainID;
>             mapBuf[bufIndex++] = mapping->newTag.index;
>         }
> 
>         *buf     = mapBuf;
>         *bufSize = mapBufSize;
> 
>         return;
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
> 
135a196,236
> /*-------------------------------------------------------------------------
>  *
>  *      Function:     UnpackPrecipitateTagMap
>  *      Description:  Unpacks the provided tag map from a remote
>  *                    domain and adds all old/new tag mappings
>  *                    to the local tag map.
>  *
>  *      Arguments:
>  *          buf       Pointer to buffer to be unpacked.  Buffer is 
>  *                    assumed to be an array of integers.
>  *
>  *------------------------------------------------------------------------*/
> static void UnpackPrecipitateTagMap(Home_t *home, int *buf)
> {
>         int     i, numMappings, bufIndex;
>         Tag_t   oldTag, newTag;
> 
>         bufIndex = 0;
> 
>         numMappings = buf[bufIndex++];
> 
>         for (i = 0; i < numMappings; i++) {
> 
>             oldTag.domainID = buf[bufIndex++];
>             oldTag.index    = buf[bufIndex++];
>             newTag.domainID = buf[bufIndex++];
>             newTag.index    = buf[bufIndex++];
> 			//printf("UnpackPrecipitate tag map (%d %d)->(%d %d)\n",oldTag.domainID,oldTag.index,newTag.domainID,newTag.index );
>             AddPrecipitateTagMapping(home, &oldTag, &newTag);
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
180a282
> 					
182a285
>                   
189a293,357
> /*-------------------------------------------------------------------------
>  *
>  *      Function:     RemapPNbrTags
>  *      Description:  Loop through all precipitate neighbors of the local nodes retagging
>  *                    the precipitate nbr if the precipitate nbr
>  *                    was retagged during initialization.  The full
>  *                    tag map (home->tagMap) must have been populated
>  *                    with new mappings for the local and all nearby
>  *                    domains, and sorted before this function is called.
>  *
>  *                    Note:  We don't alter the local node tags here
>  *                    since they were retagged (if necessary) during
>  *                    the initial node distribution.
>  *
>  *------------------------------------------------------------------------*/
> static void RemapPNbrTags(Home_t *home)
> {
>         int      i, pnbrID, maxNodeKey;
>         Node_t   *node;
>         TagMap_t key;
>         TagMap_t *mapping;
> 
>         maxNodeKey = home->newNodeKeyPtr;
> 		
>         for (i = 0; i < maxNodeKey; i++) {
> 
>             if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
>                 continue;
>             }
> 
>             for (pnbrID = 0; pnbrID < node->numPNbrs; pnbrID++) {
> 
>                 key.oldTag.domainID = node->PnbrTag[pnbrID].domainID;
>                 key.oldTag.index    = node->PnbrTag[pnbrID].index;
> 
> /*
>  *              Use the precipitateneighbors's tag to look for any new tag mapping.
>  *              If no mapping was found, the tag is fine as it is.
>  *              Otherwise, reset the arm's tag to the new value.
>  */
>                 mapping = bsearch(&key, home->precipitatetagMap, home->precipitatetagMapEnts,
>                                   sizeof(TagMap_t), TagMapCompare);
> 
>                 if (mapping != (TagMap_t *)NULL) {
> 					//if(home->myDomain==11){
> 					//printf("remap pnbr(%d,%d)->(%d,%d) node (%d,%d)\n",key.oldTag.domainID, key.oldTag.index ,mapping->newTag.domainID,mapping->newTag.index,node->myTag.domainID,node->myTag.index);
> 					//}
>                     node->PnbrTag[pnbrID].domainID = mapping->newTag.domainID;
>                     node->PnbrTag[pnbrID].index    = mapping->newTag.index;
>                   
>                 }
>             }
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
209a378,381
> 
>  //PackTagMap(home, &tagBuf, &tagBufLen);
>  
> 
213a386
>          
216c389,391
<         
---
>         //if(home->myDomain==11){   
>     //printf("Tassa1 %d %d inbufLen %d sourcerank %d  \n",home->cycle,i,remDom->inBufLen,remDomIndex);
> 	//}       
219a395,397
>          //if(home->myDomain==11){   
>     //printf("Tassa2 %d %d inbufLen %d sourcerank %d  \n",home->cycle,i,remDom->inBufLen,remDomIndex);
> 	//}   
221a400,401
> 
> 
227a408,411
> 
> 
> 
> 
234a419
>              
245a431,432
>  
>   
247c434
< 
---
> 		
283a471
>             //printf("domain %d instatus %d \n",home->myDomain,reqStatus.MPI_TAG);            
305a494,496
>   
>  //printf("before remap %d %d prtagmapents %d tagmapents %d \n",home->myDomain,home->cycle,home->precipitatetagMapEnts,home->tagMapEnts);
>  //MPI_Barrier(MPI_COMM_WORLD);
310a502
>             //RemapPNbrTags(home);
311a504,519
>         
>       
>         
>          //if(home->myDomain==11){   
>     //printf("before prremap %d %d   \n",home->myDomain,home->cycle);
> 	//}   
>      /*remap's precipitate nbrlist (AL)*/
>          if (home->precipitatetagMapEnts > 0 || home->tagMapEnts>0 ) {
>             qsort(home->precipitatetagMap, home->precipitatetagMapEnts, sizeof(TagMap_t),
>                   TagMapCompare);
> 
>             RemapPNbrTags(home);
>         }
>         
>         
>          
322a531,542
> 		free(home->precipitatetagMap);
> 
>         home->precipitatetagMap = (TagMap_t *)NULL;
>         home->precipitatetagMapSize = 0;
>         home->precipitatetagMapEnts = 0;
> 
> //printf("after remap %d %d   \n",home->myDomain,home->cycle);
>  MPI_Barrier(MPI_COMM_WORLD);
>  
> 
> 
> 
324a545,705
> 
> 
> 
> /*-------------------------------------------------------------------------
>  *
>  *      Function:     DistributePrecipitatesTagMaps
>  *      Description:  Have each domain send to all its near-neighbors
>  *                    a mapping between old and new tags for all nodes
>  *                    the domain retagged during initialization.  And
>  *                    of course, receive tag mappings from the remote
>  *                    domains.
>  *
>  *------------------------------------------------------------------------*/
> void DistributePrecipitateTagMaps(Home_t *home) 
> {
> #ifdef PARALLEL
>         int             i, remDomIndex, reqIndex = 0, tagBufLen;
>         int             *tagBuf, *remBuf;
>         RemoteDomain_t  *remDom;
>         MPI_Status      reqStatus;
> 
> /*
>  *	Pre-issue receives of tag map buffer lengths from each neighbor
>  */
>         for (i = 0; i < home->remoteDomainCount; i++) {
>             remDomIndex = home->remoteDomains[i];
>             remDom = home->remoteDomainKeys[remDomIndex];
>         
>             MPI_Irecv(&remDom->inPBufLen, 1, MPI_INT, remDomIndex,
>                       MSG_PR_TAGREMAP_LEN, MPI_COMM_WORLD,
>                       &home->inPRequests[i]);
>         }
> 
> /*
>  *	Pack a buffer with mappings between old and new tags for any
>  *      local nodes this domain has retagged.
>  */
>         PackPrecipitateTagMap(home, &tagBuf, &tagBufLen);
> 
> /*
>  *	Send the length length (in bytes) of the tag map buffer to
>  *      the neighboring domains who will be receiving the buffer.
>  */
>         for (i = 0; i < home->remoteDomainCount; i++) {
>             remDomIndex = home->remoteDomains[i];
>             remDom = home->remoteDomainKeys[remDomIndex];
>             MPI_Isend(&tagBufLen, 1, MPI_INT, remDomIndex,
>                       MSG_PR_TAGREMAP_LEN, MPI_COMM_WORLD,
>                       &home->outPRequests[i]);
>         }
> 
> /*
>  *	Wait until we've received the tag map buffer lengths from all
>  *	neighboring domains, allocate buffers for the maps and pre-
>  *	issue receives for the data.
>  */
>         MPI_Waitall(home->remoteDomainCount, home->inPRequests, home->inPStatus);
>         
>        
> 
>         for (i = 0; i < home->remoteDomainCount; i++) {
> 
>             remDomIndex = home->remoteDomains[i];
>             remDom = home->remoteDomainKeys[remDomIndex];
> 
>             remDom->inPBuf = (char *) malloc(remDom->inPBufLen);
> 
>             MPI_Irecv(remDom->inPBuf, remDom->inPBufLen / sizeof(int), MPI_INT,
>                       remDomIndex, MSG_PR_TAGREMAP, MPI_COMM_WORLD,
>                       &home->inPRequests[i]);
>         }
> 
> /*
>  *	Wait for all length sends to complete
>  */
>         MPI_Waitall(home->remoteDomainCount, home->outPRequests,home->outPStatus);
> 		 
> /*
>  *	Send the local tag-map to all neighboring domains
>  */
>         for (i = 0; i < home->remoteDomainCount; i++) {
> 
>             remDomIndex = home->remoteDomains[i];
>             remDom = home->remoteDomainKeys[remDomIndex];
> 		
>             MPI_Isend(tagBuf, tagBufLen / sizeof(int), MPI_INT, remDomIndex,
>                       MSG_PR_TAGREMAP, MPI_COMM_WORLD, &home->outPRequests[i]);
>         } 
> 
> /*
>  *	Wait for tag-map from any remote domain as long as there are any
>  *	outstanding receives.
>  */
>         for (i = 0; i < home->remoteDomainCount; i++) {
>             MPI_Waitany(home->remoteDomainCount, home->inPRequests,
>                         &reqIndex, &reqStatus);
>             remDomIndex = home->remoteDomains[reqIndex];
>             remBuf = (int *)home->remoteDomainKeys[remDomIndex]->inPBuf;
>             UnpackPrecipitateTagMap(home, remBuf);
>             free(remBuf);
>             home->remoteDomainKeys[remDomIndex]->inPBuf = (char *)NULL;
>         }
> 
> /*
>  *	Now just wait for all sends to complete (should be long done by
>  *	now) and free up any temporary buffers.
>  */
>         MPI_Waitall(home->remoteDomainCount, home->outPRequests,
>                     home->outPStatus);
> 
>         free(tagBuf);
> #endif
> 
> /*
>  *      All new tag mappings are known, so sort the tag mappings
>  *      and retag all arms that need it (if there are any in this
>  *      domain).
>  */
>         //if (home->precipitatetagMapEnts > 0) {
>             //qsort(home->precipitatetagMap, home->precipitatetagMapEnts, sizeof(TagMap_t),
>                   //TagMapCompare);
> 
>            
>         //}
> 
> /*
>  *      No longer need the tagMap stuff in the home structure, so free
>  *      the array and reinitialize associated values. Commented this part away because we free the precipitate tagmap stuff in DistributeTagMaps() function. (AL)
>  */
>         //free(home->precipitatetagMap);
> 
>         //home->precipitatetagMap = (TagMap_t *)NULL;
>         //home->precipitatetagMapSize = 0;
>         //home->precipitatetagMapEnts = 0;
>         
>         
>         
>         
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
> 
> 
> 
> 
> 
> 
> 
> 
> 
