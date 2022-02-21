28a29
> #include "InPrecipitateData.h"
735a737,739
>                 
>                 
>                 
768,769c772,774
< 
<                     AllocNodeArms(node, numNbrs); 
---
> 					
>                     AllocNodeArms(node, numNbrs);
>                     AllocNodePrecipitates(node,0); 
786a792,794
> //printf("readnodaldata  %e %d \n", node->x, node->nbrTag[iNbr].index );
> 
> 
849a858,859
> 
> 
916a927,929
>  
>  
>  
941a955,1143
>  
>         return;
> }
> 
> 
> /*---------------------------------------------------------------------------
>  *
>  *      Function:       AssignPrecipitatesToDomains
>  *      Description:    Using the provided domain decomposition data,
>  *                      loop though all precipitates in the inprecipitateData precipitate array
>  *                      and assign each precipitate to the domain which encompasses
>  *                      the precipitate's coordinates.
>  *
>  *      Arguments:
>  *          precipitateCount   number of precipitates in the inprecipitateData precipitate array.
>  *          precipitateLists   Location in which to return to the caller the
>  *                      array of lists of precipitates to be sent to the
>  *                      remote domains.  The precipitateLists array has
>  *                      a list pointer for every domain in the problem.
>  *                      (each list is an array of indices into the
>  *                      inprecipitateData precipitate array)
>  *          listCounts  Location in which to return to the caller an
>  *                      array of integers indicating the count of precipitates
>  *                      on the corresponding list in <precipitateLists>.
>  *
>  *-------------------------------------------------------------------------*/
> void AssignPrecipitatesToDomains(Home_t *home, InPrecipitateData_t *inprecipitateData, int precipitateCount,
>                           int ***precipitateLists, int **listCounts)
> {
>         int      i, j;
>         int      nXdoms, nYdoms, nZdoms;
>         int      domIndex, nDoms, len, nexti;
>         int      *qHead, *qCount, *list;
>         int      **listArray, *countArray;
>         Param_t  *param;
>         Precipitate_t   *precipitate;
> 
> 
>         if (precipitateCount == 0) {
>             *precipitateLists = (int **)NULL;
>             *listCounts = (int *)NULL;
>             return;
>         }
> 
>         param = home->param;
> 
>         nXdoms = param->nXdoms;
>         nYdoms = param->nYdoms;
>         nZdoms = param->nZdoms;
> 
>         nDoms = nXdoms * nYdoms * nZdoms;
> 
>         listArray  = (int **) calloc(1, nDoms * sizeof(int *));
>         countArray = (int *) calloc(1, nDoms * sizeof(int));
> 
> /*
>  *      Allocate and initialize an array of integers (1 per domain)
>  *      to be used as pointers to the head of the queue of precipitates
>  *      assigned to the domains.
>  */
>         qHead  = (int *)malloc(nDoms * sizeof(int));
>         qCount = (int *)malloc(nDoms * sizeof(int));
> 
>         for (i = 0; i < nDoms; i++) {
>             qHead[i] = -1;
>             qCount[i] = 0;
>         }
> 
> /*
>  *      Loop through all the precipitates on the current precipitate list, find
>  *      the proper domain for the precipitate based on the precipitate coordinates
>  *      and add the precipitate to the domain's queue.
>  */
>         for (i = 0; i < precipitateCount; i++) {
> 
>             precipitate = &inprecipitateData->precipitate[i];
>             domIndex = FindCoordDomain(home, 0, &precipitate->x, &precipitate->y, &precipitate->z);
>             //printf(" %d %d %d \n",i,domIndex,precipitate->myTag.index);
> 
> /*
>  *          Each precipitate on a queue contains the index of the next precipitate
>  *          on the queue.  To add a precipitate to the queue, we just set
>  *          that precipitate's 'next precipitate' pointer to the current head of
>  *          the queue, and then point the queue head to the new precipitate.
>  */
>             nexti = qHead[domIndex];
>             if (nexti < 0) {
>                 precipitate->next = (Precipitate_t *)NULL;
>             } else {
>                 precipitate->next = &inprecipitateData->precipitate[qHead[domIndex]];
>             }
>             qHead[domIndex] = i;
>             qCount[domIndex]++;
>         }
>         
>          
> 
> /*
>  *      For each domain, generate the list of indices in the precipitate array
>  *      for precipitates to be sent to the domain.
>  */
>         for (domIndex = 0; domIndex < nDoms; domIndex++) {
> 
>             list  = (int *)NULL;
>             len   = qCount[domIndex];
>             nexti = qHead[domIndex];
> 
>             if (len > 0) list = (int *)malloc(len * sizeof(int));
> 
>             for (j = 0; j < len; j++) {
>                 list[j] = nexti;
> /*
>  *              Do some pointer arithmetic to get convert pointer
>  *              addresses into array indices.
>  */
>                 if (inprecipitateData->precipitate[nexti].next == (Precipitate_t *)NULL) {
>                     nexti = -1;
>                 } else {
>                     nexti = inprecipitateData->precipitate[nexti].next - inprecipitateData->precipitate;
>                 }
>                 
>               
>             }
> 
>             if (nexti != -1) {
>                 Fatal("Queue error. domain %d, queue len %d", domIndex, len);
>             }
> 
>             listArray[domIndex] = list;
>             countArray[domIndex] = len;
>         }
> 
> /*
>  *      Free up the now unneeded arrays.  (the lists of precipitate indices
>  *      will be freed elsewhere when they are no longer needed)
>  */
>         free(qHead);
>         free(qCount);
> 
>         *precipitateLists  = listArray;
>         *listCounts = countArray;
> 
>         return;
> }
> 
> 
> 
> 
> 
> 
> /*---------------------------------------------------------------------------
>  *
>  *      Function:       FreePrecipitateLists
>  *      Description:    Releases temporary storage associated with 
>  *                      the arrays of precipitate indices identifying precipitates
>  *                      to be sent to remote domains.
>  *
>  *      Arguments:
>  *          precipitateLists   Address of the the array of lists of precipitates to
>  *                      that were sent to the remote domains.  The precipitateLists
>  *                      array has a list pointer for every domain in the
>  *                      problem.  On return, the contents of this adress
>  *                      will be zeroed.
>  *          listCounts  Address of the array of integers indicating the
>  *                      count of precipitates on the corresponding list in
>  *                      <precipitateLists>.  On return, the contents of this
>  *                      adress will be zeroed.
>  *
>  *-------------------------------------------------------------------------*/
> void FreePrecipitateLists(Home_t *home, int ***precipitateLists, int **listCounts)
> {
>         int dom;
> 
>         if (*listCounts != (int *)NULL) {
>             free(*listCounts);
>             *listCounts = (int *)NULL;
>         }
> 
>         if (*precipitateLists != (int **)NULL) {
>             for (dom = 0; dom < home->numDomains; dom++) {
>                 if ((*precipitateLists)[dom] != (int *)NULL) {
>                     free((*precipitateLists)[dom]);
>                     (*precipitateLists)[dom] = (int *)NULL;
>                 }
>             }
>             free(*precipitateLists);
>             *precipitateLists = (int **)NULL;
>         }
> 
943a1146,1725
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
> void ReadPrecipitateDataFile(Home_t *home, InPrecipitateData_t *inprecipitateData,InData_t *inData, char *dataFile)
> {
>         int      i, dom, iNbr;
>         int      maxTokenLen, tokenType, pIndex, okay;
>         int      valType, numVals;
>         int      numDomains, numReadTasks;
>         int      nextFileSeg, maxFileSeg, segsPerReader, taskIsReader;
>         int      distIncomplete, readCount, numNbrs;
>         int      fileSegCount = 1, fileSeqNum = 0;
>         int      *globalMsgCnt, *localMsgCnt, **precipitateLists, *listCounts;
>         int      miscData[2];
>         int      localPrecipitateCount, globalPrecipitateCount;
>         int      binRead = 0;
>         int      nextAvailableTag = 0;
>         real8    burgSumX, burgSumY, burgSumZ;
>         void     *valList;
>         char     inLine[500], token[256];
>         char     baseFileName[256], tmpFileName[256];
>         FILE     *fpSeg;
>         Precipitate_t   *precipitate;
>         Param_t  *param;
>         ParamList_t *precipitateParamList;
> 
>         param      = home->param;
>         numDomains = home->numDomains;
>         fpSeg      = (FILE *)NULL;
> 
>         precipitateParamList = home-> precipitateParamList;
> 
>         globalMsgCnt = (int *)malloc((numDomains+1) * sizeof(int));
>         localMsgCnt  = (int *)malloc((numDomains+1) * sizeof(int));
> 
>         memset(inLine, 0, sizeof(inLine));
>         maxTokenLen = sizeof(token);
> 
> /*
>  *      Only domain zero reads the initial stuff...
>  */
>         if (home->myDomain == 0) {
> 
>             if (dataFile == (char *)NULL) {
>                 Fatal("ReadPrecipitateDataFile: No data file provided");
>             }
> 
>             snprintf(tmpFileName, sizeof(tmpFileName), "%s", dataFile);
> 
>             if (!strcmp(&tmpFileName[strlen(tmpFileName)-2], ".0")) {
>                     tmpFileName[strlen(tmpFileName)-2] = 0;
>             }
> 
>             snprintf(baseFileName, sizeof(baseFileName), "%s", tmpFileName);
> 
> 	
> 
> 
> /*
>  *          Try to open the nodal data file.  If the specified file
>  *          can not be opened, try looking for a segmented data file
>  *          (i.e. look for the first file segment <dataFile>.0). If
>  *          that can't be opened either, then exit with an error.
>  */
> 
>             if ((fpSeg = fopen(tmpFileName, "r")) == (FILE *)NULL) {
>                 snprintf(tmpFileName, sizeof(tmpFileName), "%s.%d",
>                          baseFileName, fileSeqNum);
>                 if ((fpSeg = fopen(tmpFileName, "r")) == (FILE *)NULL) {
>                     Fatal("Error %d opening file %s to read nodal data",
>                           errno, dataFile);
>                 }
>             }
>             
>             
> 
> /*
>  *          Get the first token.  This should either be a known
>  *          parameter identifier, or a file version number.
>  */
>             tokenType = GetNextToken(fpSeg, token, maxTokenLen);
>             pIndex = LookupParam(precipitateParamList, token);
> 
>             if (pIndex < 0) {
> /*
>  *              If the token does not correspond to a known parameter, it
>  *              should be the version number, so convert it to an integer
>  *              and parse the rest of the file the old way (i.e. values
>  *              are positional, not identified by name/value pairs)
>  */
>                 param->dataFileVersion = atoi(token);
> 
>                 if ((param->dataFileVersion < 1) ||
>                     (param->dataFileVersion > 3)) {
>                     Fatal("ReadPrecipitateDatFile: Unsupported file version %d",
>                           param->dataFileVersion);
>                 }
> 
>                 ReadPreV4DataParams(home, fpSeg, &inprecipitateData->decomp);
> 
>             } else {
> /*
>  *              Just go through the nodal data file reading all
>  *              the associated parameters.  Need to do special
>  *              processing of the domain decomposition, and when
>  *              when we hit the 'nodaldata' identifier, just break
>  *              out of this loop.
>  */
>  
>  
> 
>                 while ((tokenType != TOKEN_ERR) && (tokenType != TOKEN_NULL)) {
> 
>                     if (pIndex >= 0) {
> /*
>  *                      Token represents a known parameter identifier, so
>  *                      read the associated value(s).
>  */
>                         valType = precipitateParamList->varList[pIndex].valType;
>                         numVals = precipitateParamList->varList[pIndex].valCnt;
>                         valList = precipitateParamList->varList[pIndex].valList;
>                         okay = GetParamVals(fpSeg, valType, numVals, valList);
>                         if (!okay) {
>                             Fatal("Parsing Error obtaining values for "
>                                   "parameter %s\n",
>                                   precipitateParamList->varList[pIndex].varName);
>                         }
>     
>                     } else {
> /*
>  *                      Token does not represent one of the simple 
>  *                      parameters.  If it's not one of the identifiers
>  *                      that needs special handling, skip it.
>  */
>                         if (strcmp(token, "domainDecomposition") == 0) {
> /*
>  *                          The minSide{XYZ} and maxSide{XYZ} values are
>  *                          now redundant but until the rest of the code
>  *                          is modified to remove references to these
>  *                          values, we need to explicitly set them now.
>  */
>                             param->minSideX = param->minCoordinates[X];
>                             param->minSideY = param->minCoordinates[Y];
>                             param->minSideZ = param->minCoordinates[Z];
> 
>                             param->maxSideX = param->maxCoordinates[X];
>                             param->maxSideY = param->maxCoordinates[Y];
>                             param->maxSideZ = param->maxCoordinates[Z];
> 
>                             tokenType = GetNextToken(fpSeg, token, maxTokenLen);
> /*
>  *                          Do a quick verification of the decomposition type
>  */
>                             if ((param->dataDecompType < 1) ||
>                                 (param->dataDecompType > 2)) {
>                                 Fatal("dataDecompType=%d is invalid.  Type must be 1 or 2\n",
>                                       param->dataDecompType);
>                             }
> 							 
>                             ReadDecompBounds(home, (void **)&fpSeg, binRead,
>                                              param->decompType,
>                                              &inprecipitateData->decomp);
> 
>                         } else if (strcmp(token, "precipitateData") == 0) {
> 							
> 						
> /*
>  *                          When we hit the nodal data, we can break
>  *                          out of the loop because we are assuming
>  *                          all other data file parameters have been
>  *                          processed.  If they have not, we have a
>  *                          problem since processing of the nodal data
>  *                          requires the other parameters.
>  *                          Note: Remainder of the file should just
>  *                          contain " = " followed by the nodal data,
>  *                          so be sure to skip the next token before
>  *                          processing the nodal data.
>  */
>                             tokenType = GetNextToken(fpSeg, token, maxTokenLen);
>                             break;
>                         } else {
> /*
>  *                          If the parameter is not recognized, skip the
>  *                          parameter and any associated value(s).
>  */
>                             printf("Ignoring unknown data file parameter %s\n",
>                                    token);
>                             valType = V_NULL;
>                             numVals = 0;
>                             valList = (void *)NULL;
>                             okay = GetParamVals(fpSeg, valType, numVals,
>                                                 valList);
>                         }
>                     }
> 
>                     tokenType = GetNextToken(fpSeg, token, maxTokenLen);
> 
>                     if ((tokenType == TOKEN_NULL)||(tokenType == TOKEN_ERR)) {
>                         Fatal("Parsing error on file %s\n", tmpFileName);
>                     }
>                     pIndex = LookupParam(precipitateParamList, token);
>                 }
>             }
> 
> /*
>  *          Need to set some of the values that are dependent on
>  *          the simulation size before we go any further.
>  */
>             param->Lx = param->maxSideX - param->minSideX;
>             param->Ly = param->maxSideY - param->minSideY;
>             param->Lz = param->maxSideZ - param->minSideZ;
> 
>             param->invLx = 1.0 / param->Lx;
>             param->invLy = 1.0 / param->Ly;
>             param->invLz = 1.0 / param->Lz;
>             
>      
> 
> /*
>  *          If we did not get a domain decomposition from the restart
>  *          file (whether because of a mismatch in domain geometry or
>  *          domain decomposition type between the current run and
>  *          that from the restart file) we'll need to create a new
>  *          uniform decomposition to start off.
>  * 
>  */
>  
>  
>  
>             if (inprecipitateData->decomp == (void *)NULL) {
>                 //printf("Generating uniform domain decomposition.\n");
>                UniformDecomp(home, &inprecipitateData->decomp);
>                 //inprecipitateData->decomp=inData->decomp;
>             }
> 
>         }  /* if (home->myDomain == 0) */
> 
>         if (numDomains >= param->numFileSegments) {
>             numReadTasks = param->numFileSegments;
>         } else {
>             numReadTasks = MIN(((param->numFileSegments + (numDomains - 1)) /
>                                numDomains), numDomains);
>         }
>         
>         
>             
>         
> #ifdef PARALLEL
> /*
>  *      Domain zero now needs to pass various pieces of data to
>  *      the remote domains before any processes begin to read
>  *      and distribute the nodal data... start with the param structure.
>  */
>         MPI_Bcast((char *)param, sizeof(Param_t), MPI_CHAR, 0, MPI_COMM_WORLD);
> 
> /*
>  *      Now the base nodal data file name.
>  */
>         MPI_Bcast((char *)baseFileName, sizeof(baseFileName),
>                   MPI_CHAR, 0, MPI_COMM_WORLD);
> 
> /*
>  *      Also send off the count of domains that will be involved in
>  *      the parallel read of nodal data and the number of files segments
>  *      comprising the nodal data.
>  */
>         miscData[0] = numReadTasks;
>         miscData[1] = param->numFileSegments;
> 
>         MPI_Bcast(miscData, 2, MPI_INT, 0, MPI_COMM_WORLD);
> 
>         numReadTasks = miscData[0];
>         fileSegCount = miscData[1];
> 
> #endif
> 
> 
>  
> 
> /*
>  *      Lastly, only domain zero knows the current domain decomposition.
>  *      Invoke a function to distribute that data to remote domains (if
>  *      any).
>  */
>         BroadcastDecomp(home, inprecipitateData->decomp);
>         
>       
> 
> /*
>  *      Have each domain determine which (if any) of the
>  *      nodal data file segments it is to read in.
>  */
>         segsPerReader = (fileSegCount + (numReadTasks - 1)) / numReadTasks;
>         nextFileSeg   = home->myDomain * segsPerReader;
>         maxFileSeg    = MIN((nextFileSeg+(segsPerReader-1)),(fileSegCount-1));
>         taskIsReader  = (nextFileSeg <= maxFileSeg);
>    
> /*
>  *      All processes loop until all nodal data has been read in
>  *      and distributed to the appropriate domains.
>  */
>         distIncomplete = 1;
>         readCount = 0;
> 
>         while (distIncomplete) {
> 
> /*
>  *          Have each reader task allocate a buffer for the next
>  *          block of precipitates to be read in.  Then have the readers
>  *          read their next blocks of precipitates.
>  */
>  
>  //Aiheutti jotain paskaa(AL)
>   //inprecipitateData->precipitate = (Precipitate_t *)calloc(1, MAX_PRECIPITATES_PER_BLOCK *
>                                                 //sizeof(Precipitate_t));
>  
>  
>  
>   
>             if (taskIsReader) {
> 				
>                 inprecipitateData->precipitate = (Precipitate_t *)calloc(1, MAX_PRECIPITATES_PER_BLOCK *
>                                                 sizeof(Precipitate_t));
>             }
>  
>             while (taskIsReader) {
> 
> /*
>  *              If the current domain doesn't have a data file segment
>  *              opened up yet, open the next segment now.
>  */
>                 if ((fpSeg == (FILE *)NULL) && (nextFileSeg <= maxFileSeg)) {
>                
>                     sprintf(tmpFileName, "%s.%d", baseFileName, nextFileSeg);
>                     fpSeg = fopen(tmpFileName, "r");
> 
>                     if (fpSeg == (FILE *)NULL) {
>                         Fatal("Task %d: Error %d opening %s", home->myDomain,
>                               errno, tmpFileName);
>                     }
>                     
>                 }
>                 
>             
>   
> /*
>  *              If we have a valid open file pointer now, continue
>  *              reading nodal data.  If we don't, it's because the
>  *              task has no more data to read so terminate its
>  *              'reader' status.
>  */
>                 if (fpSeg != (FILE *)NULL) {
> 
>                     Getline(inLine, sizeof(inLine), fpSeg);
> 
> /*
>  *                  If we hit EOF, close the current segment and
>  *                  loop back to open the next segment (if any)
>  *                  this domain is responsible for reading.
>  */
>                     if (inLine[0] == 0) {
>                         nextFileSeg++;
>                         fclose(fpSeg);
>                         fpSeg = (FILE *)NULL;
>                         continue;  
>                     }
> 
> /*
>  *                  Got information on the next precipitate, so deal with it.
>  */
>                     precipitate = &inprecipitateData->precipitate[readCount++];
> 					
>                     sscanf(inLine, "%d,%d %lf %lf %lf %lf %lf %d",
>                            &precipitate->myTag.domainID, &precipitate->myTag.index,
>                            &precipitate->x, &precipitate->y, &precipitate->z,
>                            &precipitate->forcep,&precipitate->r, &precipitate->constraint);
> 					
>                     /*AllocPrecipitateArms(precipitate, numNbrs);*/ 
> /*
>  *                  Read in data for each arm of the precipitate
>  */
>                     burgSumX = 0.0;
>                     burgSumY = 0.0;
>                     burgSumZ = 0.0;
> 				
>                     //for (iNbr = 0; iNbr < precipitate->numNbrs; iNbr++) {
> 
>                         //Getline(inLine, sizeof(inLine), fpSeg); 
>                         //sscanf(inLine, "%d,%d %lf %lf %lf",
>                                //&precipitate->nbrTag[iNbr].domainID,
>                                //&precipitate->nbrTag[iNbr].index,
>                                //&precipitate->burgX[iNbr],
>                                //&precipitate->burgY[iNbr],
>                                //&precipitate->burgZ[iNbr]);
> 
>                         //Getline(inLine, sizeof(inLine), fpSeg);
>                         //sscanf(inLine, "%lf %lf %lf", &precipitate->nx[iNbr],
>                                //&precipitate->ny[iNbr], &precipitate->nz[iNbr]);
> 
>                         //Normalize(&precipitate->nx[iNbr], &precipitate->ny[iNbr], &precipitate->nz[iNbr]);
> 
>                         //burgSumX += precipitate->burgX[iNbr];
>                         //burgSumY += precipitate->burgY[iNbr];
>                         //burgSumZ += precipitate->burgZ[iNbr];
>                     //}
>  
>                     FoldBox(param, &precipitate->x, &precipitate->y, &precipitate->z);
>                     
>                   
> 
> /*
>  *                  Just a quick sanity check, to make sure burgers
>  *                  vector is conserved for all unconstrained precipitates.
>  */
>                     if (precipitate->constraint == UNCONSTRAINED) {
> 
>                         if ((fabs(burgSumX) > 0.0001) ||
>                             (fabs(burgSumY) > 0.0001) ||
>                             (fabs(burgSumZ) > 0.0001)) {
> 
>                             printf("Error: precipitate (%d,%d)\n",
>                                    precipitate->myTag.domainID, precipitate->myTag.index);
> 
>                             //for (iNbr=0; iNbr < precipitate->numNbrs; iNbr++) {
>                                 //printf("  arm[%d] burg = %e %e %e\n",
>                                        //iNbr, precipitate->burgX[iNbr],
>                                        //precipitate->burgY[iNbr],
>                                        //precipitate->burgZ[iNbr]);
>                             //}
> 
>                             //Fatal("Burger's vector not conserved!");
>                         }
>                     }
>                     
>                     
>                     
>                 } else {
> /*
>  *                  No more files to be read by this task
>  */
>                     taskIsReader = 0;
>                 }
> 
> /*
>  *              We read the nodal data in blocks rather than all at once.
>  *              If this domain has read in the maximum number of precipitates
>  *              allowed in a block, stop reading data and get ready
>  *              to distribute it to the appropriate remote domains.
>  */
>    
>  
>                 if (readCount >= MAX_PRECIPITATES_PER_BLOCK) {
>                     break;
>                 }
> 
>             }  /* while (taskIsReader) */
> 
> /*
>  *          Determine the domains to which to send any precipitates
>  *          read in by this domain.
>  */
> 
>   
> 		    
>                   
>             AssignPrecipitatesToDomains(home, inprecipitateData, readCount, &precipitateLists,
>                                  &listCounts);
>                                
> 
> /*
>  *          Set up an array (1 element per domain).  Each reader
>  *          task sets to 1 the array entry for each remote domain
>  *          to which it will be sending data.  When we do a global
>  *          reduction to sum up the arrays, the resulting array
>  *          contains the number of messages that will be sent to
>  *          each domain during this communication.
>  *
>  *          Plus 1 extra element set to 1 if ANY process is sending
>  */
>             memset(globalMsgCnt, 0, (numDomains+1) * sizeof(int));
> 
> #ifdef PARALLEL
> 
>             memset(localMsgCnt, 0, (numDomains+1) * sizeof(int));
> 
>             if (listCounts != (int *)NULL) {
>                 for (dom = 0; dom < numDomains; dom++) {
>                     if (listCounts[dom] > 0) {
>                        localMsgCnt[dom] = 1;
>                        localMsgCnt[numDomains] = 1;
>                     }
>                 }
>             }
> 
>             MPI_Allreduce(localMsgCnt, globalMsgCnt, numDomains + 1,
>                           MPI_INT, MPI_SUM, MPI_COMM_WORLD);
> 
> #else
>             if (listCounts != (int *)NULL) {
>                 for (dom = 0; dom < numDomains; dom++) {
>                     if (listCounts[dom] > 0) {
>                        globalMsgCnt[dom] = 1;
>                        globalMsgCnt[numDomains] = 1;
>                     }
>                 }
>             }
> #endif
> 
> /*
>  *          If no domains have any more data to send out, we're
>  *          done reading/distributing the data.
>  */
>             if (globalMsgCnt[numDomains] == 0) {
>                 distIncomplete = 0;
>                 FreeInPrecipitateArray(inprecipitateData, readCount);
>                 FreePrecipitateLists(home, &precipitateLists, &listCounts);
>                 continue;
>             }
>             
>               
> 
> /*
>  *          Do next send/receive of nodal data
>  */
>             SendInitialPrecipitateData(home, inprecipitateData, globalMsgCnt, precipitateLists,
>                                 listCounts, &nextAvailableTag); 
>   
>             FreeInPrecipitateArray(inprecipitateData, readCount);
>             FreePrecipitateLists(home, &precipitateLists, &listCounts);
>             readCount = 0;
> 
>         }  /* if (distIncomplete) */
> 
> /*
>  *      This is a good place for a quick sanity check that the sum of
>  *      precipitates on all domains equals the total precipitate count from the
>  *      data file.
>  */
>         localPrecipitateCount = 0;
>         globalPrecipitateCount = 0;
>         
>          
> 
>         for (i = 0; i < home->newPrecipitateKeyPtr; i++) {
> 			
>             if (home->precipitateKeys[i] != (Precipitate_t *)NULL) {
>                 localPrecipitateCount++;
>                 
>             }
>         }
> 
> #ifdef PARALLEL
>         MPI_Reduce(&localPrecipitateCount, &globalPrecipitateCount, 1, MPI_INT, MPI_SUM,
>                    0, MPI_COMM_WORLD);
> 
> #else
>         globalPrecipitateCount = localPrecipitateCount;
> #endif
>  
>         if ((home->myDomain == 0) && (param->precipitateCount != globalPrecipitateCount)) {
>             Fatal("ReadPrecipitatedataFile: Read %d precipitates, expected %d!",
>                   globalPrecipitateCount, param->precipitateCount);
>         }
> 
>         free(localMsgCnt);
>         free(globalMsgCnt);
> 
>         return;
> }
> 
