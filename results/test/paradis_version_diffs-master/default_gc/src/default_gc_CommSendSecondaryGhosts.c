103c103
< 
---
> 		//printf("add secondary request %d remdom %d nodetag (%d,%d) nbrtag (%d,%d) \n",home->myDomain,remDom->domainIdx,node->myTag.domainID,node->myTag.index,node->nbrTag[segID].domainID,node->nbrTag[segID].index);
160c160
<         neededSpace = FLTS_PER_GHOST2_NODE + FLTS_PER_GHOST_ARM * node->numNbrs;
---
>         neededSpace = FLTS_PER_GHOST2_NODE + FLTS_PER_GHOST_ARM * node->numNbrs+FLTS_PER_GHOST_PNRB* node->numPNbrs;
163c163,164
<             allocSizeIncr = 5 * (FLTS_PER_GHOST2_NODE + (2*FLTS_PER_GHOST_ARM));
---
> 			//allocSizeIncr = 5 * (FLTS_PER_GHOST2_NODE + (2*FLTS_PER_GHOST_ARM)+(3*FLTS_PER_GHOST_PNRB));
>             allocSizeIncr = 5 * (FLTS_PER_GHOST2_NODE + (2*FLTS_PER_GHOST_ARM)+(50*FLTS_PER_GHOST_PNRB));
170a172
>         buf[bufOffset++] = (real8)node->numPNbrs;
175a178,181
>         
>         buf[bufOffset++] = node->nodedx;
>         buf[bufOffset++] = node->nodedy;
>         buf[bufOffset++] = node->nodedz;
216a223,232
>         
>         for (i = 0; i < node->numPNbrs; i++) {
>             buf[bufOffset++] = (real8)node->PnbrTag[i].domainID;
>             buf[bufOffset++] = (real8)node->PnbrTag[i].index;
> 
>         }
>         
>         
>         
>         
265c281,282
<         allocatedSize = 5 * (FLTS_PER_GHOST2_NODE + (2*FLTS_PER_GHOST_ARM)) + 1;
---
> 		//allocatedSize = 5 * (FLTS_PER_GHOST2_NODE + (2*FLTS_PER_GHOST_ARM)+(3*FLTS_PER_GHOST_PNRB)) + 1;
>         allocatedSize = 5 * (FLTS_PER_GHOST2_NODE + (2*FLTS_PER_GHOST_ARM)+(50*FLTS_PER_GHOST_PNRB)) + 1;
275c292,295
< 
---
>         
> 		//if(home->myDomain == 17){
> 		//PrintAllprecipitateNBR(home);}
> 		
277c297
< 
---
> 			
280c300
< 
---
> 			
282c302
< 
---
> 			//PrintAllprecipitateNBR(home);
283a304,307
> 				
> 				PrintAllprecipitateNBR(home);
> 				printf("Fatal mydomain %d tag (%d,%d) \n",home->myDomain,tag.domainID,tag.index);
> 				PrintNodelist(home);
320,321c344,345
<         int            i, j, remDomID, armID, bufOffset;
<         int            numNodes, numArms, nodeVals, totRemDomCount;
---
>         int            i, j, remDomID, armID,iPNbr, bufOffset;
>         int            numNodes, numArms,numPNbrs , nodeVals, totRemDomCount;
338c362,363
<  */
---
>  */			
> 			
341c366
< 
---
> 			//printf("unpack secondary ghost mydomain %d tag (%d,%d) \n",home->myDomain,tmpTag.domainID,tmpTag.index);
348a374
>  //tassa?				
349a376
>                 numPNbrs = (int)inBuf[bufOffset+3];
351c378
<                            FLTS_PER_GHOST_ARM * numArms;
---
>                            FLTS_PER_GHOST_ARM * numArms+FLTS_PER_GHOST_PNRB*numPNbrs;
352a380,382
>                 //printf("duplicatenode tag(%d %d) numArms %d nodeVals %d \n",tmpNode->myTag.domainID,tmpNode->myTag.index,numArms,nodeVals);
>                 
>                 
373a404
>             numPNbrs = (int)inBuf[bufOffset++];
377a409,412
>             
>             node->nodedx = inBuf[bufOffset++];
>             node->nodedy = inBuf[bufOffset++];
>             node->nodedz = inBuf[bufOffset++];
404c439,440
< 
---
> 			AllocNodePrecipitates(node,numPNbrs );
> 			
409c445,447
< 
---
> 				if(node->nbrTag[armID].domainID < 0 || node->nbrTag[armID].index < 0){
> 				 printf("unpacksecondary myDomain %d tag ( %d %d ) nbrtag ( %d %d) \n",home->myDomain,node->myTag.domainID,node->myTag.index,node->nbrTag[armID].domainID,node->nbrTag[armID].index)	
> 					;}
421a460,470
>             
>              /* Precipitate neighbors (AL)*/
>                      for (iPNbr = 0; iPNbr < node->numPNbrs; iPNbr++) {
>                         
>                         node->PnbrTag[iPNbr].domainID=inBuf[bufOffset++];
>                         node->PnbrTag[iPNbr].index=inBuf[bufOffset++];
>  
>                         
>                     }
>             
>             
551a601
>  //Tarkasta tämä, jossain getneighbornodessa on vika
700a751
>         
