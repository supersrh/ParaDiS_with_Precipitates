103c103,105
< 
---
>         //if(node->myTag.domainID != remDom->domainIdx ){
> 		//printf("add secondary request %d remdom %d nodetag (%d,%d) nbrtag (%d,%d) \n",home->myDomain,remDom->domainIdx,node->myTag.domainID,node->myTag.index,node->nbrTag[segID].domainID,node->nbrTag[segID].index);
> 	     //}
160c162
<         neededSpace = FLTS_PER_GHOST2_NODE + FLTS_PER_GHOST_ARM * node->numNbrs;
---
>         neededSpace = FLTS_PER_GHOST2_NODE + FLTS_PER_GHOST_ARM * node->numNbrs+FLTS_PER_GHOST_PNRB* node->numPNbrs;
163c165,166
<             allocSizeIncr = 5 * (FLTS_PER_GHOST2_NODE + (2*FLTS_PER_GHOST_ARM));
---
> 			//allocSizeIncr = 5 * (FLTS_PER_GHOST2_NODE + (2*FLTS_PER_GHOST_ARM)+(3*FLTS_PER_GHOST_PNRB));
>             allocSizeIncr = 5 * (FLTS_PER_GHOST2_NODE + (2*FLTS_PER_GHOST_ARM)+(50*FLTS_PER_GHOST_PNRB));
170a174
>         buf[bufOffset++] = (real8)node->numPNbrs;
175a180,183
>         
>         buf[bufOffset++] = node->nodedx;
>         buf[bufOffset++] = node->nodedy;
>         buf[bufOffset++] = node->nodedz;
216a225,234
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
265c283,284
<         allocatedSize = 5 * (FLTS_PER_GHOST2_NODE + (2*FLTS_PER_GHOST_ARM)) + 1;
---
> 		//allocatedSize = 5 * (FLTS_PER_GHOST2_NODE + (2*FLTS_PER_GHOST_ARM)+(3*FLTS_PER_GHOST_PNRB)) + 1;
>         allocatedSize = 5 * (FLTS_PER_GHOST2_NODE + (2*FLTS_PER_GHOST_ARM)+(50*FLTS_PER_GHOST_PNRB)) + 1;
275c294,296
< 
---
>         
> 		
> 		
277c298
< 
---
> 			
280c301
< 
---
> 			
282c303,305
< 
---
> 			//PrintAllprecipitateNBR(home);
> 			
> 			//printf("mydomain %d asked domain and node tag (%d,%d) \n",home->myDomain,tag.domainID,tag.index);
283a307,311
> 				
> 				//PrintAllprecipitateNBR(home);
> 				printf("Fatal mydomain %d tag (%d,%d) \n",home->myDomain,tag.domainID,tag.index);
> 				
> 				//PrintNodelist(home);
320,321c348,349
<         int            i, j, remDomID, armID, bufOffset;
<         int            numNodes, numArms, nodeVals, totRemDomCount;
---
>         int            i, j, remDomID, armID,iPNbr, bufOffset;
>         int            numNodes, numArms,numPNbrs , nodeVals, totRemDomCount;
338c366,367
<  */
---
>  */			
> 			
341c370
< 
---
> 			//printf("unpack secondary ghost mydomain %d tag (%d,%d) \n",home->myDomain,tmpTag.domainID,tmpTag.index);
348a378
>  //tassa?				
349a380
>                 numPNbrs = (int)inBuf[bufOffset+3];
351c382
<                            FLTS_PER_GHOST_ARM * numArms;
---
>                            FLTS_PER_GHOST_ARM * numArms+FLTS_PER_GHOST_PNRB*numPNbrs;
352a384,386
>                 //printf("duplicatenode tag(%d %d) numArms %d nodeVals %d \n",tmpNode->myTag.domainID,tmpNode->myTag.index,numArms,nodeVals);
>                 
>                 
373a408
>             numPNbrs = (int)inBuf[bufOffset++];
377a413,416
>             
>             node->nodedx = inBuf[bufOffset++];
>             node->nodedy = inBuf[bufOffset++];
>             node->nodedz = inBuf[bufOffset++];
404c443,444
< 
---
> 			AllocNodePrecipitates(node,numPNbrs );
> 			
409c449,451
< 
---
> 				if(node->nbrTag[armID].domainID < 0 || node->nbrTag[armID].index < 0){
> 				 printf("unpacksecondary myDomain %d tag ( %d %d ) nbrtag ( %d %d) \n",home->myDomain,node->myTag.domainID,node->myTag.index,node->nbrTag[armID].domainID,node->nbrTag[armID].index)	
> 					;}
421a464,474
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
551a605
>  //Tarkasta tämä, jossain getneighbornodessa on vika
700a755
>         
