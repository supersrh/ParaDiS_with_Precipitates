153a154,180
> /*---------------------------------------------------------------------------
>  *
>  *      Function:     AllocNodePrecipitateArrays
>  *      Description:  This function allocates all arrays in the node
>  *                    structure that are dependent on the number of 
>  *                    segments attached to the node.
>  *
>  *                    WARNING!  This must be either kept in sync with the
>  *                    functions in Util.c for managing node arms, or
>  *                    merged with that function.
>  *
>  *      Arguments:
>  *          node      pointer to the node for which arrays are to be
>  *                    allocated.
>  *          armCount  number of segments attached to the node
>  *
>  *-------------------------------------------------------------------------*/
> static void AllocNodePrecipitateArrays(Node_t *node, int numPNbrs)
> {
>         node->PnbrTag   = (Tag_t *)malloc(numPNbrs * sizeof(Tag_t));
>         
> 
>         
> 
>         return;
> }
> 
192a220,230
> void FreeNodePnbrArrays(Node_t *node)
> {
>             free(node->PnbrTag);
>             
> 
>         return;
> }
> 
> 
> 
> 
213c251
<         int    i, armCount;
---
>         int    i, armCount,prCount;
215a254
>         prCount	= source->numPNbrs;	
240a280,285
>         
>         for (i = 0; i < prCount; i++) {
> 			dest->PnbrTag[i].domainID = source->PnbrTag[i].domainID;
>             dest->PnbrTag[i].index    = source->PnbrTag[i].index;
> 			
> 		}
272a318
>  //printf("Backupnode %d \n",home->myDomain);
274a321
>         AllocNodePrecipitateArrays(bkupNode, origNode->numPNbrs);
309a357,358
>  
>   //printf("Restorenode %d \n",home->myDomain);
310a360
>         FreeNodePnbrArrays(origNode);
312a363
>         AllocNodePrecipitateArrays(origNode, origNode->numPNbrs);
1838a1890
>                     FreeNodePnbrArrays(&bkupNodeList[k]);
1918c1970
< 	int	i, j, k, tarm, *tarmList, *tarmList2;
---
> 	int	i, j, k, tarm, *tarmList, *tarmList2,nop;
2047c2099
< 
---
> 	FreeNodePNbrs(newNode);
2050c2102,2110
< 
---
> 	
> 	nop=node->numPNbrs;
> 	
> 	//this assigns  numPNbrs to newNode (AL)
> 	AllocNodePrecipitates(newNode,nop);
> 	
> 	for (i =0; i < nop; i++) {
> 		newNode->PnbrTag[i]=node->PnbrTag[i];
> 	}
2087a2148,2152
>         
>         (*splitNode1)->nodedx =node->nodedx+(pos1[X]-node->x);
>         (*splitNode1)->nodedy =node->nodedy+(pos1[Y]-node->y);
>         (*splitNode1)->nodedz =node->nodedz+(pos1[Z]-node->z);
>         
2095a2161,2165
>         
>         (*splitNode2)->nodedx =node->nodedx+(pos2[X]-node->x);
>         (*splitNode2)->nodedy =node->nodedy+(pos2[Y]-node->y);
>         (*splitNode2)->nodedz =node->nodedz+(pos2[Z]-node->z);
>         
2158,2159c2228,2232
<                 newNode->vZ);           /* commandeer the AddOp normal*/
<                                         /* arguments for velocity data*/
---
>                 newNode->vZ,			/* plane normal is not needed */
>                 newNode->nodedx,         /* but node velocity is, so we*/   
>                 newNode->nodedy,         /* commandeer the AddOp normal*/   
>                 newNode->nodedz);         /* arguments for velocity data*/  
>                                         
2167c2240,2241
<                 0.0, 0.0, 0.0);         /* plane normal not needed */
---
>                 0.0, 0.0, 0.0,			 /* plane normal not needed */
>                  node->nodedx, node->nodedy, node->nodedz);        
2520a2595,2600
>             
>             targetNode->nodedx = targetNode->nodedx+(position[X]-targetNode->x);
>             targetNode->nodedy = targetNode->nodedy+(position[Y]-targetNode->y);
>             targetNode->nodedz = targetNode->nodedz+(position[Z]-targetNode->z);
>             
>             
2533c2613,2614
<                     0.0, 0.0, 0.0);         /* plane normal not needed */
---
>                     0.0, 0.0, 0.0,			/* plane normal not needed */
>                     targetNode->nodedx, targetNode->nodedy, targetNode->nodedz);         
