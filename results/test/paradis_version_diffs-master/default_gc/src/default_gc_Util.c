113a114
> #include "InPrecipitateData.h"
755c756
< 
---
> 			
1190a1192,1206
> static void InitNodePNbr(Node_t *node, int pnbr)
> {
>         node->PnbrTag[pnbr].domainID = -1;
>         node->PnbrTag[pnbr].index = -1;
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
1251c1267
<             InitNodeArm(node, i);
---
>              InitNodeArm(node, i);
1255a1272,1320
> /*-------------------------------------------------------------------------
>  *
>  *      Function:    AllocNodePrecipitates
>  *      Description: Allocate memory for 'n' precipitate neighbors in the specified 
>  *                   node structure.
>  *
>  *      Arguments:
>  *          node       Pointer to node for which to allocate precipitates
>  *
>  *------------------------------------------------------------------------*/
> void AllocNodePrecipitates(Node_t *node, int n)
> {
>         int i;
> 
> /*
>  *      The node structure *may* have some arm arrays already allocated.
>  *      If this is the case but the number of currently allocated
>  *      arms is not what we need, we have to free the previously
>  *      allocated arrays before allocating new ones.
>  */
> 	// if (node->numPNbrs != n) {
> 
> 	  // if (node->numPNbrs > 0) {
> 
>                 free(node->PnbrTag);
>                 
> 		// }
>              
>            //if (node->numPNbrs == 0) {
>              //n=n+1
>              
> 		    //}
>                 
>             node->numPNbrs = n;
>         
>             node->PnbrTag = (Tag_t *)malloc(n*sizeof(Tag_t));
>             
> 	    // }
>         
>  for (i = 0; i < node->numPNbrs; i++) {
>             InitNodePNbr(node, i);
>         }
> 
>         return;
> }
> 
> 
> 
> 
1281a1347
> 	FreeNodePNbrs(home->nodeKeys[index]);
1282a1349,1351
>     
> 	return;
> }
1283a1353,1361
> void FreePrecipitate(Home_t *home, int index)
> {
> 	
> 	RemovePrecipitateFromCellQ(home, home->precipitateKeys[index]);
> 	RemovePrecipitateFromCell2Q(home, home->precipitateKeys[index]);
> 	RecyclePrecipitateTag(home, index);
> 	PushFreePrecipitateQ (home, home->precipitateKeys[index]);
> 	home->precipitateKeys[index] = (Precipitate_t *)NULL;
> 	
1287a1366
> 
1317c1396
< 
---
> 		
1319c1398
< 
---
> 	
1323a1403,1413
> void FreeNodePNbrs(Node_t *node)
> {
> 	 //if (node->numPNbrs == 0) {
>            // return;
>         //}
>       free(node->PnbrTag);      node->PnbrTag = (Tag_t *)NULL;  
>       node->numPNbrs = 0; 
> 	  return;
> }
> 
> 
1354c1444,1448
< 
---
>         
>          if (node->nbrTag[n].domainID < 0 || node->nbrTag[n].index < 0) {
>         PrintNode(node);
> 		printf(" getneighbornode  mydomain %d tag (%d,%d) native %d cycle %d \n",home->myDomain,node->nbrTag[n].domainID,node->nbrTag[n].index,node->native,home->cycle);
> 		}
1370a1465,1470
> 				
> 				 if (node->nbrTag[n].domainID < 0 || node->nbrTag[n].index < 0) {
>          PrintNode(node);
> 		printf("mydomain %d node (%d,%d) nbr tag (%d,%d) native %d  n %d\n",home->myDomain,node->myTag.domainID,node->myTag.index,node->nbrTag[n].domainID,node->nbrTag[n].index,node->native,n);
> 		}
> 				
1384a1485,1543
> Precipitate_t *GetNeighborPrecipitate(Home_t *home, Node_t *node, int n)
> {
>         int    i, j;
>         Precipitate_t *neighbor;
> #if 0
> /*
>  *      New version which assumes no invalid arms in arm list
>  */
>         if (n >= node->numNbrs) {
>             printf("GetNeighborNode: Error finding neighbor %d\n", n);
>             PrintNode(node);
>             return((Precipitate_t *)NULL);
>         }
> 
>         neighbor = GetPrecipitateFromTag(home, node->PnbrTag[n]);
> 
>         return(neighbor);
> #else
> /*
>  *      Old version which assumes the arm list may be sparsely
>  *      populated and returns the n'th valid neighbor, which may
>  *      not be at index n.  
>  */
>         j = -1;
>  
>         for (i = 0; i < node->numNbrs; i++) {
> 
>             if (node->nbrTag[i].domainID >= 0) j++;
> 
>             if (j == n) {
>                 neighbor = GetPrecipitateFromTag(home, node->PnbrTag[i]);
>                 return(neighbor);
>             }
>         }
> 
>         printf("GetNeighborNode returning NULL for node (%d,%d) nbr %d\n",
>                node->myTag.domainID, node->myTag.index, n);
>         PrintNode(node);
> 
>         return((Precipitate_t *)NULL);
> #endif
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
1468a1628
> 		int				i;
1472a1633,1635
> 			
> 			printf("getnode %d cycle %d \n",home->myDomain,home->cycle);
> 			
1474a1638,1647
> 			
> 			for (i=0; i<home->newNodeKeyPtr;i++) {
> 
> 		if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
> 			continue;
> 			}
>         PrintNode(node);
> 			}
> 			
> 			
1475a1649,1650
>         
>          
1483a1659
> 				//printf(" tag.index >= home->newNodeKeyPtr %d %d \n",tag.domainID,tag.index);
1487a1664
> 				//printf("node = home->nodeKeys[tag.index]) == (Node_t *)NULL %d %d \n",tag.domainID,tag.index);
1497c1674
<  *          either the remote doamin or the remote node.  Hence, it
---
>  *          either the remote domain or the remote node.  Hence, it
1502a1680
> 				//printf("remDom == NULL %d %d \n",tag.domainID,tag.index);
1506a1685
> 				//printf("tag.index >= remDom->maxTagIndex %d %d \n",tag.domainID,tag.index);
1515a1695,1909
> 
> /*-------------------------------------------------------------------------
>  *
>  *      Function:     GetPrecipitateFromIndex
>  *      Description:  Given the domain ID and index return
>  *                    a pointer to the corresponding precipitate.
>  *
>  *      Arguments:
>  *          domID   ID of the domain owning the precipitate in question
>  *          index   Index in the <precipitateKeys> array for domain <domID>
>  *                  of the node in question.
>  *
>  *      Returns:    Pointer to the requested precipitate if found.
>  *                  NULL in all other cases.
>  *
>  *------------------------------------------------------------------------*/
> Precipitate_t *GetPrecipitateFromIndex(Home_t *home, int domID, int index)
> {
>         Precipitate_t         *precipitate;
>         RemoteDomain_t *remDom;
> 
>    
>         if ((domID < 0) || (index < 0)) {
>             return((Precipitate_t *)NULL);
>         }
> 
> /*
>  *      If the node in question is owned by the current domain,
>  *      look up the node in the local <nodeKeys> array.
>  */
>         if (domID == home->myDomain) {
> 
>             if (index >= home->newPrecipitateKeyPtr) {
>                 return((Precipitate_t *)NULL);
>             }
> 
>             if ((precipitate = home->precipitateKeys[index]) == (Precipitate_t *)NULL) {
>                 return((Precipitate_t *)NULL);
>             }
> 
>             return(precipitate);
> 
>         } else {
> /*
>  *          Node is owned by a remote domain, so look up the 
>  *          node in the appropriate <remoteDomainKeys> array.
>  */
>             remDom = home->remoteDomainKeys[domID];
> 
>             if (remDom == (RemoteDomain_t *)NULL) {
>                 return((Precipitate_t *)NULL);
>             }
> 
>             if (index >= remDom->maxprecipitateTagIndex) {
>                 return((Precipitate_t *)NULL);
>             }
> 
>             if ((precipitate = remDom->precipitateKeys[index]) == (Precipitate_t *)NULL) {
>                 return((Precipitate_t *)NULL);
>             }
> 
>             return(precipitate);
>         }
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
> 
> 
> 
> 
> 
> /*-------------------------------------------------------------------------
>  *
>  *      Function:     GetPrecipitateFromTag
>  *      Description:  Given a precipitate tag, returns a pointer to the
>  *                    corresponding precipitate structure.
>  *
>  *                    NOTE: If the specified tag is for a local precipitate
>  *                    and the local precipitate is not found, the function
>  *                    will force a code abort with a fatal error.
>  *
>  *      Arguments:
>  *          tag    Tag identifying the desired precipitate
>  *
>  *      Returns:    Pointer to the requested precipitate if found.
>  *                  NULL in all other cases.
>  *
>  *------------------------------------------------------------------------*/
> Precipitate_t *GetPrecipitateFromTag (Home_t *home, Tag_t tag)
> {
>         Precipitate_t         *precipitate;
>         RemoteDomain_t *remDom;
>    
>         if (tag.domainID < 0 || tag.index < 0) {
> 			printf("getprcipitate fatal (%d %d)",tag.domainID,tag.index);
>             Fatal("GetPrecipitateFromTag: invalid tag (%d,%d)",
>                   tag.domainID, tag.index);
>         }
> 
> /*
>  *      If the tag is for a local domain, look up the precipitate in
>  *      the local <precipitateKeys> array.
>  */
>         if (tag.domainID == home->myDomain) {
> 
>             if (tag.index >= home->newPrecipitateKeyPtr) {
> 			printf("Precipitate null1 (%d %d)",tag.domainID,tag.index);
>                 return((Precipitate_t *)NULL);
>             }
> 
>             if ((precipitate = home->precipitateKeys[tag.index]) == (Precipitate_t *)NULL) {
> 				printf("Precipitate null2 (%d %d)",tag.domainID,tag.index);
>                 return((Precipitate_t *)NULL);
>                 
>             }
> 
>             return(precipitate);
> 
>         } else {
> /*
>  *          If the precipitate is in a remote domain, there are valid situations
>  *          in which the current domain does not have information on
>  *          either the remote domain or the remote precipitate.  Hence, it
>  *          is not an error to return a NULL pointer.
>  */
>             remDom = home->remoteDomainKeys[tag.domainID];
> 
>             if (remDom == NULL) {
> 				printf("Precipitate null3 (%d %d)",tag.domainID,tag.index);
>                 return((Precipitate_t *)NULL);
>             }
> 
>             if (tag.index >= remDom->maxprecipitateTagIndex) {
> 				printf("Precipitate null4 (%d %d)",tag.domainID,tag.index);
>                 return((Precipitate_t *)NULL);
>             }
>      // printf("getprecipitate before remote (%d %d)",tag.domainID,tag.index);
>             precipitate = remDom->precipitateKeys[tag.index];
>              //printf("getprecipitate after remote (%d %d)",tag.domainID,tag.index);     
>             return(precipitate);
>         }
> }
> 
> 
> 
> 
> 
> 
> 
> void PrintAllprecipitateNBR(Home_t *home)
> {
> 	Node_t *node;
> 	int i,j;
> 	
> 	printf("Mydomain %d\n",home->myDomain);
> 	
> 	for (i=0; i<home->newNodeKeyPtr;i++) {
> 
> 		if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
> 			continue;
>         }
> 	
> 		PrintNode(node);
> 		//printf("Node (%d,%d)  ",node->myTag.domainID,node->myTag.index);
> 		//for (j = 0; j < node->numPNbrs; j++) {
> 		//	printf(" Pnbr (%d,%d) \n", node->PnbrTag[j].domainID,node->PnbrTag[j].index);
>          
> 		//	}
> 	}
> 	return;
> }
> 
> 
> void PrintNodelist(Home_t *home)
> {	Node_t *node;
> 	int i,j;
> 
> 
> for (i=0; i<home->newNodeKeyPtr;i++) {
> 	
> 	
> 
> 		if ((node = home->nodeKeys[i]) == (Node_t *)NULL ) {
> 			continue;
>         }
> 	
> 		
> 		printf("Nativenode %d (%d,%d) native %d \n ",home->myDomain,node->myTag.domainID,node->myTag.index,node->native);
> 		
> 	}
> 
> 	
> 
> 
> }
> 
> 
> 
> 
1537a1932,1937
>         
>         for (i = 0; i < node->numPNbrs; i++) {
>             printf("(Pnbr %d,%d) ", node->PnbrTag[i].domainID,
>                    node->PnbrTag[i].index);
>         }
>         
1541a1942,1943
> 
> 
1542a1945
>         printf("\n");
1602a2006,2087
> void PrintPrecipitate(Precipitate_t *precipitate)
> {
>  
> /*
>  *      Print the precipitate Data
>  */
>         printf("  %d,%d %.15e %.15e %.15e %.15e %.15e  %d  \n",
>                precipitate->myTag.domainID, precipitate->myTag.index, 
>                precipitate->x, precipitate->y, precipitate->z,
>                precipitate->forcep,precipitate->r,precipitate->constraint);
>         return;
> }
> 
> 
> //this is for simple proximity check for updating neighborlist
> int NodeinPrecipitateforNbrList(real8 R,real8 xp,real8 yp,real8 zp, real8 x1,real8 y1,real8 z1)
> {
> 		real8   dx, dy, dz;
> 		
> /* we solve x^2+y^2+z^2=R^2, x=x1+tdx..,dx=x2-x1 */		
> 		dx=xp-x1;
> 		dy=yp-y1;
> 		dz=zp-z1;
>         
> 		
> 		
> 		
> 		if(sqrt(dx*dx+dy*dy+dz*dz)>R){ 
> 			//printf(" Node not inprecipitate %e (%e %e %e) \n ",sqrt(dx*dx+dy*dy+dz*dz),x1,y1,z1);			
> 			return(0);
> 			
> 			}
> 		else{ 
> 			//printf(" Node  inprecipitate %e (%e %e %e) \n ",sqrt(dx*dx+dy*dy+dz*dz),x1,y1,z1)
> 			return(1);
> 	}
>         
>         
> }
> 
> 
> 
> int NodeinPrecipitate(real8 R,real8 xp,real8 yp,real8 zp, real8 x1,real8 y1,real8 z1,real8 x2,real8 y2,real8 z2)
> {
> 		real8   dx, dy, dz;
> 		real8	a,b,c,t1,t2,det;
> /* we solve x^2+y^2+z^2=R^2, x=x1+tdx..,dx=x2-x1 */		
> 		dx=x2-x1;
> 		dy=y2-y1;
> 		dz=z2-z1;
>         a=dx*dx+dy*dy+dz*dz;
> 		b=2.0*((x1-xp)*dx+(y1-yp)*dy +(z1-zp)*dz);
> 		c=xp*xp+yp*yp+zp*zp+x1*x1+y1*y1+z1*z1-2.0*(x1*xp+y1*yp+z1*zp)-R*R;
> 		
> 		det=b*b-4.0*a*c;
> 		
> 		if(det<0){ 
> 					
> 			return(0);
> 			
> 			}
> 		else{ 
> 			
> 		//t1=(-b+sqrt(det))/(2.0*a);
> 		//t2=(-b-sqrt(det))/(2.0*a);
> 		//if((t1<=1 && t1>0) || (t2<=1 && t2>0) ){
> 		////printf(" det,t1,t2 %e %e %e  \n ",det, t1,t2);	
> 		//return(1);	
> 		//}
> 		//else{
> 			
>         //return(0);
> 	//}
> 	  return(1);
> 	}
>         
>         
> }
> 
> 
> 
> 
1679a2165,2201
>  *      Function:     GetFreePrecipitateTag
>  *      Description:  Return a free node tag.  If possible, return a
>  *                    recycled tag. If there are no recycled tags available
>  *                    return a new tag at the end of the nodeKeys array.
>  *                    If there are no tags available there either, expand
>  *                    nodeKeys first, then return new tag.
>  *
>  *------------------------------------------------------------------------*/
> int GetFreePrecipitateTag(Home_t *home) 
> {
>         int tagIndex;
> 
>         tagIndex = GetRecycledPrecipitateTag(home);
> 
>         if (tagIndex < 0) {
> 
>             if (home->newPrecipitateKeyPtr == home->newPrecipitateKeyMax) {
>                 home->newPrecipitateKeyMax += NEW_PRECIPITATEKEY_INC;
>                 home->precipitateKeys = (Precipitate_t **) realloc(home->precipitateKeys,
>                 home->newPrecipitateKeyMax * sizeof(Precipitate_t *));
>             }
> 
>             tagIndex = home->newPrecipitateKeyPtr++;
>         }
> 
>         return(tagIndex);
> }
> 
> 
> 
> 
> 
> 
> 
> 
> /*-------------------------------------------------------------------------
>  *
1721a2244,2292
> 
> 
> 
> /*-------------------------------------------------------------------------
>  *
>  *      Function:    GetRecycledPrecipitateTag
>  *      Description: Return the lowest available recycled precipitate tag
>  *                   from the recycled tag heap.  
>  *
>  *      Returns:  -1 if no recycled tag is available, otherwise returns the
>  *                lowest available recycled tag index.
>  *                
>  *------------------------------------------------------------------------*/
> int GetRecycledPrecipitateTag(Home_t *home)
> {
>         int tagIndex;
> 
>         tagIndex = HeapRemove(home->recycledPrecipitateHeap,
>                               &home->recycledPrecipitateHeapEnts);
> 
>         return(tagIndex);
> }
> 
> 
> /*-------------------------------------------------------------------------
>  *
>  *      Function:    RecyclePrecipitateTag
>  *      Description: Add the specified tag to the list of available
>  *                   recycled Precipitates.  List is maintained as a heap.
>  *
>  *      Arguments:
>  *          tagIndex  Index of the local node tag to be recycled.
>  *
>  *------------------------------------------------------------------------*/
> void RecyclePrecipitateTag(Home_t *home, int tagIndex)
> {
> /*
>  *      Add the indicated tag to the recycle heap.  If the heap is
>  *      not large enough it's size will automatically be increased
>  */
>         HeapAdd(&home->recycledPrecipitateHeap, &home->recycledPrecipitateHeapSize,
>                 &home->recycledPrecipitateHeapEnts, tagIndex);
> 
>         return;
> }
> 
> 
> 
> 
2020c2591,2592
<                 nx, ny, nz);
---
>                 nx, ny, nz,
>                  0,0,0);
2119d2690
< 
2180a2752
>                   0.0,0.0,0.0,
2186c2758
< 
---
> //printf("Change connection %d tag ( %d %d ), (%d %d) to (%d %d)   \n",home->myDomain,node1->myTag.domainID,node1->myTag.index,tag2->domainID,tag2->index,tag3->domainID,tag3->index);	
2404c2976,2977
<                   nx,ny,nz);
---
>                  nx,ny,nz,
>                   0,0,0);   /* Integration values (AL) */
2447d3019
< 
2482c3054,3055
<             0.0,0.0,0.0);  /* nx, ny, nz */
---
>             0.0,0.0,0.0,	/* nx, ny, nz */
>             0.0, 0.0, 0.0);  /* nodedx, nodedy, nodedz AL */
2487a3061,3063
> 
> 
> 
2600c3176,3177
<                 0.0,0.0,0.0);  /* nx, ny, nz */
---
>                 0.0,0.0,0.0,/* nx, ny, nz */
>                 0.0,0.0,0.0);  
2676,2677c3253,3255
<                 f2x, f2y, f2z);  /* nx, ny, nz comandeered for */
<                                  /* 2nd set of force values    */
---
>                 f2x, f2y, f2z,/* nx, ny, nz comandeered for */
>                 0.0, 0.0, 0.0); /* 2nd set of force values    */
>                                  
2837c3415,3416
<                 newPlane[X], newPlane[Y], newPlane[Z]);
---
>                 newPlane[X], newPlane[Y], newPlane[Z],
>                  0, 0, 0);
2901a3481
>                 0.0, 0.0, 0.0,
2907a3488
> 		printf(" repositionnode  mydomain %d tag (%d,%d) \n",home->myDomain,tag->domainID,tag->index);
2919a3501
> 
2952c3534,3535
<             real8 nx, real8 ny, real8 nz) 
---
>             real8 nx, real8 ny, real8 nz,
>              real8 nodedx, real8 nodedy, real8 nodedz) 
2982a3566,3568
>         op->nodedx = nodedx;
>         op->nodedy = nodedy;
>         op->nodedz = nodedz;  
2986a3573
> 
