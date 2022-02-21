90c90
<                          int dislocType)
---
>                          int dislocType, int nodezero)
108a109,129
>  /*int      ic, np, ip, id0, i;
>         int      nplane, indp, indb, indf, inds, indr;
>         int      newNodeIndex, lastBlock;
>         int      startRemeshIndex = 0;
>         real8    inv3, inv6, invsq2, sq2over3;
>         real8    xp[6], yp[6], zp[6], cubeSize;
>         real8    tnx[4], tny[4], tnz[4], burg[12][3];
>         Node_t   *node;
>         Param_t  *param;*/
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
525c546
<                 WriteInitialNodeData(home, inData, lastBlock);
---
>                 WriteInitialNodeData(home, inData, lastBlock,nodezero);
531a553,635
>         
>         
>         
>         
>       //for (ic = 0; ic < numChains; ic++) {
> 
>             //np = 3;    /* number of points along 1 chain */
>             //newNodeIndex = inData->nodeCount;
>             //inData->nodeCount += np;
>             //inData->node = (Node_t *)realloc(inData->node,
>                            //inData->nodeCount * sizeof(Node_t));
>             //memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);
> 
> ///*
>  //*          plane normal cycle 4, BV cycle 12, 60deg line sense 24,
>  //*          reflection 48
>  //*/
>             //indp = ic%4;                /* index of plane */
>             //indb = indp*3+ic%3;         /* index of BV */
>             //indf = ((ic-ic%12)/12)%2;   /* index of alternative line sense */
>             //inds = indp*3+(ic+indf+1)%3;        /* index of line sense */
>             //indr = 1-2*(((ic-ic%24)/24)%2);     /* sign of reflection */
> 
>             //xp[0] = (randm(&seed)-0.5)*cubeSize;
>             //yp[0] = (randm(&seed)-0.5)*cubeSize;
>             //zp[0] = (randm(&seed)-0.5)*cubeSize;
> 
> ///*
>  //*          shift along neighboring 60 degree BV
>  //*/
>             //xp[1] = xp[0] + indr * cubeSize * burg[inds][0];
>             //yp[1] = yp[0] + indr * cubeSize * burg[inds][1];
>             //zp[1] = zp[0] + indr * cubeSize * burg[inds][2];
> 
>             //xp[2] = xp[1] + indr * cubeSize * burg[inds][0];
>             //yp[2] = yp[1] + indr * cubeSize * burg[inds][1];
>             //zp[2] = zp[1] + indr * cubeSize * burg[inds][2];
> 
> ///*
>  //*          PBC
>  //*/        
>             //for (i = 0; i < 3; i++) {
>                 //if (xp[i] < -cubeSize/2) xp[i] += cubeSize;
>                 //if (yp[i] < -cubeSize/2) yp[i] += cubeSize;
>                 //if (zp[i] < -cubeSize/2) zp[i] += cubeSize;
>             //}
> ///*
>             //printf("ic indp indb indf inds indr %d %d %d %d %d %d\n",
>                    //ic,indp,indb,indf,inds,indr);
> //*/
>             //for (ip = 0; ip < np; ip++) {
>                 //node = &inData->node[ip+id0];
> 
>                 //node->x = xp[ip];
>                 //node->y = yp[ip];
>                 //node->z = zp[ip];
>                 //node->constraint = PINNED_NODE; 
>                
>                 //node->myTag.domainID = dislocType;
>                 //node->myTag.index = ip+id0;
> 
>                 //AllocNodeArms(node, 2);
> 
>                 //node->nbrTag[0].domainID = dislocType;
>                 //node->nbrTag[0].index = (ip-1+np)%np+id0;
>                 //node->burgX[0] = burg[indb][0];
>                 //node->burgY[0] = burg[indb][1];
>                 //node->burgZ[0] = burg[indb][2];
>                 //node->nx[0] = tnx[indp];
>                 //node->ny[0] = tny[indp];
>                 //node->nz[0] = tnz[indp];
> 
>                 //node->nbrTag[1].domainID = dislocType;
>                 //node->nbrTag[1].index = (ip+1+np)%np+id0;
>                 //node->burgX[1] = -burg[indb][0];
>                 //node->burgY[1] = -burg[indb][1];
>                 //node->burgZ[1] = -burg[indb][2];
>                 //node->nx[1] = tnx[indp];
>                 //node->ny[1] = tny[indp];
>                 //node->nz[1] = tnz[indp];
>             //}   
>         
>       
532a637,659
> ///*
>  //*          The initial segments created are not necessarily limited to
>  //*          param->maxSegLen, so a call to InitRemesh() is needed to
>  //*          chop any excessively long segments into proper lengths.
>  //*          Then, if the count of nodes currently contained in memoru
>  //*          exceeds the threshhold write the current block of nodal data
>  //*          to the file.
>  //*/
>             //InitRemesh(inData, dislocType, startRemeshIndex);
> 
> ///*
>  //*          We won't write nodal data out here, but rather after the
>  //*          hexagonal loops are added below.
>  //*/
>             //id0 = inData->nodeCount;
>             //startRemeshIndex = id0;
> 
>         //} /* end of chains */
>    
>    
>    
>    
>    
533a661,664
>         
>        
>         
>       
551c682
<                        int dislocType)
---
>                        int dislocType, int nodezero)
734a866,868
>             //xp[1] = 0.5*cubeSize;
>             //yp[1] = 0.5*cubeSize;
>             //zp[1] = 0.5*cubeSize;
745a880,888
> 
> 
> //Let's modify this to produce more uniform wall of screw dislocations (like in 2D)
> 
> 
> 
> 
> 
> 
800c943
<                 WriteInitialNodeData(home, inData, lastBlock);
---
>                 WriteInitialNodeData(home, inData, lastBlock,nodezero);
832c975
<                           real8 *totDislocLen, int dislocType)
---
>                           real8 *totDislocLen, int dislocType,int nodezero)
965c1108,1109
<                 node->constraint = PINNED_NODE;
---
>                 node->constraint = PINNED_NODE; 
>                
1067a1212
>                 /*  node->constraint = PINNED_NODE; (AL)  */
1115c1260
<                 WriteInitialNodeData(home, inData, lastBlock);
---
>                 WriteInitialNodeData(home, inData, lastBlock,nodezero);
1144c1289
<                      int dislocType)
---
>                      int dislocType,int nodezero)
1324c1469
<                 WriteInitialNodeData(home, inData, lastBlock);
---
>                 WriteInitialNodeData(home, inData, lastBlock,nodezero);
1337d1481
< 
1340,1344c1484,1485
<  *            @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
<  *            @                                    @
<  *            @ THIS FUNCTION IS NOT YET COMPLETE! @
<  *            @                                    @
<  *            @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
---
>  *      Function:     CreateFCC2DConfig
>  *      Description:  Generates random paraller edge dislocations in random places .  (AL)
1346,1348c1487,1731
<  *      Function:     CreateFCCPerfectLoop
<  *      Description:  Generate an initial set of nodal data consisting
<  *                    of Glissile perfect loops for FCC.  (Masato Hiratani)
---
>  *      Arguments:
>  *          cubeLength  Length of cubic problem space in a single
>  *                      dimension (units of b)
>  *          numChains   Number of chains to create
>  *          seed        Seed value for random number generator
>  *
>  *-------------------------------------------------------------------------*/
> //void CreateFCC2DConfig(Home_t *home, InData_t *inData, int cubeLength,
>                  //int numChains, int seed, real8 *totDislocLen,
>                  //int dislocType,int nodezero)
> //{
>         //int           i,ic, ip, np, id0, lastBlock, signFactor;
>         //int           newNodeIndex, startRemeshIndex;
>         //int           ldIndex, gpIndex, burgIndex, nbr1Index, nbr2Index;
>         //real8         posFactor, cubeSize;
>         //real8         xp[3], yp[3], zp[3];
>         //Param_t       *param;
>         //Node_t        *node;
>         //// only one slip system b=[-1 1 0], n=[1 1 1], t=[-1 -1 2] FCC (AL)
>         //static real8  lineDir[3][3] = {
>                       //{-0.40825, -0.40825, 0.81650},
>                       //{-0.40825, -0.40825, 0.81650},
>                       //{-0.40825, -0.40825, 0.81650}};
>         //static real8  burg[4][3] = {
>                       //{ -0.70711,  0.70711,  0.0},
>                       //{ -0.70711,  0.70711,  0.0},
>                       //{ -0.70711,  0.70711,  0.0},
>                       //{-0.70711,  0.70711,  0.0}};
>         //static real8  glidePlane[12][3] = {
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [001] b [111]  */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [001] b [11-1] */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [001] b [1-11] */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [001] b [-111] */
>                       //{ 0.57735, 0.57735,  0.57735},  /* ldir [010] b [111]  */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [010] b [11-1] */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [010] b [1-11] */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [010] b [-111] */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [100] b [111]  */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [100] b [11-1] */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [100] b [1-11] */
>                       //{  0.57735, 0.57735,  0.57735}}; /* ldir [100] b [-111] */
> 
> 
>         //if (numChains <= 0) {
>             //Fatal("%s: numChains is %d, but must be > 0.\n",
>                   //"CreateEdges", numChains);
>         //}
> 
>         //param = inData->param;
>         //cubeSize = (real8)cubeLength;
>         //id0 = 0;
>         //inData->nodeCount = 0;
>         //startRemeshIndex = 0;
>         //dislocType=4;
>         ////posFactor = 0.333 * cubeSize;
>         //posFactor = 0.333 * cubeSize;
> 
> //
>  //*      Create the specified number of chains.
>  //*/
>         //for (ic = 0; ic < numChains; ic++) {
> 
>             //gpIndex = ic % 12; 
>             //burgIndex = gpIndex % 4;
>             //ldIndex = gpIndex / 4;
>             
> //
>  //*          First 12 burgers vector/normal sets use positive line
>  //*          sense, next set of 12 uses opposite line sense, and
>  //*          so on.
>  //*/
>             //signFactor = ((ic / 12) & 0x01) ? -1 : 1;
>             
>             //np = 3;
>             //newNodeIndex = inData->nodeCount;
>             //inData->nodeCount += np;
> 
>             //inData->node = (Node_t *)realloc(inData->node,
>                            //inData->nodeCount * sizeof(Node_t));
>             //memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);
> 
> //
>  //*          Set up 3 initial points for the line.  Point 1 is a base position
>  //*          at a random location, point 0 is in the negative direction along
>  //*          the line and point 2 is in the positive direction along the line.
>  //*/
>             //xp[1] = (randm(&seed)-0.5)*cubeSize;
>             //yp[1] = (randm(&seed)-0.5)*cubeSize;
>             //zp[1] = (randm(&seed)-0.5)*cubeSize;
> 
>             //xp[0] = xp[1] - (posFactor * signFactor * lineDir[ldIndex][X]);
>             //yp[0] = yp[1] - (posFactor * signFactor * lineDir[ldIndex][Y]);
>             //zp[0] = zp[1] - (posFactor * signFactor * lineDir[ldIndex][Z]);
> 
>             //xp[2] = xp[1] + (posFactor * signFactor * lineDir[ldIndex][X]);
>             //yp[2] = yp[1] + (posFactor * signFactor * lineDir[ldIndex][Y]);
>             //zp[2] = zp[1] + (posFactor * signFactor * lineDir[ldIndex][Z]);
>             
>             
>                        
> //
>  //*          PBC
>  //*/        
>             ////for (i = 0; i < 3; i++) {
>                 ////if (xp[i]<-cubeSize/2) xp[i] += cubeSize;
>                 ////if (yp[i]<-cubeSize/2) yp[i] += cubeSize;
>                 ////if (zp[i]<-cubeSize/2) zp[i] += cubeSize;
>                 ////printf("x,y,z %e %e %e \n",xp[i],yp[i],zp[i]);
>             ////}
>             
>             
>              ////for (i = 0; i < 3; i++) {
>                 ////if (xp[i]>cubeSize/2) xp[i] += -cubeSize;
>                 ////if (yp[i]>cubeSize/2) yp[i] += -cubeSize;
>                 ////if (zp[i]>cubeSize/2) zp[i] += -cubeSize;
>                 ////printf("x,y,z %e %e %e \n",xp[i],yp[i],zp[i]);
>             ////}
>             
>             
>             
> 
> //
>  //*          Loop over the points and set up the nodes, link them to 
>  //*          the neighbor nodes, etc.
>  //*/
>             ////for (ip = 0; ip < np; ip++) {
> 
>                 ////node = &inData->node[ip+id0];
> 
>                 ////node->x = xp[ip];
>                 ////node->y = yp[ip];
>                 ////node->z = zp[ip];
> 
>                 ////node->constraint = UNCONSTRAINED;
>                 ////node->myTag.domainID = dislocType;
>                 ////node->myTag.index = ip+id0;
> 
>                 ////AllocNodeArms(node, 2);
> 
>                 ////if ((nbr1Index = ip + 1) >= np) nbr1Index = 0;
>                 ////if ((nbr2Index= ip - 1) < 0) nbr2Index = np - 1;
> 
>                 ////node->nbrTag[0].domainID = dislocType;
>                 ////node->nbrTag[0].index = id0 + nbr1Index;
>                 ////node->burgX[0] = burg[burgIndex][0];
>                 ////node->burgY[0] = burg[burgIndex][1];
>                 ////node->burgZ[0] = burg[burgIndex][2];
>                 ////node->nx[0] = glidePlane[gpIndex][X];
>                 ////node->ny[0] = glidePlane[gpIndex][Y];
>                 ////node->nz[0] = glidePlane[gpIndex][Z];
>             
>                 ////node->nbrTag[1].domainID = dislocType;
>                 ////node->nbrTag[1].index = id0 + nbr2Index;
>                 ////node->burgX[1] = -burg[burgIndex][0];
>                 ////node->burgY[1] = -burg[burgIndex][1];
>                 ////node->burgZ[1] = -burg[burgIndex][2];
>                 ////node->nx[1] = glidePlane[gpIndex][X];
>                 ////node->ny[1] = glidePlane[gpIndex][Y];
>                 ////node->nz[1] = glidePlane[gpIndex][Z];
> 
>             ////}
>             
>             
>             //for (ip = 0; ip < np; ip++) {
> 
>                 //node = &inData->node[ip+id0];
> 
>                 //node->x = xp[ip];
>                 //node->y = yp[ip];
>                 //node->z = zp[ip];
>                 //node->constraint = UNCONSTRAINED;
>                 //node->myTag.domainID = dislocType;
>                 //node->myTag.index = ip+id0;
> 
>                 //AllocNodeArms(node, 2);
> 
>                 //node->nbrTag[0].domainID = dislocType;
>                 //node->nbrTag[0].index = (ip-1+np)%np+id0;
> 				//node->burgX[0] = burg[burgIndex][0];
>                 //node->burgY[0] = burg[burgIndex][1];
>                 //node->burgZ[0] = burg[burgIndex][2];
>                 //node->nx[0] = glidePlane[gpIndex][X];
>                 //node->ny[0] = glidePlane[gpIndex][Y];
>                 //node->nz[0] = glidePlane[gpIndex][Z];
> 
>                 //node->nbrTag[1].domainID = dislocType;
>                 //node->nbrTag[1].index = (ip+1+np)%np+id0;
>                 //node->burgX[1] = -burg[burgIndex][0];
>                 //node->burgY[1] = -burg[burgIndex][1];
>                 //node->burgZ[1] = -burg[burgIndex][2];
>                 //node->nx[1] = glidePlane[gpIndex][X];
>                 //node->ny[1] = glidePlane[gpIndex][Y];
>                 //node->nz[1] = glidePlane[gpIndex][Z];
> //
>                 //printf("node(%d,%d) burg=(%f %f %f) (%f %f %f)\n",
>                        //node->myTag.domainID, node->myTag.index, 
>                        //node->burgX[0],node->burgY[0],node->burgZ[0],
>                        //node->burgX[1],node->burgY[1],node->burgZ[1]);
>                 //printf("node(%d,%d) normal=(%f %f %f) (%f %f %f)\n",
>                        //node->myTag.domainID, node->myTag.index, 
>                        //node->nx[0],node->ny[0],node->nz[0],
>                        //node->nx[1],node->ny[1],node->nz[1]);
>  //*/           
>             //}
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
> //
>  //*          The initial segments created are not necessarily limited to
>  //*          param->maxSegLen, so a call to InitRemesh() is needed to
>  //*          chop any excessively long segments into proper lengths.
>  //*          When we've generated the nodal data for the final chain,
>  //*          write the current block of nodal data to the file.
>  //*/
>             //InitRemesh(inData, dislocType, startRemeshIndex);
> 
>             //lastBlock = (ic == numChains - 1);
>             //if (lastBlock) {
>                 //IncDislocationDensity(inData, totDislocLen);
>                 //param->nodeCount = inData->nodeCount;
>                 //WriteInitialNodeData(home, inData, lastBlock,nodezero);
>                 //FreeInNodeArray(inData, inData->nodeCount);
>                 //inData->nodeCount = 0;
>             //}
> 
>             //id0 = inData->nodeCount;
>             //startRemeshIndex = id0;
>         //}
>         //return;
> //}
> 
> 
> /*---------------------------------------------------------------------------
>  *
>  *      Function:     CreateFCC2DConfig
>  *      Description:  Generates random paraller edge dislocations in random places .  (AL)
1351,1354c1734,1737
<  *          cubeLength Length of cubic problem space in a single
<  *                     dimension (units of b)
<  *          numChains  Number of chains to create
<  *          seed       Seed value for random number generator
---
>  *          cubeLength  Length of cubic problem space in a single
>  *                      dimension (units of b)
>  *          numChains   Number of chains to create
>  *          seed        Seed value for random number generator
1357,1359c1740,1742
< void CreateFCCPerfectLoop(Home_t *home, InData_t *inData, int cubeLength,
<                           int numChains, int seed, real8 *totDislocLen,
<                           int dislocType)
---
> void CreateFCC2DConfig(Home_t *home, InData_t *inData, int cubeLength,
>                      int numChains, int seed, real8 *totDislocLen,
>                      int dislocType,int nodezero)
1361,1363c1744,1745
<         int      nplane, indexplane, ic, np, ip, id0, i;
<         int      inds11, inds12, inds21, inds22;
<         int      indp, indp1, indp2, indb1, indb2, indr;
---
>         int      ic, np, ip, id0, i;
>         int      nplane, indp, indb, indf, inds, indr;
1366,1371c1748,1750
<         real8    prisml, invsq2, cubeSize;
<         real8    xp[4], yp[4], zp[4];
<         real8    tempbx,tempby,tempbz,tnx[4], tny[4], tnz[4];
<         real8    sxp,syp,szp, bvx, bvy, bvz, siz;
<         real8    burg[12][3];
<         real8    vec1[3], vec2[3], tmpVec[3];
---
>         real8    invsq2;
>         real8    xp[3], yp[3], zp[3];
>         real8    tnx[4], tny[4], tnz[4], burg[12][3], cubeSize;
1374a1754,1757
> static real8  lineDir[3][3] = {
>                       {-0.40825, -0.40825, 0.81650},
>                       {-0.40825, -0.40825, 0.81650},
>                       {-0.40825, -0.40825, 0.81650}};
1376,1377c1759
<         invsq2 = 0.70710678;
<         prisml = 0.3535533*cubeLength;        /* size of loop */
---
>         invsq2 = 0.70710678118;
1382c1764
<                   "CreateFCCPerfectLoop", numChains);
---
>                   "CreateFCCConfig", numChains);
1384c1766
<     
---
>      
1392a1775
> 		dislocType=4;
1400c1783
<         tnx[0] = -1;
---
>         tnx[0] = 1;
1402c1785
<         tnz[0] = -1;
---
>         tnz[0] = 1;
1405,1406c1788,1789
<         tny[1] = -1;
<         tnz[1] = -1; 
---
>         tny[1] = 1;
>         tnz[1] = 1; 
1408,1409c1791,1792
<         tnx[2] = -1;
<         tny[2] = -1;
---
>         tnx[2] = 1;
>         tny[2] = 1;
1414a1798,1818
>         
>    //burg[4][3] = {
>                       //{ -0.70711,  0.70711,  0.0},
>                       //{ -0.70711,  0.70711,  0.0},
>                       //{ -0.70711,  0.70711,  0.0},
>                       //{-0.70711,  0.70711,  0.0}};
>         //static real8  glidePlane[12][3] = {
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [001] b [111]  */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [001] b [11-1] */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [001] b [1-11] */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [001] b [-111] */
>                       //{ 0.57735, 0.57735,  0.57735},  /* ldir [010] b [111]  */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [010] b [11-1] */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [010] b [1-11] */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [010] b [-111] */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [100] b [111]  */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [100] b [11-1] */
>                       //{  0.57735, 0.57735,  0.57735},  /* ldir [100] b [1-11] */
>                       //{  0.57735, 0.57735,  0.57735}}     
>         
>         
1421c1825
<             burg[3*i][0] = 0; 
---
>             burg[3*i][0] = -invsq2*tnz[i]; 
1423c1827
<             burg[3*i][2] = -invsq2*tnz[i]; 
---
>             burg[3*i][2] = 0; 
1426,1427c1830,1831
<             burg[1+3*i][1] = 0;
<             burg[1+3*i][2] = invsq2*tnz[i]; 
---
>             burg[1+3*i][1] = invsq2*tnz[i];
>             burg[1+3*i][2] = 0; 
1429,1430c1833,1834
<             burg[2+3*i][0] = invsq2*tnx[i];
<             burg[2+3*i][1] = -invsq2*tny[i];
---
>             burg[2+3*i][0] = -invsq2*tnx[i];
>             burg[2+3*i][1] = invsq2*tny[i];
1441c1845,1846
<             np = 4;    /* number of points along 1 chain */
---
>             np = 3;    /* number of points along 1 chain */
> 
1449c1854,1855
<  *          BV cycle 6, determination of zone axis
---
>  *          plane normal cycle 4, BV cycle 12, 60deg line sense 24,
>  *          reflection 48
1451,1502c1857,1861
<             indp = ic%6;  /* zone axis */
<             if (indp == 0) { indp1 = 0; indp2 = 1;}
<             if (indp == 1) { indp1 = 0; indp2 = 2;}
<             if (indp == 2) { indp1 = 0; indp2 = 3;}
<             if (indp == 3) { indp1 = 1; indp2 = 2;}
<             if (indp == 4) { indp1 = 1; indp2 = 3;}
<             if (indp == 5) { indp1 = 2; indp2 = 3;}
< /*
<             printf("n1= %e %e %e\n", tnx[indp1], tny[indp1], tnz[indp1]);
<             printf("n2= %e %e %e\n", tnx[indp2], tny[indp2], tnz[indp2]);
< */
<             vec1[0] = tnx[indp1];
<             vec1[1] = tny[indp1];
<             vec1[2] = tnz[indp1];
< 
<             vec2[0] = tnx[indp2];
<             vec2[1] = tny[indp2];
<             vec2[2] = tnz[indp2];
< 
<             NormalizedCrossVector(vec1, vec2, tmpVec);
< 
<             tempbx = tmpVec[0];
<             tempby = tmpVec[1];
<             tempbz = tmpVec[2];
< 
<             if (tempbx == 0) {
<                 indb1  = indp1*3;
<                 inds11 = indb1+1;
<                 inds12 = indb1+2;
<                 indb2  = indp2*3;
<                 inds21 = indb2+1;
<                 inds22 = indb2+2;
<             } else if (tempby == 0) {
<                 indb1  = indp1*3+1;
<                 inds11 = indb1+1;
<                 inds12 = indb1-1;
<                 indb2  = indp2*3+1;
<                 inds21 = indb2+1;
<                 inds22 = indb2-1;
<             } else if (tempbz == 0) {
<                 indb1  = indp1*3+2;
<                 inds11 = indb1-2;
<                 inds12 = indb1-1;
<                 indb2  = indp2*3+2;
<                 inds21 = indb2-2;
<                 inds22 = indb2-1;
<             } else {    
<                 printf("n1xn2= %e %e %e\n", tempbx, tempby, tempbz);
<                 Fatal("CreateFCCPerfectLoop: wrong zone BV");
<             }
< 
<             indr = 1-2*(((ic-ic%6)/6)%2); /* sign change: I-loop & V-loop */
---
>             indp = ic%4;              /* index of plane */
>             indb = indp*3+ic%3;       /* index of BV */
>             indf = ((ic-ic%12)/12)%2; /* index of alternative line sense */
>             inds = indp*3+(ic+indf+1)%3;    /* index of line sense */
>             indr = 1-2*(((ic-ic%24)/24)%2); /* sign of reflection */
1504,1506c1863,1865
<             xp[0] = (randm(&seed)-0.5)*cubeSize;
<             yp[0] = (randm(&seed)-0.5)*cubeSize;
<             zp[0] = (randm(&seed)-0.5)*cubeSize;
---
>             xp[0] = 0.01;//(randm(&seed)-0.5)*cubeSize;
>             yp[0] = 0.01;//(randm(&seed)-0.5)*cubeSize;
>             zp[0] = 0.01;//(randm(&seed)-0.5)*cubeSize;
1510,1548c1869,1886
<  */
<             if (indr > 0) {
<                 xp[1] = xp[0] + indr*prisml*burg[inds11][0];
<                 yp[1] = yp[0] + indr*prisml*burg[inds11][1];
<                 zp[1] = zp[0] + indr*prisml*burg[inds11][2];
<                 xp[2] = xp[1] + indr*prisml*burg[inds12][0];
<                 yp[2] = yp[1] + indr*prisml*burg[inds12][1];
<                 zp[2] = zp[1] + indr*prisml*burg[inds12][2];
<                 sxp   = xp[2] + indr*prisml*burg[inds21][0];
<                 syp   = yp[2] + indr*prisml*burg[inds21][1];
<                 szp   = zp[2] + indr*prisml*burg[inds21][2];
<                 xp[3] = xp[2] - indr*prisml*burg[inds11][0];
<                 yp[3] = yp[2] - indr*prisml*burg[inds11][1];
<                 zp[3] = zp[2] - indr*prisml*burg[inds11][2];
<                 indexplane = indp1;
<              } else {
<                 xp[1] = xp[0] + indr*prisml*burg[inds21][0];
<                 yp[1] = yp[0] + indr*prisml*burg[inds21][1];
<                 zp[1] = zp[0] + indr*prisml*burg[inds21][2];
<                 xp[2] = xp[1] + indr*prisml*burg[inds22][0];
<                 yp[2] = yp[1] + indr*prisml*burg[inds22][1];
<                 zp[2] = zp[1] + indr*prisml*burg[inds22][2];
<                 sxp   = xp[2] + indr*prisml*burg[inds11][0];
<                 syp   = yp[2] + indr*prisml*burg[inds11][1];
<                 szp   = zp[2] + indr*prisml*burg[inds11][2];
<                 xp[3] = xp[2] - indr*prisml*burg[inds21][0];
<                 yp[3] = yp[2] - indr*prisml*burg[inds21][1];
<                 zp[3] = zp[2] - indr*prisml*burg[inds21][2];
<                 indexplane = indp2;
<             }
< 
<             GetUnitVector(1, sxp, syp, szp, xp[1], yp[1], zp[1],
<                           &bvx, &bvy, &bvz, &siz); 
< 
<             printf("ic,indb1 indr indexplane=\n");
<             printf("%d %d %d %d \n", ic, indb1, indr, indexplane);    
<             printf("n= %e %e %e\n", tnx[indexplane],
<                    tny[indexplane], tnz[indexplane]);    
<             printf("b= %e %e %e\n", bvx,bvy,bvz);
---
>  */								
>             //xp[1] = xp[0] +2.0* indr*cubeSize*lineDir[0][0];
>             //yp[1] = yp[0] +2.0*indr*cubeSize*lineDir[1][1];
>             //zp[1] = zp[0] +2.0*indr*cubeSize*lineDir[2][2];
> 
>             //xp[2] = xp[0] -2.0*indr*cubeSize*lineDir[0][0];
>             //yp[2] = yp[0] -2.0*indr*cubeSize*lineDir[1][1];
>             //zp[2] = zp[0] -2.0*indr*cubeSize*lineDir[2][2];
>             
>             xp[1] = xp[0] +indr*cubeSize*(-0.40825)*6.0;
>             yp[1] = yp[0] +indr*cubeSize*(-0.40825)*6.0;;
>             zp[1] = zp[0] +indr*cubeSize*(0.81650)*6.0;;
> 
>             xp[2] = xp[1] +indr*cubeSize*(-0.40825)*6.0;;
>             yp[2] = yp[1] +indr*cubeSize*(-0.40825)*6.0;;
>             zp[2] = zp[1] +indr*cubeSize*(0.81650)*6.0;;
>             
>             
1553c1891
<             for (i = 0; i < np; i++) {
---
>             for (i = 0; i < 3; i++) {
1557a1896
> 
1559,1560c1898,1901
<             printf("ic indp indr %d %d %d\n",ic,indp, indr);
<  */
---
>             printf("ic indp indb indf inds indr %d %d %d %d %d %d\n",
>                    ic,indp, indb, indf, inds, indr);
> */
> 
1576,1584c1917,1922
< 
<                 node->burgX[0] = bvx;
<                 node->burgY[0] = bvy;
<                 node->burgZ[0] = bvz;
< 
<                 node->nx[0] = 0.0;
<                 node->ny[0] = 0.0;
<                 node->nz[0] = 0.0;
< 
---
>                 node->burgX[0] = burg[indb][0];
>                 node->burgY[0] = burg[indb][1];
>                 node->burgZ[0] = burg[indb][2];
>                 node->nx[0] = tnx[indp]; 
>                 node->ny[0] = tny[indp];
>                 node->nz[0] = tnz[indp];
1588,1596c1926,1931
< 
<                 node->burgX[1] = -bvx;
<                 node->burgY[1] = -bvy;
<                 node->burgZ[1] = -bvz;
< 
<                 node->nx[1] = 0.0;
<                 node->ny[1] = 0.0;
<                 node->nz[1] = 0.0;
< 
---
>                 node->burgX[1] = -burg[indb][0];
>                 node->burgY[1] = -burg[indb][1];
>                 node->burgZ[1] = -burg[indb][2];
>                 node->nx[1] = tnx[indp]; 
>                 node->ny[1] = tny[indp];
>                 node->nz[1] = tnz[indp];
1606c1941
<  */   
---
>  */           
1607a1943
> 
1621c1957
<                 WriteInitialNodeData(home, inData, lastBlock);
---
>                 WriteInitialNodeData(home, inData, lastBlock,nodezero);
1632c1968,1971
< }
---
> }                     
> 
> 
> 
1637,1641c1976,1984
<  *      Function:     CreateFiniteMixedConfig
<  *      Description:  Creates random screw and edge dislocations,
<  *                    but assumes a finite problem space in at
<  *                    least 1 dimension and terminates the dislocations
<  *                    at the free surface(s) 
---
>  *            @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
>  *            @                                    @
>  *            @ THIS FUNCTION IS NOT YET COMPLETE! @
>  *            @                                    @
>  *            @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
>  *
>  *      Function:     CreateFCCPerfectLoop
>  *      Description:  Generate an initial set of nodal data consisting
>  *                    of Glissile perfect loops for FCC.  (Masato Hiratani)
1644,1650c1987,1990
<  *          cubeLength    Length of cubic problem space in a single
<  *                        dimension (units of b)
<  *          numChains     Number of chains to create
<  *          seed          Seed value for random number generator
<  *          totDislocLen  Pointer to location at which to return
<  *                        to caller the total combined length of
<  *                        all created dislocations.
---
>  *          cubeLength Length of cubic problem space in a single
>  *                     dimension (units of b)
>  *          numChains  Number of chains to create
>  *          seed       Seed value for random number generator
1653,1655c1993,1995
< void CreateFiniteMixedConfig(Home_t *home, InData_t *inData, int cubeLength,
<                              int numChains, int seed, real8 *totDislocLen,
<                              int dislocType)
---
> void CreateFCCPerfectLoop(Home_t *home, InData_t *inData, int cubeLength,
>                           int numChains, int seed, real8 *totDislocLen,
>                           int dislocType,int nodezero)
1657,1659c1997,2000
<         int      i, lineIndex, chain, baseNodeID, numPoints;
<         int      newNodeIndex, lastBlock, burgIndex;
<         int      isScrew, gpIndex, gpBurgIndex;
---
>         int      nplane, indexplane, ic, np, ip, id0, i;
>         int      inds11, inds12, inds21, inds22;
>         int      indp, indp1, indp2, indb1, indb2, indr;
>         int      newNodeIndex, lastBlock;
1661,1674c2002,2007
<         int      signFact; 
<         int      nbr1Index, nbr2Index, numConnections, numSegs;
<         int      intersectIndex1, intersectIndex2;
<         int      pbc[3];
<         real8    burgSign, vecLen;
<         real8    surfCoord, len, minLen1, minLen2;
<         real8    intersectSurfCoord1, intersectSurfCoord2;
<         real8    p0[3], newPos[3];
<         real8    range[3], rangeCntr[3];
<         real8    minCoord[3], maxCoord[3];
<         real8    intersectPos1[3], intersectPos2[3], vector[3];
<         real8    normalizedBurg[3];
<         real8    linedir[16][3], glidePlane[4][6][3], plane[3];
<         Node_t   *node;
---
>         real8    prisml, invsq2, cubeSize;
>         real8    xp[4], yp[4], zp[4];
>         real8    tempbx,tempby,tempbz,tnx[4], tny[4], tnz[4];
>         real8    sxp,syp,szp, bvx, bvy, bvz, siz;
>         real8    burg[12][3];
>         real8    vec1[3], vec2[3], tmpVec[3];
1676,1682c2009
<         
<         if (numChains <= 0) {
<             Fatal("%s: numChains is %d, but must be > 0.\n",
<                   "CreateFiniteMixedConfig", numChains);
<         }
<         
<         param = inData->param;
---
>         Node_t   *node;
1684,1686c2011,2014
<         pbc[X] = (param->xBoundType == Periodic);
<         pbc[Y] = (param->yBoundType == Periodic);
<         pbc[Z] = (param->zBoundType == Periodic);
---
>         param = inData->param;
>         invsq2 = 0.70710678;
>         prisml = 0.3535533*cubeLength;        /* size of loop */
>         cubeSize = (real8)cubeLength;
1688,1693c2016,2018
< /*
<  *      Just a sanity check...
<  */
<         if (pbc[X]*pbc[Y]*pbc[Z] != 0) {
<             Fatal("CreateFiniteMixedConfig: Requires free surfaces in at"
<                   "least one dimension\n    but all boundaries are periodic!");
---
>         if (numChains <= 0) {
>             Fatal("%s: numChains is %d, but must be a multiple of 12.\n",
>                   "CreateFCCPerfectLoop", numChains);
1695c2020
< 
---
>     
1697,1704c2022,2027
<  *      Set up an array of dislocation line directions that
<  *      are used to create the new dislocations.  There are
<  *      essentially 4 sets of line directions (1 per burgers
<  *      vector) with 4 line directions per set; the first for
<  *      the screw and the following three for edge.
<  *      The burger's vector for all 4 dislocation lines in
<  *      a set is the same as the line direction of the screw
<  *      dislocation in the group.
---
>  *      The nodal data is generated in a single block.  All nodes will
>  *      be assigned a tag.domainID equal to the dislocation type, and
>  *      the index for every node will be the node's actual index into the
>  *      block of nodes. The InitRemesh() calls also accept the dislocation
>  *      type to set the tag.domainID value for any new nodes it adds to the
>  *      node array.
1705a2029
>         inData->nburg = 12;
1708,1789c2032
<  *      Type  [1 1 1] burgers vector
<  */
<         linedir[0][0] =  0.5773503;  
<         linedir[0][1] =  0.5773503;
<         linedir[0][2] =  0.5773503; 
<         
<         linedir[1][0] = -0.8164966;
<         linedir[1][1] =  0.4082483;
<         linedir[1][2] =  0.4082483; 
<         
<         linedir[2][0] =  0.4082483;
<         linedir[2][1] = -0.8164966;
<         linedir[2][2] =  0.4082483; 
<         
<         linedir[3][0] =  0.4082483;
<         linedir[3][1] =  0.4082483; 
<         linedir[3][2] = -0.8164966;
<         
< /*
<  *      Type [-1 1 1] burgers vector
<  */
<         linedir[4][0] = -0.5773503;
<         linedir[4][1] =  0.5773503;
<         linedir[4][2] =  0.5773503; 
<         
<         linedir[5][0] =  0.8164966;
<         linedir[5][1] =  0.4082483;
<         linedir[5][2] =  0.4082483; 
<         
<         linedir[6][0] =  0.4082483;
<         linedir[6][1] =  0.8164966;
<         linedir[6][2] = -0.4082483; 
<         
<         linedir[7][0] =  0.4082483;
<         linedir[7][1] = -0.4082483; 
<         linedir[7][2] =  0.8164966;
<         
< /*
<  *      Type [1 -1 1] burgers vector
<  */
<         linedir[8][0] =  0.5773503;
<         linedir[8][1] = -0.5773503;
<         linedir[8][2] =  0.5773503; 
<         
<         linedir[9][0] =  0.4082483;
<         linedir[9][1] =  0.8164966;
<         linedir[9][2] =  0.4082483; 
<         
<         linedir[10][0] =  0.8164966;
<         linedir[10][1] =  0.4082483;
<         linedir[10][2] = -0.4082483; 
<         
<         linedir[11][0] = -0.4082483;
<         linedir[11][1] =  0.4082483; 
<         linedir[11][2] =  0.8164966;
<         
< /*
<  *      Type [1 1 -1] burgers vector
<  */
<         linedir[12][0] =  0.5773503;
<         linedir[12][1] =  0.5773503;
<         linedir[12][2] = -0.5773503; 
<         
<         linedir[13][0] = -0.4082483;
<         linedir[13][1] =  0.8164966;
<         linedir[13][2] =  0.4082483; 
<         
<         linedir[14][0] =  0.8164966;
<         linedir[14][1] = -0.4082483;
<         linedir[14][2] =  0.4082483; 
<         
<         linedir[15][0] =  0.4082483;
<         linedir[15][1] =  0.4082483; 
<         linedir[15][2] =  0.8164966;
<         
< /*
<  *      Set up the valid glide planes for each screw burgers vector,
<  *      six glide planes per burgers vector.  For edges, glide plane
<  *      will simply be cross product between burgers vector and 
<  *      linedir.
<  *
<  *      glide planes for [1 1 1]
---
>  *      alpha, beta, gamma, and delta planes
1791,1797c2034
<         glidePlane[0][0][0] =  0.7071068;
<         glidePlane[0][0][1] = -0.7071068;
<         glidePlane[0][0][2] =  0.0000000;
< 
<         glidePlane[0][1][0] =  0.7071068;
<         glidePlane[0][1][1] =  0.0000000;
<         glidePlane[0][1][2] = -0.7071068;
---
>         nplane =  4;
1799,1801c2036,2038
<         glidePlane[0][2][0] =  0.0000000;
<         glidePlane[0][2][1] =  0.7071068;
<         glidePlane[0][2][2] = -0.7071068;
---
>         tnx[0] = -1;
>         tny[0] =  1;
>         tnz[0] = -1;
1803,1805c2040,2042
<         glidePlane[0][3][0] = -0.7071068;
<         glidePlane[0][3][1] =  0.7071068;
<         glidePlane[0][3][2] = -0.0000000;
---
>         tnx[1] =  1;
>         tny[1] = -1;
>         tnz[1] = -1; 
1807,1809c2044,2046
<         glidePlane[0][4][0] = -0.7071068;
<         glidePlane[0][4][1] = -0.0000000;
<         glidePlane[0][4][2] =  0.7071068;
---
>         tnx[2] = -1;
>         tny[2] = -1;
>         tnz[2] =  1;
1811,1813c2048,2050
<         glidePlane[0][5][0] = -0.0000000;
<         glidePlane[0][5][1] = -0.7071068;
<         glidePlane[0][5][2] =  0.7071068;
---
>         tnx[3] =  1;
>         tny[3] =  1;
>         tnz[3] =  1; 
1816c2053
<  *      glide planes for [-1 1 1]
---
>  *      BV in counter-clockwise on each plane
1818,1820c2055
<         glidePlane[1][0][0] =  0.0000000;
<         glidePlane[1][0][1] =  0.7071068;
<         glidePlane[1][0][2] = -0.7071068;
---
>         for (i = 0; i < nplane; i++) {
1822,1824c2057,2059
<         glidePlane[1][1][0] =  0.7071068;
<         glidePlane[1][1][1] =  0.0000000;
<         glidePlane[1][1][2] =  0.7071068;
---
>             burg[3*i][0] = 0; 
>             burg[3*i][1] = invsq2*tny[i];
>             burg[3*i][2] = -invsq2*tnz[i]; 
1826c2061,2913
<         glidePlane[1][2][0] =  0.0000000;
---
>             burg[1+3*i][0] = -invsq2*tnx[i];
>             burg[1+3*i][1] = 0;
>             burg[1+3*i][2] = invsq2*tnz[i]; 
> 
>             burg[2+3*i][0] = invsq2*tnx[i];
>             burg[2+3*i][1] = -invsq2*tny[i];
>             burg[2+3*i][2] = 0; 
> 
>             Normalize(&tnx[i],&tny[i],&tnz[i]);   
>         }
>     
>         id0 = 0;
>         inData->nodeCount = 0;
> 
>         for (ic = 0; ic < numChains; ic++) {
> 
>             np = 4;    /* number of points along 1 chain */
>             newNodeIndex = inData->nodeCount;
>             inData->nodeCount += np;
>             inData->node = (Node_t *) realloc(inData->node,
>                            inData->nodeCount * sizeof(Node_t));
>             memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);
> 
> /*
>  *          BV cycle 6, determination of zone axis
>  */
>             indp = ic%6;  /* zone axis */
>             if (indp == 0) { indp1 = 0; indp2 = 1;}
>             if (indp == 1) { indp1 = 0; indp2 = 2;}
>             if (indp == 2) { indp1 = 0; indp2 = 3;}
>             if (indp == 3) { indp1 = 1; indp2 = 2;}
>             if (indp == 4) { indp1 = 1; indp2 = 3;}
>             if (indp == 5) { indp1 = 2; indp2 = 3;}
> /*
>             printf("n1= %e %e %e\n", tnx[indp1], tny[indp1], tnz[indp1]);
>             printf("n2= %e %e %e\n", tnx[indp2], tny[indp2], tnz[indp2]);
> */
>             vec1[0] = tnx[indp1];
>             vec1[1] = tny[indp1];
>             vec1[2] = tnz[indp1];
> 
>             vec2[0] = tnx[indp2];
>             vec2[1] = tny[indp2];
>             vec2[2] = tnz[indp2];
> 
>             NormalizedCrossVector(vec1, vec2, tmpVec);
> 
>             tempbx = tmpVec[0];
>             tempby = tmpVec[1];
>             tempbz = tmpVec[2];
> 
>             if (tempbx == 0) {
>                 indb1  = indp1*3;
>                 inds11 = indb1+1;
>                 inds12 = indb1+2;
>                 indb2  = indp2*3;
>                 inds21 = indb2+1;
>                 inds22 = indb2+2;
>             } else if (tempby == 0) {
>                 indb1  = indp1*3+1;
>                 inds11 = indb1+1;
>                 inds12 = indb1-1;
>                 indb2  = indp2*3+1;
>                 inds21 = indb2+1;
>                 inds22 = indb2-1;
>             } else if (tempbz == 0) {
>                 indb1  = indp1*3+2;
>                 inds11 = indb1-2;
>                 inds12 = indb1-1;
>                 indb2  = indp2*3+2;
>                 inds21 = indb2-2;
>                 inds22 = indb2-1;
>             } else {    
>                 printf("n1xn2= %e %e %e\n", tempbx, tempby, tempbz);
>                 Fatal("CreateFCCPerfectLoop: wrong zone BV");
>             }
> 
>             indr = 1-2*(((ic-ic%6)/6)%2); /* sign change: I-loop & V-loop */
> 
>             xp[0] = (randm(&seed)-0.5)*cubeSize;
>             yp[0] = (randm(&seed)-0.5)*cubeSize;
>             zp[0] = (randm(&seed)-0.5)*cubeSize;
> 
> /*
>  *          shift along neighboring 60 degree BV
>  */
>             if (indr > 0) {
>                 xp[1] = xp[0] + indr*prisml*burg[inds11][0];
>                 yp[1] = yp[0] + indr*prisml*burg[inds11][1];
>                 zp[1] = zp[0] + indr*prisml*burg[inds11][2];
>                 xp[2] = xp[1] + indr*prisml*burg[inds12][0];
>                 yp[2] = yp[1] + indr*prisml*burg[inds12][1];
>                 zp[2] = zp[1] + indr*prisml*burg[inds12][2];
>                 sxp   = xp[2] + indr*prisml*burg[inds21][0];
>                 syp   = yp[2] + indr*prisml*burg[inds21][1];
>                 szp   = zp[2] + indr*prisml*burg[inds21][2];
>                 xp[3] = xp[2] - indr*prisml*burg[inds11][0];
>                 yp[3] = yp[2] - indr*prisml*burg[inds11][1];
>                 zp[3] = zp[2] - indr*prisml*burg[inds11][2];
>                 indexplane = indp1;
>              } else {
>                 xp[1] = xp[0] + indr*prisml*burg[inds21][0];
>                 yp[1] = yp[0] + indr*prisml*burg[inds21][1];
>                 zp[1] = zp[0] + indr*prisml*burg[inds21][2];
>                 xp[2] = xp[1] + indr*prisml*burg[inds22][0];
>                 yp[2] = yp[1] + indr*prisml*burg[inds22][1];
>                 zp[2] = zp[1] + indr*prisml*burg[inds22][2];
>                 sxp   = xp[2] + indr*prisml*burg[inds11][0];
>                 syp   = yp[2] + indr*prisml*burg[inds11][1];
>                 szp   = zp[2] + indr*prisml*burg[inds11][2];
>                 xp[3] = xp[2] - indr*prisml*burg[inds21][0];
>                 yp[3] = yp[2] - indr*prisml*burg[inds21][1];
>                 zp[3] = zp[2] - indr*prisml*burg[inds21][2];
>                 indexplane = indp2;
>             }
> 
>             GetUnitVector(1, sxp, syp, szp, xp[1], yp[1], zp[1],
>                           &bvx, &bvy, &bvz, &siz); 
> 
>             printf("ic,indb1 indr indexplane=\n");
>             printf("%d %d %d %d \n", ic, indb1, indr, indexplane);    
>             printf("n= %e %e %e\n", tnx[indexplane],
>                    tny[indexplane], tnz[indexplane]);    
>             printf("b= %e %e %e\n", bvx,bvy,bvz);
> 
> /*
>  *          PBC
>  */        
>             for (i = 0; i < np; i++) {
>                 if (xp[i]<-cubeSize/2) xp[i] += cubeSize;
>                 if (yp[i]<-cubeSize/2) yp[i] += cubeSize;
>                 if (zp[i]<-cubeSize/2) zp[i] += cubeSize;
>             }
> /*
>             printf("ic indp indr %d %d %d\n",ic,indp, indr);
>  */
>             for (ip = 0; ip < np; ip++) {
> 
>                 node = &inData->node[ip+id0];
> 
>                 node->x = xp[ip];
>                 node->y = yp[ip];
>                 node->z = zp[ip];
>                 node->constraint = UNCONSTRAINED;
>                 node->myTag.domainID = dislocType;
>                 node->myTag.index = ip+id0;
> 
>                 AllocNodeArms(node, 2);
> 
>                 node->nbrTag[0].domainID = dislocType;
>                 node->nbrTag[0].index = (ip-1+np)%np+id0;
> 
>                 node->burgX[0] = bvx;
>                 node->burgY[0] = bvy;
>                 node->burgZ[0] = bvz;
> 
>                 node->nx[0] = 0.0;
>                 node->ny[0] = 0.0;
>                 node->nz[0] = 0.0;
> 
> 
>                 node->nbrTag[1].domainID = dislocType;
>                 node->nbrTag[1].index = (ip+1+np)%np+id0;
> 
>                 node->burgX[1] = -bvx;
>                 node->burgY[1] = -bvy;
>                 node->burgZ[1] = -bvz;
> 
>                 node->nx[1] = 0.0;
>                 node->ny[1] = 0.0;
>                 node->nz[1] = 0.0;
> 
> /*
>                 printf("node(%d,%d) burg=(%f %f %f) (%f %f %f)\n",
>                        node->myTag.domainID, node->myTag.index, 
>                        node->burgX[0],node->burgY[0],node->burgZ[0],
>                        node->burgX[1],node->burgY[1],node->burgZ[1]);
>                 printf("node(%d,%d) normal=(%f %f %f) (%f %f %f)\n",
>                        node->myTag.domainID, node->myTag.index, 
>                        node->nx[0],node->ny[0],node->nz[0],
>                        node->nx[1],node->ny[1],node->nz[1]);
>  */   
>             }
> /*
>  *          The initial segments created are not necessarily limited to
>  *          param->maxSegLen, so a call to InitRemesh() is needed to
>  *          chop any excessively long segments into proper lengths.
>  *          When we've generated the nodal data for the final chain,
>  *          write the block of nodal data to the file.
>  */
>             InitRemesh(inData, dislocType, startRemeshIndex);
> 
>             lastBlock = (ic == (numChains - 1));
>             if (lastBlock) {
>                 IncDislocationDensity(inData, totDislocLen);
>                 param->nodeCount = inData->nodeCount;
>                 WriteInitialNodeData(home, inData, lastBlock,nodezero);
>                 FreeInNodeArray(inData, inData->nodeCount);
>                 inData->nodeCount = 0;
>             }
> 
>             id0 = inData->nodeCount;
>             startRemeshIndex = id0;
> 
>         } /* loop over chains */
> 
>         return;
> }
> 
> 
> /*---------------------------------------------------------------------------
>  *
>  *      Function:     CreateFiniteMixedConfig
>  *      Description:  Creates random screw and edge dislocations,
>  *                    but assumes a finite problem space in at
>  *                    least 1 dimension and terminates the dislocations
>  *                    at the free surface(s) 
>  *
>  *      Arguments:
>  *          cubeLength    Length of cubic problem space in a single
>  *                        dimension (units of b)
>  *          numChains     Number of chains to create
>  *          seed          Seed value for random number generator
>  *          totDislocLen  Pointer to location at which to return
>  *                        to caller the total combined length of
>  *                        all created dislocations.
>  *
>  *-------------------------------------------------------------------------*/
> void CreateFiniteMixedConfig(Home_t *home, InData_t *inData, int cubeLength,
>                              int numChains, int seed, real8 *totDislocLen,
>                              int dislocType,int nodezero)
> {
>         int      i, lineIndex, chain, baseNodeID, numPoints;
>         int      newNodeIndex, lastBlock, burgIndex;
>         int      isScrew, gpIndex, gpBurgIndex;
>         int      startRemeshIndex = 0;
>         int      signFact; 
>         int      nbr1Index, nbr2Index, numConnections, numSegs;
>         int      intersectIndex1, intersectIndex2;
>         int      pbc[3];
>         real8    burgSign, vecLen;
>         real8    surfCoord, len, minLen1, minLen2;
>         real8    intersectSurfCoord1, intersectSurfCoord2;
>         real8    p0[3], newPos[3];
>         real8    range[3], rangeCntr[3];
>         real8    minCoord[3], maxCoord[3];
>         real8    intersectPos1[3], intersectPos2[3], vector[3];
>         real8    normalizedBurg[3];
>         real8    linedir[16][3], glidePlane[4][6][3], plane[3];
>         Node_t   *node;
>         Param_t  *param;
>         
>         if (numChains <= 0) {
>             Fatal("%s: numChains is %d, but must be > 0.\n",
>                   "CreateFiniteMixedConfig", numChains);
>         }
>         
>         param = inData->param;
> 
>         pbc[X] = (param->xBoundType == Periodic);
>         pbc[Y] = (param->yBoundType == Periodic);
>         pbc[Z] = (param->zBoundType == Periodic);
> 
> /*
>  *      Just a sanity check...
>  */
>         if (pbc[X]*pbc[Y]*pbc[Z] != 0) {
>             Fatal("CreateFiniteMixedConfig: Requires free surfaces in at"
>                   "least one dimension\n    but all boundaries are periodic!");
>         }
> 
> /*
>  *      Set up an array of dislocation line directions that
>  *      are used to create the new dislocations.  There are
>  *      essentially 4 sets of line directions (1 per burgers
>  *      vector) with 4 line directions per set; the first for
>  *      the screw and the following three for edge.
>  *      The burger's vector for all 4 dislocation lines in
>  *      a set is the same as the line direction of the screw
>  *      dislocation in the group.
>  */
> 
> /*
>  *      Type  [1 1 1] burgers vector
>  */
>         linedir[0][0] =  0.5773503;  
>         linedir[0][1] =  0.5773503;
>         linedir[0][2] =  0.5773503; 
>         
>         linedir[1][0] = -0.8164966;
>         linedir[1][1] =  0.4082483;
>         linedir[1][2] =  0.4082483; 
>         
>         linedir[2][0] =  0.4082483;
>         linedir[2][1] = -0.8164966;
>         linedir[2][2] =  0.4082483; 
>         
>         linedir[3][0] =  0.4082483;
>         linedir[3][1] =  0.4082483; 
>         linedir[3][2] = -0.8164966;
>         
> /*
>  *      Type [-1 1 1] burgers vector
>  */
>         linedir[4][0] = -0.5773503;
>         linedir[4][1] =  0.5773503;
>         linedir[4][2] =  0.5773503; 
>         
>         linedir[5][0] =  0.8164966;
>         linedir[5][1] =  0.4082483;
>         linedir[5][2] =  0.4082483; 
>         
>         linedir[6][0] =  0.4082483;
>         linedir[6][1] =  0.8164966;
>         linedir[6][2] = -0.4082483; 
>         
>         linedir[7][0] =  0.4082483;
>         linedir[7][1] = -0.4082483; 
>         linedir[7][2] =  0.8164966;
>         
> /*
>  *      Type [1 -1 1] burgers vector
>  */
>         linedir[8][0] =  0.5773503;
>         linedir[8][1] = -0.5773503;
>         linedir[8][2] =  0.5773503; 
>         
>         linedir[9][0] =  0.4082483;
>         linedir[9][1] =  0.8164966;
>         linedir[9][2] =  0.4082483; 
>         
>         linedir[10][0] =  0.8164966;
>         linedir[10][1] =  0.4082483;
>         linedir[10][2] = -0.4082483; 
>         
>         linedir[11][0] = -0.4082483;
>         linedir[11][1] =  0.4082483; 
>         linedir[11][2] =  0.8164966;
>         
> /*
>  *      Type [1 1 -1] burgers vector
>  */
>         linedir[12][0] =  0.5773503;
>         linedir[12][1] =  0.5773503;
>         linedir[12][2] = -0.5773503; 
>         
>         linedir[13][0] = -0.4082483;
>         linedir[13][1] =  0.8164966;
>         linedir[13][2] =  0.4082483; 
>         
>         linedir[14][0] =  0.8164966;
>         linedir[14][1] = -0.4082483;
>         linedir[14][2] =  0.4082483; 
>         
>         linedir[15][0] =  0.4082483;
>         linedir[15][1] =  0.4082483; 
>         linedir[15][2] =  0.8164966;
>         
> /*
>  *      Set up the valid glide planes for each screw burgers vector,
>  *      six glide planes per burgers vector.  For edges, glide plane
>  *      will simply be cross product between burgers vector and 
>  *      linedir.
>  *
>  *      glide planes for [1 1 1]
>  */
>         glidePlane[0][0][0] =  0.7071068;
>         glidePlane[0][0][1] = -0.7071068;
>         glidePlane[0][0][2] =  0.0000000;
> 
>         glidePlane[0][1][0] =  0.7071068;
>         glidePlane[0][1][1] =  0.0000000;
>         glidePlane[0][1][2] = -0.7071068;
> 
>         glidePlane[0][2][0] =  0.0000000;
>         glidePlane[0][2][1] =  0.7071068;
>         glidePlane[0][2][2] = -0.7071068;
> 
>         glidePlane[0][3][0] = -0.7071068;
>         glidePlane[0][3][1] =  0.7071068;
>         glidePlane[0][3][2] = -0.0000000;
> 
>         glidePlane[0][4][0] = -0.7071068;
>         glidePlane[0][4][1] = -0.0000000;
>         glidePlane[0][4][2] =  0.7071068;
> 
>         glidePlane[0][5][0] = -0.0000000;
>         glidePlane[0][5][1] = -0.7071068;
>         glidePlane[0][5][2] =  0.7071068;
> 
> /*
>  *      glide planes for [-1 1 1]
>  */
>         glidePlane[1][0][0] =  0.0000000;
>         glidePlane[1][0][1] =  0.7071068;
>         glidePlane[1][0][2] = -0.7071068;
> 
>         glidePlane[1][1][0] =  0.7071068;
>         glidePlane[1][1][1] =  0.0000000;
>         glidePlane[1][1][2] =  0.7071068;
> 
>         glidePlane[1][2][0] =  0.0000000;
>         glidePlane[1][2][1] =  0.7071068;
>         glidePlane[1][2][2] = -0.7071068;
> 
>         glidePlane[1][3][0] = -0.0000000;
>         glidePlane[1][3][1] = -0.7071068;
>         glidePlane[1][3][2] =  0.7071068;
> 
>         glidePlane[1][4][0] = -0.7071068;
>         glidePlane[1][4][1] = -0.0000000;
>         glidePlane[1][4][2] = -0.7071068;
> 
>         glidePlane[1][5][0] = -0.0000000;
>         glidePlane[1][5][1] = -0.7071068;
>         glidePlane[1][5][2] =  0.7071068;
> 
> /*
>  *      glide planes for [1 -1 1]
>  */
>         glidePlane[2][0][0] =  0.7071068;
>         glidePlane[2][0][1] =  0.7071068;
>         glidePlane[2][0][2] =  0.0000000;
> 
>         glidePlane[2][1][0] =  0.7071068;
>         glidePlane[2][1][1] =  0.0000000;
>         glidePlane[2][1][2] = -0.7071068;
> 
>         glidePlane[2][2][0] =  0.0000000;
>         glidePlane[2][2][1] =  0.7071068;
>         glidePlane[2][2][2] =  0.7071068;
> 
>         glidePlane[2][3][0] = -0.7071068;
>         glidePlane[2][3][1] = -0.7071068;
>         glidePlane[2][3][2] = -0.0000000;
> 
>         glidePlane[2][4][0] = -0.7071068;
>         glidePlane[2][4][1] = -0.0000000;
>         glidePlane[2][4][2] =  0.7071068;
> 
>         glidePlane[2][5][0] = -0.0000000;
>         glidePlane[2][5][1] = -0.7071068;
>         glidePlane[2][5][2] = -0.7071068;
> 
> /*
>  *      glide planes for [1 1 -1]
>  */
>         glidePlane[3][0][0] =  0.7071068;
>         glidePlane[3][0][1] = -0.7071068;
>         glidePlane[3][0][2] =  0.0000000;
> 
>         glidePlane[3][1][0] =  0.0000000;
>         glidePlane[3][1][1] =  0.7071068;
>         glidePlane[3][1][2] =  0.7071068;
> 
>         glidePlane[3][2][0] =  0.7071068;
>         glidePlane[3][2][1] =  0.0000000;
>         glidePlane[3][2][2] =  0.7071068;
> 
>         glidePlane[3][3][0] = -0.7071068;
>         glidePlane[3][3][1] =  0.7071068;
>         glidePlane[3][3][2] = -0.0000000;
> 
>         glidePlane[3][4][0] = -0.0000000;
>         glidePlane[3][4][1] = -0.7071068;
>         glidePlane[3][4][2] = -0.7071068;
> 
>         glidePlane[3][5][0] = -0.7071068;
>         glidePlane[3][5][1] = -0.0000000;
>         glidePlane[3][5][2] = -0.7071068;
> 
>         minCoord[X] = param->xBoundMin;
>         maxCoord[X] = param->xBoundMax;
>         minCoord[Y] = param->yBoundMin;
>         maxCoord[Y] = param->yBoundMax;
>         minCoord[Z] = param->zBoundMin;
>         maxCoord[Z] = param->zBoundMax;
> 
>         for (i = 0; i < 3; i++) {
>             range[i] = maxCoord[i] - minCoord[i];
>             rangeCntr[i] = 0.5 * (minCoord[i] + maxCoord[i]);
>         }
> 
> /*
>  *      Create the specified number of chains.  Anytime the number of
>  *      nodes maintained in memory exceeds the threshhold, write the
>  *      block of nodal data out to the data file.
>  */
>         baseNodeID = 0;
>         inData->nodeCount = 0;
> 
>         for (chain = 0; chain < numChains; chain++) {
>         
>             numPoints = 3;
>             lineIndex = chain % 16; 
>             burgIndex = 4 * (lineIndex / 4);
>             gpBurgIndex = lineIndex / 4;
>             gpIndex = (chain / 16) % 6;
>             isScrew = (chain % 4) == 0;
> 
>             normalizedBurg[X] = linedir[burgIndex][X];
>             normalizedBurg[Y] = linedir[burgIndex][Y];
>             normalizedBurg[Z] = linedir[burgIndex][Z];
> 
>             Normalize(&normalizedBurg[X], &normalizedBurg[Y],
>                       &normalizedBurg[Z]);
> 
> /*
>  *          Select an initial point (p0) that is within the boundaries
>  *          of the simulation and then calculate the positions at which
>  *          a dislocation line with the given line direction would
>  *          interesect the nearest free surface in each direction.
>  */
>             for (i = 0; i < 3; i++) {
>                 p0[i] = (randm(&seed)-0.5) * (0.5 * range[i]) + rangeCntr[i];
>             }
> 
>             minLen1 = 1.0e+20;
>             minLen2 = 1.0e+20;
> 
>             for (i = 0; i < 3; i++) {
>                 if (pbc[i] == 0) {
>                     signFact = (linedir[burgIndex][i] < 0.0 ? -1 : 1);
>                     surfCoord = signFact > 0 ? maxCoord[i] : minCoord[i];
>                     len = fabs((surfCoord - p0[i]) / normalizedBurg[i]);
>                     if (len < minLen1) {
>                         minLen1 = len;
>                         intersectIndex1 = i;
>                         intersectSurfCoord1 = surfCoord;
>                     }
> 
>                     signFact = -signFact;
>                     surfCoord = signFact > 0 ? maxCoord[i] : minCoord[i];
>                     len = fabs((surfCoord - p0[i]) / normalizedBurg[i]);
>                     if (len < minLen2) {
>                         minLen2 = len;
>                         intersectIndex2 = i;
>                         intersectSurfCoord2 = surfCoord;
>                     }
>                 }
>             }
> 
> /*
>  *          We know how far the dislocation can extend in each direction, now
>  *          calculate the exact intersect point in both directions
>  */
>             for (i = 0; i < 3; i++) {
> 
>                 if (i == intersectIndex1) {
>                     intersectPos1[i] = intersectSurfCoord1;
>                 } else {
>                     intersectPos1[i] = p0[i] + (minLen1 * normalizedBurg[i]);
>                 }
> 
>                 if (i == intersectIndex2) {
>                     intersectPos2[i] = intersectSurfCoord2;
>                 } else {
>                     intersectPos2[i] = p0[i] - (minLen2 * normalizedBurg[i]);
>                 }
>             }
> 
> /*
>  *          Find a vector from the first intersection point to the second,
>  *          calculate how many segments the line should be broken into based
>  *          on the <maxSeg> value.
>  */
>             for (i = 0; i < 3; i++) {
>                 vector[i] = intersectPos2[i] - intersectPos1[i];
>             }
> 
>             vecLen = sqrt(vector[0]*vector[0] +
>                           vector[1]*vector[1] +
>                           vector[2]*vector[2]);
> 
>             numSegs = (int)(vecLen / (.95 * param->maxSeg)) + 1;
>             numPoints = numSegs + 1;
> 
>             for (i = 0; i < 3; i++) {
>                 vector[i] /= (real8)numSegs;
>             }
> 
> /*
>  *          Reallocate the node array with sufficient size to add
>  *          all the new nodes defining this chain.
>  */
>             newNodeIndex = inData->nodeCount;
>             inData->nodeCount += numPoints;
>             inData->node = (Node_t *)realloc(inData->node, inData->nodeCount
>                                              * sizeof(Node_t));
>             memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numPoints);
>         
> 
> /*
>  *          Starting with the first intersection point, create a
>  *          series of dislocation segments ending at the second point.
>  */
>             newPos[X] = intersectPos1[X];
>             newPos[Y] = intersectPos1[Y];
>             newPos[Z] = intersectPos1[Z];
> 
>             FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);
> 
>             for (i = 0; i < numPoints; i++) {
>                 if (i == 0) {
>                     numConnections = 1;
>                     burgSign = -1.0;
>                     nbr1Index = 1;
>                 } else if (i == (numPoints - 1)) {
>                     numConnections = 1;
>                     burgSign = 1.0;
>                     nbr1Index = i - 1;
> /*
>  *                  Make sure the final node is at the surface
>  *                  intersection point
>  */
>                     newPos[X] = intersectPos2[X];
>                     newPos[Y] = intersectPos2[Y];
>                     newPos[Z] = intersectPos2[Z];
>                     FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);
>                 } else {
>                     numConnections = 2;
>                     burgSign = 1.0;
>                     nbr1Index = i - 1;
>                     nbr2Index = i + 1;
>                 }
> 
>                 node = &inData->node[baseNodeID+i];
>                 node->myTag.domainID = dislocType;
>                 node->myTag.index    = baseNodeID+i;
> 
>                 node->x = newPos[X];
>                 node->y = newPos[Y];
>                 node->z = newPos[Z];
> 
>                 node->constraint = UNCONSTRAINED;
> 
>                 AllocNodeArms(node, numConnections);
> 
>                 node->nbrTag[0].domainID = dislocType;
>                 node->nbrTag[0].index    = baseNodeID + nbr1Index;
> 
>                 node->burgX[0] = burgSign * linedir[burgIndex][X];
>                 node->burgY[0] = burgSign * linedir[burgIndex][Y];
>                 node->burgZ[0] = burgSign * linedir[burgIndex][Z];
> 
> /*
>  *              For screw dislocations, use the glide plane from the table.
>  *              For edge dislocations, glide plane is the cross product
>  *              of the burgers vector and line direction.
>  */
>                 if (isScrew) {
>                     node->nx[0] = glidePlane[gpBurgIndex][gpIndex][X];
>                     node->ny[0] = glidePlane[gpBurgIndex][gpIndex][Y];
>                     node->nz[0] = glidePlane[gpBurgIndex][gpIndex][Z];
>                     if (numConnections == 2) {
>                         node->nx[1] = glidePlane[gpBurgIndex][gpIndex][X];
>                         node->ny[1] = glidePlane[gpBurgIndex][gpIndex][Y];
>                         node->nz[1] = glidePlane[gpBurgIndex][gpIndex][Z];
>                     }
>                 } else {
>                     cross(linedir[burgIndex], linedir[lineIndex], plane);
>                     plane[X] = (floor(plane[X] * 1.0e+07)) * 1.0e-07;
>                     plane[Y] = (floor(plane[Y] * 1.0e+07)) * 1.0e-07;
>                     plane[Z] = (floor(plane[Z] * 1.0e+07)) * 1.0e-07;
>                     node->nx[0] = plane[X];
>                     node->ny[0] = plane[Y];
>                     node->nz[0] = plane[Z];
>                     if (numConnections == 2) {
>                         node->nx[1] = plane[X];
>                         node->ny[1] = plane[Y];
>                         node->nz[1] = plane[Z];
>                     }
>                 }
> 
> /*
>  *              Calculate the next node's position relative to this one.
>  */
>                 newPos[X] += vector[X];
>                 newPos[Y] += vector[Y];
>                 newPos[Z] += vector[Z];
> 
>                 FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);
> 
>                 if (numConnections == 1) {
>                     node->constraint = SURFACE_NODE;
>                     continue;
>                 }
> 
>                 node->nbrTag[1].domainID = dislocType;
>                 node->nbrTag[1].index    = baseNodeID + nbr2Index;
> 
>                 node->burgX[1] = -burgSign * linedir[burgIndex][X];
>                 node->burgY[1] = -burgSign * linedir[burgIndex][Y];
>                 node->burgZ[1] = -burgSign * linedir[burgIndex][Z];
>             }
> 
> /*
>  *          The initial segments created are not necessarily limited to
>  *          param->maxSegLen, so a call to InitRemesh() is needed to
>  *          chop any excessively long segments into proper lengths.
>  *          When we've generated the nodal data for the final chain,
>  *          write the block of nodal data to the file.
>  */
>             InitRemesh(inData, dislocType, startRemeshIndex);
> 
>             lastBlock = (chain == (numChains - 1));
> 
>             if (lastBlock) {
>                     IncDislocationDensity(inData, totDislocLen);
>                     param->nodeCount = inData->nodeCount;
>                     WriteInitialNodeData(home, inData, lastBlock,nodezero);
>                     FreeInNodeArray(inData, inData->nodeCount);
>                     inData->nodeCount = 0;
>             }
> 
>             baseNodeID = inData->nodeCount;
>             startRemeshIndex = baseNodeID;
> 
>         }  /* for (chain = 0; chain < numChains; ...) */
> 
>         return;
> }
> 
> 
> /*---------------------------------------------------------------------------
>  *
>  *      Function:     CreateFRSource
>  *      Description:  Generate a configuration consisting of one or more
>  *                    frank-read sources of random lengths between the
>  *                    specified minimum and maximum lengths.  
>  *
>  *                    NOTE:  This function currently assumes periodic
>  *                    boundary conditions!
>  *
>  *      Arguments:
>  *          cubeLength    Length of cubic problem space in a single
>  *                        dimension (units of b)
>  *          numSources    Number of frank read sources to create
>  *          srcLenMin     Minimum length of any created frank-read source
>  *          srcLenMax     Maximum length of any created frank-read source
>  *          seed          Seed value for random number generator
>  *          totDislocLen  Pointer to location at which to return
>  *                        to caller the total combined length of
>  *                        all created dislocations.
>  *
>  *-------------------------------------------------------------------------*/
> void CreateFRSource(Home_t *home, InData_t *inData, int cubeLength,
>                     int numSources, int srcLenMin, int srcLenMax, int seed,
>                     real8 *totDislocLen, int dislocType,int nodezero)
> {
>         int      i, lineIndex, chain, baseNodeID, numPoints;
>         int      newNodeIndex, lastBlock, burgIndex;
>         int      isScrew, gpIndex, gpBurgIndex;
>         int      lenRange;
>         int      startRemeshIndex = 0;
>         real8    srcLen;
>         real8    p0[3], p1[3], p2[3];
>         real8    linedir[4][3], glidePlane[4][6][3], plane[3];
>         Node_t   *node;
>         Param_t  *param;
> 
>         if (numSources <= 0) {
>             Fatal("%s: numSources is %d, but must be > 0.\n",
>                   "CreateFRSource", numSources);
>         }
> 
>         lenRange = srcLenMax - srcLenMin;
>         param = inData->param;
> 
> 
> /*
>  *      Set up an array of dislocation line directions that
>  *      are used to create the new dislocations.  There are
>  *      essentially 4 sets of line directions (1 per burgers
>  *      vector) with 4 line directions per set; the first for
>  *      the screw and the following three for edge.
>  *      The burger's vector for all 4 dislocation lines in
>  *      a set is the same as the line direction of the screw
>  *      dislocation in the group.
>  */
> 
> /*
>  *      Type  [1 1 1] burgers vector
>  */
>         linedir[0][0] =  0.5;
>         linedir[0][1] =  0.5;
>         linedir[0][2] =  0.5;
> 
> /*
>  *      Type [-1 1 1] burgers vector
>  */
>         linedir[1][0] = -0.5;
>         linedir[1][1] =  0.5;
>         linedir[1][2] =  0.5;
> 
> /*
>  *      Type [1 -1 1] burgers vector
>  */
>         linedir[2][0] =  0.5;
>         linedir[2][1] = -0.5;
>         linedir[2][2] =  0.5;
> 
> /*
>  *      Type [1 1 -1] burgers vector
>  */
>         linedir[3][0] =  0.5;
>         linedir[3][1] =  0.5;
>         linedir[3][2] = -0.5;
> 
> /*
>  *      Set up the valid glide planes for each screw burgers vector,
>  *      six glide planes per burgers vector.  For edges, glide plane
>  *      will simply be cross product between burgers vector and
>  *      linedir.
>  *
>  *      glide planes for [1 1 1]
>  */
>         glidePlane[0][0][0] =  0.7071068;
>         glidePlane[0][0][1] = -0.7071068;
>         glidePlane[0][0][2] =  0.0000000;
> 
>         glidePlane[0][1][0] =  0.7071068;
>         glidePlane[0][1][1] =  0.0000000;
>         glidePlane[0][1][2] = -0.7071068;
> 
>         glidePlane[0][2][0] =  0.0000000;
>         glidePlane[0][2][1] =  0.7071068;
>         glidePlane[0][2][2] = -0.7071068;
> 
>         glidePlane[0][3][0] = -0.7071068;
>         glidePlane[0][3][1] =  0.7071068;
>         glidePlane[0][3][2] = -0.0000000;
> 
>         glidePlane[0][4][0] = -0.7071068;
>         glidePlane[0][4][1] = -0.0000000;
>         glidePlane[0][4][2] =  0.7071068;
> 
>         glidePlane[0][5][0] = -0.0000000;
>         glidePlane[0][5][1] = -0.7071068;
>         glidePlane[0][5][2] =  0.7071068;
> 
> /*
>  *      glide planes for [-1 1 1]
>  */
>         glidePlane[1][0][0] =  0.0000000;
>         glidePlane[1][0][1] =  0.7071068;
>         glidePlane[1][0][2] = -0.7071068;
> 
>         glidePlane[1][1][0] =  0.7071068;
>         glidePlane[1][1][1] =  0.0000000;
>         glidePlane[1][1][2] =  0.7071068;
> 
>         glidePlane[1][2][0] =  0.0000000;
1896,1907d2982
<         minCoord[X] = param->xBoundMin;
<         maxCoord[X] = param->xBoundMax;
<         minCoord[Y] = param->yBoundMin;
<         maxCoord[Y] = param->yBoundMax;
<         minCoord[Z] = param->zBoundMin;
<         maxCoord[Z] = param->zBoundMax;
< 
<         for (i = 0; i < 3; i++) {
<             range[i] = maxCoord[i] - minCoord[i];
<             rangeCntr[i] = 0.5 * (minCoord[i] + maxCoord[i]);
<         }
< 
1916,1927c2991
<         for (chain = 0; chain < numChains; chain++) {
<         
<             numPoints = 3;
<             lineIndex = chain % 16; 
<             burgIndex = 4 * (lineIndex / 4);
<             gpBurgIndex = lineIndex / 4;
<             gpIndex = (chain / 16) % 6;
<             isScrew = (chain % 4) == 0;
< 
<             normalizedBurg[X] = linedir[burgIndex][X];
<             normalizedBurg[Y] = linedir[burgIndex][Y];
<             normalizedBurg[Z] = linedir[burgIndex][Z];
---
>         for (chain = 0; chain < numSources; chain++) {
1929,1930c2993,2997
<             Normalize(&normalizedBurg[X], &normalizedBurg[Y],
<                       &normalizedBurg[Z]);
---
>             numPoints = 3;
>             lineIndex = chain % 4;
>             burgIndex = lineIndex;
>             gpBurgIndex = lineIndex;
>             gpIndex = (chain / 4) % 6;
1933,1936c3000,3001
<  *          Select an initial point (p0) that is within the boundaries
<  *          of the simulation and then calculate the positions at which
<  *          a dislocation line with the given line direction would
<  *          interesect the nearest free surface in each direction.
---
>  *          Reallocate the node array with sufficient size to add
>  *          all the new nodes defining this chain.
1938,1965c3003,3007
<             for (i = 0; i < 3; i++) {
<                 p0[i] = (randm(&seed)-0.5) * (0.5 * range[i]) + rangeCntr[i];
<             }
< 
<             minLen1 = 1.0e+20;
<             minLen2 = 1.0e+20;
< 
<             for (i = 0; i < 3; i++) {
<                 if (pbc[i] == 0) {
<                     signFact = (linedir[burgIndex][i] < 0.0 ? -1 : 1);
<                     surfCoord = signFact > 0 ? maxCoord[i] : minCoord[i];
<                     len = fabs((surfCoord - p0[i]) / normalizedBurg[i]);
<                     if (len < minLen1) {
<                         minLen1 = len;
<                         intersectIndex1 = i;
<                         intersectSurfCoord1 = surfCoord;
<                     }
< 
<                     signFact = -signFact;
<                     surfCoord = signFact > 0 ? maxCoord[i] : minCoord[i];
<                     len = fabs((surfCoord - p0[i]) / normalizedBurg[i]);
<                     if (len < minLen2) {
<                         minLen2 = len;
<                         intersectIndex2 = i;
<                         intersectSurfCoord2 = surfCoord;
<                     }
<                 }
<             }
---
>             newNodeIndex = inData->nodeCount;
>             inData->nodeCount += numPoints;
>             inData->node = (Node_t *)realloc(inData->node, inData->nodeCount
>                                                * sizeof(Node_t));
>             memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numPoints);
1968,1969c3010,3011
<  *          We know how far the dislocation can extend in each direction, now
<  *          calculate the exact intersect point in both directions
---
>  *          Length of the frank read source should be a random length
>  *          between srcLenMin and srcLenMax.
1971,1977c3013,3017
<             for (i = 0; i < 3; i++) {
< 
<                 if (i == intersectIndex1) {
<                     intersectPos1[i] = intersectSurfCoord1;
<                 } else {
<                     intersectPos1[i] = p0[i] + (minLen1 * normalizedBurg[i]);
<                 }
---
>             if (lenRange > 0) {
>                 srcLen = (srcLenMin + (randm(&seed) * lenRange)) / sqrt(3.0);
>             } else {
>                 srcLen = srcLenMin / sqrt(3.0);
>             }
1979,1983c3019,3022
<                 if (i == intersectIndex2) {
<                     intersectPos2[i] = intersectSurfCoord2;
<                 } else {
<                     intersectPos2[i] = p0[i] - (minLen2 * normalizedBurg[i]);
<                 }
---
>             for (i = 0; i < 3; i++) {
>                 p0[i] = (randm(&seed)-0.5) * cubeLength * 2.0;
>                 p1[i] = p0[i] - (srcLen * linedir[lineIndex][i]);
>                 p2[i] = p0[i] + (srcLen * linedir[lineIndex][i]);
1987,1989c3026,3027
<  *          Find a vector from the first intersection point to the second,
<  *          calculate how many segments the line should be broken into based
<  *          on the <maxSeg> value.
---
>  *          Set up the 3 nodes we're using to define the dislocation line.
>  *          First do point p0.
1991,1993c3029
<             for (i = 0; i < 3; i++) {
<                 vector[i] = intersectPos2[i] - intersectPos1[i];
<             }
---
>             node = &inData->node[baseNodeID];
1995,1997c3031,3035
<             vecLen = sqrt(vector[0]*vector[0] +
<                           vector[1]*vector[1] +
<                           vector[2]*vector[2]);
---
>             node->myTag.domainID = dislocType;
>             node->myTag.index    = baseNodeID;
>             node->x = p0[0];
>             node->y = p0[1];
>             node->z = p0[2];
1999,2000c3037
<             numSegs = (int)(vecLen / (.95 * param->maxSeg)) + 1;
<             numPoints = numSegs + 1;
---
>             node->constraint = UNCONSTRAINED;
2002,2004c3039
<             for (i = 0; i < 3; i++) {
<                 vector[i] /= (real8)numSegs;
<             }
---
>             AllocNodeArms(node, 2);
2006,2015c3041,3052
< /*
<  *          Reallocate the node array with sufficient size to add
<  *          all the new nodes defining this chain.
<  */
<             newNodeIndex = inData->nodeCount;
<             inData->nodeCount += numPoints;
<             inData->node = (Node_t *)realloc(inData->node, inData->nodeCount
<                                              * sizeof(Node_t));
<             memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numPoints);
<         
---
>             node->nbrTag[0].domainID = dislocType;
>             node->nbrTag[0].index    = baseNodeID + 1;
>             node->nbrTag[1].domainID = dislocType;
>             node->nbrTag[1].index    = baseNodeID + 2;
> 
>             node->burgX[0] = linedir[burgIndex][0];
>             node->burgY[0] = linedir[burgIndex][1];
>             node->burgZ[0] = linedir[burgIndex][2];
> 
>             node->burgX[1] = -linedir[burgIndex][0];
>             node->burgY[1] = -linedir[burgIndex][1];
>             node->burgZ[1] = -linedir[burgIndex][2];
2018,2019c3055,3057
<  *          Starting with the first intersection point, create a
<  *          series of dislocation segments ending at the second point.
---
>  *          For screw dislocations, use the glide plane from the table.
>  *          For edge dislocations, glide plane is the cross product
>  *          of the burgers vector and line direction.
2021,2025c3059,3077
<             newPos[X] = intersectPos1[X];
<             newPos[Y] = intersectPos1[Y];
<             newPos[Z] = intersectPos1[Z];
< 
<             FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);
---
>             if (isScrew) {
>                 node->nx[0] = glidePlane[gpBurgIndex][gpIndex][0];
>                 node->ny[0] = glidePlane[gpBurgIndex][gpIndex][1];
>                 node->nz[0] = glidePlane[gpBurgIndex][gpIndex][2];
>                 node->nx[1] = glidePlane[gpBurgIndex][gpIndex][0];
>                 node->ny[1] = glidePlane[gpBurgIndex][gpIndex][1];
>                 node->nz[1] = glidePlane[gpBurgIndex][gpIndex][2];
>             } else {
>                 cross(linedir[burgIndex], linedir[lineIndex], plane);
>                 plane[0] = (floor(plane[0] * 1.0e+07)) * 1.0e-07;
>                 plane[1] = (floor(plane[1] * 1.0e+07)) * 1.0e-07;
>                 plane[2] = (floor(plane[2] * 1.0e+07)) * 1.0e-07;
>                 node->nx[0] = plane[0];
>                 node->ny[0] = plane[1];
>                 node->nz[0] = plane[2];
>                 node->nx[1] = plane[0];
>                 node->ny[1] = plane[1];
>                 node->nz[1] = plane[2];
>             }
2027,2035d3078
<             for (i = 0; i < numPoints; i++) {
<                 if (i == 0) {
<                     numConnections = 1;
<                     burgSign = -1.0;
<                     nbr1Index = 1;
<                 } else if (i == (numPoints - 1)) {
<                     numConnections = 1;
<                     burgSign = 1.0;
<                     nbr1Index = i - 1;
2037,2038c3080
<  *                  Make sure the final node is at the surface
<  *                  intersection point
---
>  *          Now point p1...
2040,2049c3082
<                     newPos[X] = intersectPos2[X];
<                     newPos[Y] = intersectPos2[Y];
<                     newPos[Z] = intersectPos2[Z];
<                     FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);
<                 } else {
<                     numConnections = 2;
<                     burgSign = 1.0;
<                     nbr1Index = i - 1;
<                     nbr2Index = i + 1;
<                 }
---
>             node = &inData->node[baseNodeID+1];
2051,2053c3084,3085
<                 node = &inData->node[baseNodeID+i];
<                 node->myTag.domainID = dislocType;
<                 node->myTag.index    = baseNodeID+i;
---
>             node->myTag.domainID = dislocType;
>             node->myTag.index    = baseNodeID+1;
2055,2057c3087,3089
<                 node->x = newPos[X];
<                 node->y = newPos[Y];
<                 node->z = newPos[Z];
---
>             node->x = p1[0];
>             node->y = p1[1];
>             node->z = p1[2];
2059c3091
<                 node->constraint = UNCONSTRAINED;
---
>             node->constraint = PINNED_NODE;
2061c3093
<                 AllocNodeArms(node, numConnections);
---
>             AllocNodeArms(node, 1);
2063,2064c3095,3096
<                 node->nbrTag[0].domainID = dislocType;
<                 node->nbrTag[0].index    = baseNodeID + nbr1Index;
---
>             node->nbrTag[0].domainID = dislocType;
>             node->nbrTag[0].index    = baseNodeID;
2066,2068c3098,3103
<                 node->burgX[0] = burgSign * linedir[burgIndex][X];
<                 node->burgY[0] = burgSign * linedir[burgIndex][Y];
<                 node->burgZ[0] = burgSign * linedir[burgIndex][Z];
---
>             node->burgX[0] = -linedir[burgIndex][0];
>             node->burgY[0] = -linedir[burgIndex][1];
>             node->burgZ[0] = -linedir[burgIndex][2];
>             node->nx[0] = glidePlane[gpBurgIndex][gpIndex][0];
>             node->ny[0] = glidePlane[gpBurgIndex][gpIndex][1];
>             node->nz[0] = glidePlane[gpBurgIndex][gpIndex][2];
2071,2073c3106
<  *              For screw dislocations, use the glide plane from the table.
<  *              For edge dislocations, glide plane is the cross product
<  *              of the burgers vector and line direction.
---
>  *          Now point p2...
2075,2097c3108
<                 if (isScrew) {
<                     node->nx[0] = glidePlane[gpBurgIndex][gpIndex][X];
<                     node->ny[0] = glidePlane[gpBurgIndex][gpIndex][Y];
<                     node->nz[0] = glidePlane[gpBurgIndex][gpIndex][Z];
<                     if (numConnections == 2) {
<                         node->nx[1] = glidePlane[gpBurgIndex][gpIndex][X];
<                         node->ny[1] = glidePlane[gpBurgIndex][gpIndex][Y];
<                         node->nz[1] = glidePlane[gpBurgIndex][gpIndex][Z];
<                     }
<                 } else {
<                     cross(linedir[burgIndex], linedir[lineIndex], plane);
<                     plane[X] = (floor(plane[X] * 1.0e+07)) * 1.0e-07;
<                     plane[Y] = (floor(plane[Y] * 1.0e+07)) * 1.0e-07;
<                     plane[Z] = (floor(plane[Z] * 1.0e+07)) * 1.0e-07;
<                     node->nx[0] = plane[X];
<                     node->ny[0] = plane[Y];
<                     node->nz[0] = plane[Z];
<                     if (numConnections == 2) {
<                         node->nx[1] = plane[X];
<                         node->ny[1] = plane[Y];
<                         node->nz[1] = plane[Z];
<                     }
<                 }
---
>             node = &inData->node[baseNodeID+2];
2099,2104c3110,3111
< /*
<  *              Calculate the next node's position relative to this one.
<  */
<                 newPos[X] += vector[X];
<                 newPos[Y] += vector[Y];
<                 newPos[Z] += vector[Z];
---
>             node->myTag.domainID = dislocType;
>             node->myTag.index    = baseNodeID+2;
2106c3113,3115
<                 FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);
---
>             node->x = p2[0];
>             node->y = p2[1];
>             node->z = p2[2];
2108,2111c3117
<                 if (numConnections == 1) {
<                     node->constraint = SURFACE_NODE;
<                     continue;
<                 }
---
>             node->constraint = PINNED_NODE;
2113,2114c3119
<                 node->nbrTag[1].domainID = dislocType;
<                 node->nbrTag[1].index    = baseNodeID + nbr2Index;
---
>             AllocNodeArms(node, 1);
2116,2119c3121,3129
<                 node->burgX[1] = -burgSign * linedir[burgIndex][X];
<                 node->burgY[1] = -burgSign * linedir[burgIndex][Y];
<                 node->burgZ[1] = -burgSign * linedir[burgIndex][Z];
<             }
---
>             node->nbrTag[0].domainID = dislocType;
>             node->nbrTag[0].index    = baseNodeID;
> 
>             node->burgX[0] = linedir[burgIndex][0];
>             node->burgY[0] = linedir[burgIndex][1];
>             node->burgZ[0] = linedir[burgIndex][2];
>             node->nx[0] = glidePlane[gpBurgIndex][gpIndex][0];
>             node->ny[0] = glidePlane[gpBurgIndex][gpIndex][1];
>             node->nz[0] = glidePlane[gpBurgIndex][gpIndex][2];
2130c3140
<             lastBlock = (chain == (numChains - 1));
---
>             lastBlock = (chain == (numSources - 1));
2135c3145
<                     WriteInitialNodeData(home, inData, lastBlock);
---
>                     WriteInitialNodeData(home, inData, lastBlock,nodezero);
2151,2157c3161,3169
<  *      Function:     CreateFRSource
<  *      Description:  Generate a configuration consisting of one or more
<  *                    frank-read sources of random lengths between the
<  *                    specified minimum and maximum lengths.  
<  *
<  *                    NOTE:  This function currently assumes periodic
<  *                    boundary conditions!
---
>  *      Function:     CreateEdges
>  *      Description:  Creates "edge" dislocations in the [100] [010] and [001]
>  *                    directions.  Each of the three line senses has 4
>  *                    combinations of burgers vector and normals.  With the
>  *                    opposite sign of the line sense vectors considered,
>  *                    the total number of types of dislocations is 24,
>  *                    so number of chains specified should be a multiple
>  *                    of 24 to include all types of.  (Okay, they're
>  *                    not pure edge, but they are not screw either)
2160,2168c3172,3175
<  *          cubeLength    Length of cubic problem space in a single
<  *                        dimension (units of b)
<  *          numSources    Number of frank read sources to create
<  *          srcLenMin     Minimum length of any created frank-read source
<  *          srcLenMax     Maximum length of any created frank-read source
<  *          seed          Seed value for random number generator
<  *          totDislocLen  Pointer to location at which to return
<  *                        to caller the total combined length of
<  *                        all created dislocations.
---
>  *          cubeLength  Length of cubic problem space in a single
>  *                      dimension (units of b)
>  *          numChains   Number of chains to create
>  *          seed        Seed value for random number generator
2171,2173c3178,3180
< void CreateFRSource(Home_t *home, InData_t *inData, int cubeLength,
<                     int numSources, int srcLenMin, int srcLenMax, int seed,
<                     real8 *totDislocLen, int dislocType)
---
> void CreateEdges(Home_t *home, InData_t *inData, int cubeLength,
>                  int numChains, int seed, real8 *totDislocLen,
>                  int dislocType,int nodezero)
2175,2184c3182,3210
<         int      i, lineIndex, chain, baseNodeID, numPoints;
<         int      newNodeIndex, lastBlock, burgIndex;
<         int      isScrew, gpIndex, gpBurgIndex;
<         int      lenRange;
<         int      startRemeshIndex = 0;
<         real8    srcLen;
<         real8    p0[3], p1[3], p2[3];
<         real8    linedir[4][3], glidePlane[4][6][3], plane[3];
<         Node_t   *node;
<         Param_t  *param;
---
>         int           ic, ip, np, id0, lastBlock, signFactor;
>         int           newNodeIndex, startRemeshIndex;
>         int           ldIndex, gpIndex, burgIndex, nbr1Index, nbr2Index;
>         real8         posFactor, cubeSize;
>         real8         xp[3], yp[3], zp[3];
>         Param_t       *param;
>         Node_t        *node;
>         static real8  lineDir[3][3] = {
>                       {0.0, 0.0, 1.0},
>                       {0.0, 1.0, 0.0},
>                       {1.0, 0.0, 0.0}};
>         static real8  burg[4][3] = {
>                       { 0.5773503,  0.5773503,  0.5773503},
>                       { 0.5773503,  0.5773503, -0.5773503},
>                       { 0.5773503, -0.5773503,  0.5773503},
>                       {-0.5773503,  0.5773503,  0.5773503}};
>         static real8  glidePlane[12][3] = {
>                       {  0.7071068, -0.7071068,  0},  /* ldir [001] b [111]  */
>                       {  0.7071068, -0.7071068,  0},  /* ldir [001] b [11-1] */
>                       {  0.7071068,  0.7071068,  0},  /* ldir [001] b [1-11] */
>                       {  0.7071068,  0.7071068,  0},  /* ldir [001] b [-111] */
>                       {  0.7071068,  0, -0.7071068},  /* ldir [010] b [111]  */
>                       {  0.7071068,  0,  0.7071068},  /* ldir [010] b [11-1] */
>                       {  0.7071068,  0, -0.7071068},  /* ldir [010] b [1-11] */
>                       {  0.7071068,  0,  0.7071068},  /* ldir [010] b [-111] */
>                       {  0, -0.7071068,  0.7071068},  /* ldir [100] b [111]  */
>                       {  0,  0.7071068,  0.7071068},  /* ldir [100] b [11-1] */
>                       {  0,  0.7071068,  0.7071068},  /* ldir [100] b [1-11] */
>                       {  0, -0.7071068,  0.7071068}}; /* ldir [100] b [-111] */
2186,2188c3212,3215
<         if (numSources <= 0) {
<             Fatal("%s: numSources is %d, but must be > 0.\n",
<                   "CreateFRSource", numSources);
---
> 
>         if (numChains <= 0) {
>             Fatal("%s: numChains is %d, but must be > 0.\n",
>                   "CreateEdges", numChains);
2191d3217
<         lenRange = srcLenMax - srcLenMin;
2193,2211c3219,3223
< 
< 
< /*
<  *      Set up an array of dislocation line directions that
<  *      are used to create the new dislocations.  There are
<  *      essentially 4 sets of line directions (1 per burgers
<  *      vector) with 4 line directions per set; the first for
<  *      the screw and the following three for edge.
<  *      The burger's vector for all 4 dislocation lines in
<  *      a set is the same as the line direction of the screw
<  *      dislocation in the group.
<  */
< 
< /*
<  *      Type  [1 1 1] burgers vector
<  */
<         linedir[0][0] =  0.5;
<         linedir[0][1] =  0.5;
<         linedir[0][2] =  0.5;
---
>         cubeSize = (real8)cubeLength;
>         id0 = 0;
>         inData->nodeCount = 0;
>         startRemeshIndex = 0;
>         posFactor = 0.333 * cubeSize;
2214c3226
<  *      Type [-1 1 1] burgers vector
---
>  *      Create the specified number of chains.
2216,2218c3228
<         linedir[1][0] = -0.5;
<         linedir[1][1] =  0.5;
<         linedir[1][2] =  0.5;
---
>         for (ic = 0; ic < numChains; ic++) {
2219a3230,3233
>             gpIndex = ic % 12; 
>             burgIndex = gpIndex % 4;
>             ldIndex = gpIndex / 4;
>             
2221c3235,3237
<  *      Type [1 -1 1] burgers vector
---
>  *          First 12 burgers vector/normal sets use positive line
>  *          sense, next set of 12 uses opposite line sense, and
>  *          so on.
2223,2225c3239,3247
<         linedir[2][0] =  0.5;
<         linedir[2][1] = -0.5;
<         linedir[2][2] =  0.5;
---
>             signFactor = ((ic / 12) & 0x01) ? -1 : 1;
> 
>             np = 3;
>             newNodeIndex = inData->nodeCount;
>             inData->nodeCount += np;
> 
>             inData->node = (Node_t *)realloc(inData->node,
>                            inData->nodeCount * sizeof(Node_t));
>             memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);
2228c3250,3252
<  *      Type [1 1 -1] burgers vector
---
>  *          Set up 3 initial points for the line.  Point 1 is a base position
>  *          at a random location, point 0 is in the negative direction along
>  *          the line and point 2 is in the positive direction along the line.
2230,2232c3254,3264
<         linedir[3][0] =  0.5;
<         linedir[3][1] =  0.5;
<         linedir[3][2] = -0.5;
---
>             xp[1] = (randm(&seed)-0.5)*cubeSize;
>             yp[1] = (randm(&seed)-0.5)*cubeSize;
>             zp[1] = (randm(&seed)-0.5)*cubeSize;
> 
>             xp[0] = xp[1] - (posFactor * signFactor * lineDir[ldIndex][X]);
>             yp[0] = yp[1] - (posFactor * signFactor * lineDir[ldIndex][Y]);
>             zp[0] = zp[1] - (posFactor * signFactor * lineDir[ldIndex][Z]);
> 
>             xp[2] = xp[1] + (posFactor * signFactor * lineDir[ldIndex][X]);
>             yp[2] = yp[1] + (posFactor * signFactor * lineDir[ldIndex][Y]);
>             zp[2] = zp[1] + (posFactor * signFactor * lineDir[ldIndex][Z]);
2235,2240c3267,3268
<  *      Set up the valid glide planes for each screw burgers vector,
<  *      six glide planes per burgers vector.  For edges, glide plane
<  *      will simply be cross product between burgers vector and
<  *      linedir.
<  *
<  *      glide planes for [1 1 1]
---
>  *          Loop over the points and set up the nodes, link them to 
>  *          the neighbor nodes, etc.
2242,2244c3270
<         glidePlane[0][0][0] =  0.7071068;
<         glidePlane[0][0][1] = -0.7071068;
<         glidePlane[0][0][2] =  0.0000000;
---
>             for (ip = 0; ip < np; ip++) {
2246,2248c3272
<         glidePlane[0][1][0] =  0.7071068;
<         glidePlane[0][1][1] =  0.0000000;
<         glidePlane[0][1][2] = -0.7071068;
---
>                 node = &inData->node[ip+id0];
2250,2252c3274,3276
<         glidePlane[0][2][0] =  0.0000000;
<         glidePlane[0][2][1] =  0.7071068;
<         glidePlane[0][2][2] = -0.7071068;
---
>                 node->x = xp[ip];
>                 node->y = yp[ip];
>                 node->z = zp[ip];
2254,2256c3278,3280
<         glidePlane[0][3][0] = -0.7071068;
<         glidePlane[0][3][1] =  0.7071068;
<         glidePlane[0][3][2] = -0.0000000;
---
>                 node->constraint = UNCONSTRAINED;
>                 node->myTag.domainID = dislocType;
>                 node->myTag.index = ip+id0;
2258,2260c3282
<         glidePlane[0][4][0] = -0.7071068;
<         glidePlane[0][4][1] = -0.0000000;
<         glidePlane[0][4][2] =  0.7071068;
---
>                 AllocNodeArms(node, 2);
2262,2264c3284,3305
<         glidePlane[0][5][0] = -0.0000000;
<         glidePlane[0][5][1] = -0.7071068;
<         glidePlane[0][5][2] =  0.7071068;
---
>                 if ((nbr1Index = ip + 1) >= np) nbr1Index = 0;
>                 if ((nbr2Index= ip - 1) < 0) nbr2Index = np - 1;
> 
>                 node->nbrTag[0].domainID = dislocType;
>                 node->nbrTag[0].index = id0 + nbr1Index;
>                 node->burgX[0] = burg[burgIndex][0];
>                 node->burgY[0] = burg[burgIndex][1];
>                 node->burgZ[0] = burg[burgIndex][2];
>                 node->nx[0] = glidePlane[gpIndex][X];
>                 node->ny[0] = glidePlane[gpIndex][Y];
>                 node->nz[0] = glidePlane[gpIndex][Z];
>             
>                 node->nbrTag[1].domainID = dislocType;
>                 node->nbrTag[1].index = id0 + nbr2Index;
>                 node->burgX[1] = -burg[burgIndex][0];
>                 node->burgY[1] = -burg[burgIndex][1];
>                 node->burgZ[1] = -burg[burgIndex][2];
>                 node->nx[1] = glidePlane[gpIndex][X];
>                 node->ny[1] = glidePlane[gpIndex][Y];
>                 node->nz[1] = glidePlane[gpIndex][Z];
> 
>             }
2267c3308,3312
<  *      glide planes for [-1 1 1]
---
>  *          The initial segments created are not necessarily limited to
>  *          param->maxSegLen, so a call to InitRemesh() is needed to
>  *          chop any excessively long segments into proper lengths.
>  *          When we've generated the nodal data for the final chain,
>  *          write the current block of nodal data to the file.
2269,2271c3314
<         glidePlane[1][0][0] =  0.0000000;
<         glidePlane[1][0][1] =  0.7071068;
<         glidePlane[1][0][2] = -0.7071068;
---
>             InitRemesh(inData, dislocType, startRemeshIndex);
2273,2275c3316,3323
<         glidePlane[1][1][0] =  0.7071068;
<         glidePlane[1][1][1] =  0.0000000;
<         glidePlane[1][1][2] =  0.7071068;
---
>             lastBlock = (ic == numChains - 1);
>             if (lastBlock) {
>                 IncDislocationDensity(inData, totDislocLen);
>                 param->nodeCount = inData->nodeCount;
>                 WriteInitialNodeData(home, inData, lastBlock,nodezero);
>                 FreeInNodeArray(inData, inData->nodeCount);
>                 inData->nodeCount = 0;
>             }
2277,2279c3325,3329
<         glidePlane[1][2][0] =  0.0000000;
<         glidePlane[1][2][1] =  0.7071068;
<         glidePlane[1][2][2] = -0.7071068;
---
>             id0 = inData->nodeCount;
>             startRemeshIndex = id0;
>         }
>         return;
> }
2281,2283d3330
<         glidePlane[1][3][0] = -0.0000000;
<         glidePlane[1][3][1] = -0.7071068;
<         glidePlane[1][3][2] =  0.7071068;
2285,2287d3331
<         glidePlane[1][4][0] = -0.7071068;
<         glidePlane[1][4][1] = -0.0000000;
<         glidePlane[1][4][2] = -0.7071068;
2289,2291c3333,3397
<         glidePlane[1][5][0] = -0.0000000;
<         glidePlane[1][5][1] = -0.7071068;
<         glidePlane[1][5][2] =  0.7071068;
---
> /*---------------------------------------------------------------------------
>  *Modified to generate real edges! (AL) 
>  *      Function:     CreateEdges
>  *      Description:  Creates "edge" dislocations in the [100] [010] and [001]
>  *                    directions.  Each of the three line senses has 4
>  *                    combinations of burgers vector and normals.  With the
>  *                    opposite sign of the line sense vectors considered,
>  *                    the total number of types of dislocations is 24,
>  *                    so number of chains specified should be a multiple
>  *                    of 24 to include all types of.  (Okay, they're
>  *                    not pure edge, but they are not screw either)
>  *
>  *      Arguments:
>  *          cubeLength  Length of cubic problem space in a single
>  *                      dimension (units of b)
>  *          numChains   Number of chains to create
>  *          seed        Seed value for random number generator
>  *
>  *-------------------------------------------------------------------------*/
> //void CreateEdges(Home_t *home, InData_t *inData, int cubeLength,
>                  //int numChains, int seed, real8 *totDislocLen,
>                  //int dislocType,int nodezero)
> //{
>         //int           ic, ip, np, id0, lastBlock, signFactor;
>         //int           newNodeIndex, startRemeshIndex;
>         //int           ldIndex, gpIndex, burgIndex, nbr1Index, nbr2Index;
>         //real8         posFactor, cubeSize;
>         //real8         xp[3], yp[3], zp[3];
>         //Param_t       *param;
>         //Node_t        *node;
>         //static real8  lineDir[3][3] = {
>                       //{0.0, 0.0, 1.0},
>                       //{0.0, 1.0, 0.0},
>                       //{1.0, 0.0, 0.0}};
>         //static real8  burg[3][3] = {
>                       //{ 0.0,  1.0,  0.0},
>                       //{ 1.0,  0.0, 0.0},
>                       //{ 0.0, 0.0,  1.0}};
>         //static real8  glidePlane[6][3] = {
>                       //{  1.0, 0.0,  0.0},  			  /* ldir [001] b [010]  */
>                       //{  0.0, 1.0,  0.0},  			  /* ldir [001] b [100] */
>                       //{  1.0,  0.0, 0.0},  			  /* ldir [010] b [001]  */
>                       //{  0.0,  0.0,  1.0},  		  /* ldir [010] b [100] */
>                       //{  0.0,  0.0, 1.0},  			  /* ldir [100] b [010] */
>                       //{  0, 1.0,  0.0}}; 			  /* ldir [100] b [001] */
> 
> 
>         //if (numChains <= 0) {
>             //Fatal("%s: numChains is %d, but must be > 0.\n",
>                   //"CreateEdges", numChains);
>         //}
> 
>         //param = inData->param;
>         //cubeSize = (real8)cubeLength;
>         //id0 = 0;
>         //inData->nodeCount = 0;
>         //startRemeshIndex = 0;
>         //posFactor = 0.333 * cubeSize;
> 
> //
>  //*      Create the specified number of chains.
>  //*/
>  
> 		////for (ic = 5; ic < numChains; ic++) { for2d
>         //for (ic = 0; ic < numChains; ic++) {
2293,2298c3399,3503
< /*
<  *      glide planes for [1 -1 1]
<  */
<         glidePlane[2][0][0] =  0.7071068;
<         glidePlane[2][0][1] =  0.7071068;
<         glidePlane[2][0][2] =  0.0000000;
---
>             //gpIndex = ic % 6; 
>             //burgIndex = gpIndex % 3;
>             //ldIndex = gpIndex / 3;
>             
> //
>  //*          First 12 burgers vector/normal sets use positive line
>  //*          sense, next set of 12 uses opposite line sense, and
>  //*          so on.
>  //*/
>             //signFactor = ((ic / 6) & 0x01) ? -1 : 1;
> 
>             //np = 3;
>             //newNodeIndex = inData->nodeCount;
>             //inData->nodeCount += np;
> 
>             //inData->node = (Node_t *)realloc(inData->node,
>                            //inData->nodeCount * sizeof(Node_t));
>             //memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);
> 
> //
>  //*          Set up 3 initial points for the line.  Point 1 is a base position
>  //*          at a random location, point 0 is in the negative direction along
>  //*          the line and point 2 is in the positive direction along the line.
>  //*/
>             //xp[1] = (randm(&seed)-0.5)*cubeSize;
>             //yp[1] = (randm(&seed)-0.5)*cubeSize;
>             //zp[1] = (randm(&seed)-0.5)*cubeSize;
>             
>             ////for single glide plane (AL)
>             ////xp[1] = -1000.0;
>             ////yp[1] = 0.000;
>             ////zp[1] = 1000.00;
> 
>             //xp[0] = xp[1] - (posFactor * signFactor * lineDir[ldIndex][X]);
>             //yp[0] = yp[1] - (posFactor * signFactor * lineDir[ldIndex][Y]);
>             //zp[0] = zp[1] - (posFactor * signFactor * lineDir[ldIndex][Z]);
> 
>             //xp[2] = xp[1] + (posFactor * signFactor * lineDir[ldIndex][X]);
>             //yp[2] = yp[1] + (posFactor * signFactor * lineDir[ldIndex][Y]);
>             //zp[2] = zp[1] + (posFactor * signFactor * lineDir[ldIndex][Z]);
> 
> //
>  //*          Loop over the points and set up the nodes, link them to 
>  //*          the neighbor nodes, etc.
>  //*/
>             //for (ip = 0; ip < np; ip++) {
> 
>                 //node = &inData->node[ip+id0];
> 
>                 //node->x = xp[ip];
>                 //node->y = yp[ip];
>                 //node->z = zp[ip];
> 
>                 //node->constraint = UNCONSTRAINED;
>                 //node->myTag.domainID = dislocType;
>                 //node->myTag.index = ip+id0;
> 
>                 //AllocNodeArms(node, 2);
> 
>                 //if ((nbr1Index = ip + 1) >= np) nbr1Index = 0;
>                 //if ((nbr2Index= ip - 1) < 0) nbr2Index = np - 1;
> 
>                 //node->nbrTag[0].domainID = dislocType;
>                 //node->nbrTag[0].index = id0 + nbr1Index;
>                 //node->burgX[0] = burg[burgIndex][0];
>                 //node->burgY[0] = burg[burgIndex][1];
>                 //node->burgZ[0] = burg[burgIndex][2];
>                 //node->nx[0] = glidePlane[gpIndex][X];
>                 //node->ny[0] = glidePlane[gpIndex][Y];
>                 //node->nz[0] = glidePlane[gpIndex][Z];
>             
>                 //node->nbrTag[1].domainID = dislocType;
>                 //node->nbrTag[1].index = id0 + nbr2Index;
>                 //node->burgX[1] = -burg[burgIndex][0];
>                 //node->burgY[1] = -burg[burgIndex][1];
>                 //node->burgZ[1] = -burg[burgIndex][2];
>                 //node->nx[1] = glidePlane[gpIndex][X];
>                 //node->ny[1] = glidePlane[gpIndex][Y];
>                 //node->nz[1] = glidePlane[gpIndex][Z];
> 
>             //}
> 
> //
>  //*          The initial segments created are not necessarily limited to
>  //*          param->maxSegLen, so a call to InitRemesh() is needed to
>  //*          chop any excessively long segments into proper lengths.
>  //*          When we've generated the nodal data for the final chain,
>  //*          write the current block of nodal data to the file.
>  //*/
>             //InitRemesh(inData, dislocType, startRemeshIndex);
> 
>             //lastBlock = (ic == numChains - 1);
>             //if (lastBlock) {
>                 //IncDislocationDensity(inData, totDislocLen);
>                 //param->nodeCount = inData->nodeCount;
>                 //WriteInitialNodeData(home, inData, lastBlock,nodezero);
>                 //FreeInNodeArray(inData, inData->nodeCount);
>                 //inData->nodeCount = 0;
>             //}
> 
>             //id0 = inData->nodeCount;
>             //startRemeshIndex = id0;
>         //}
>         //return;
> //}
2300,2302d3504
<         glidePlane[2][1][0] =  0.7071068;
<         glidePlane[2][1][1] =  0.0000000;
<         glidePlane[2][1][2] = -0.7071068;
2304,2306c3506,3530
<         glidePlane[2][2][0] =  0.0000000;
<         glidePlane[2][2][1] =  0.7071068;
<         glidePlane[2][2][2] =  0.7071068;
---
> /*---------------------------------------------------------------------------
>  *
>  *      Function:     CreatePrecipitates
>  *      Description:  Creates spherical precipitates in random locations of simulation space. 
>  * 		Defining parameters of the precipitates are  radius r and force parameter forcep.
>  *
>  *      Arguments:
>  *          cubeLength  Length of cubic problem space in a single
>  *                      dimension (units of b)
>  *          numChains   Number of chains to create
>  *          seed        Seed value for random number generator
>  *
>  *-------------------------------------------------------------------------*/
> void CreatePrecipitates(Home_t *home, InPrecipitateData_t *inprecipitateData, int cubeLength,
>                  int numPrecipitates, int seed, real8 r,real8 forcep,
>                  int nodezero)
> {
>         int           ic, ip, np, id0, lastBlock, signFactor;
>         int           newPrecipitateIndex, startRemeshIndex;
>         int           ldIndex, gpIndex, burgIndex, nbr1Index, nbr2Index;
>         real8         posFactor, cubeSize;
>         real8         xp, yp, zp;
>         Param_t       *param;
>         Precipitate_t        *precipitate;
>         
2308,2310d3531
<         glidePlane[2][3][0] = -0.7071068;
<         glidePlane[2][3][1] = -0.7071068;
<         glidePlane[2][3][2] = -0.0000000;
2312,2314c3533,3536
<         glidePlane[2][4][0] = -0.7071068;
<         glidePlane[2][4][1] = -0.0000000;
<         glidePlane[2][4][2] =  0.7071068;
---
>         if (numPrecipitates <= 0) {
>             Fatal("%s: numPrecipitates is %d, but must be > 0.\n",
>                   "CreatePrecipitates", numPrecipitates);
>         }
2316,2318c3538,3545
<         glidePlane[2][5][0] = -0.0000000;
<         glidePlane[2][5][1] = -0.7071068;
<         glidePlane[2][5][2] = -0.7071068;
---
>         param = inprecipitateData->param;
>         cubeSize = (real8)cubeLength;
>         id0 = 0;
>         inprecipitateData->precipitateCount = 0;
>         startRemeshIndex = 0;
>         posFactor = 0.333 * cubeSize;
>         
>         
2321c3548
<  *      glide planes for [1 1 -1]
---
>  *      Create the specified number of chains.
2323,2333c3550,3553
<         glidePlane[3][0][0] =  0.7071068;
<         glidePlane[3][0][1] = -0.7071068;
<         glidePlane[3][0][2] =  0.0000000;
< 
<         glidePlane[3][1][0] =  0.0000000;
<         glidePlane[3][1][1] =  0.7071068;
<         glidePlane[3][1][2] =  0.7071068;
< 
<         glidePlane[3][2][0] =  0.7071068;
<         glidePlane[3][2][1] =  0.0000000;
<         glidePlane[3][2][2] =  0.7071068;
---
>         for (ic = 0; ic < numPrecipitates; ic++) {
> 			
>             
>               
2335,2337d3554
<         glidePlane[3][3][0] = -0.7071068;
<         glidePlane[3][3][1] =  0.7071068;
<         glidePlane[3][3][2] = -0.0000000;
2339,2341c3556,3558
<         glidePlane[3][4][0] = -0.0000000;
<         glidePlane[3][4][1] = -0.7071068;
<         glidePlane[3][4][2] = -0.7071068;
---
>             np = 3;
>             newPrecipitateIndex = inprecipitateData->precipitateCount;
>             inprecipitateData->precipitateCount += 1;
2343,2345c3560,3564
<         glidePlane[3][5][0] = -0.7071068;
<         glidePlane[3][5][1] = -0.0000000;
<         glidePlane[3][5][2] = -0.7071068;
---
>             inprecipitateData->precipitate = (Precipitate_t *)realloc(inprecipitateData->precipitate,
>                            (inprecipitateData->precipitateCount+1) * sizeof(Precipitate_t));
>              
>                            
>             memset(&inprecipitateData->precipitate[newPrecipitateIndex], 0, sizeof(Precipitate_t) );
2348,2350c3567,3569
<  *      Create the specified number of chains.  Anytime the number of
<  *      nodes maintained in memory exceeds the threshhold, write the
<  *      block of nodal data out to the data file.
---
>  *          Set up 3 initial points for the line.  Point 1 is a base position
>  *          at a random location, point 0 is in the negative direction along
>  *          the line and point 2 is in the positive direction along the line.
2352,2353c3571
<         baseNodeID = 0;
<         inData->nodeCount = 0;
---
>  
2355c3573,3580
<         for (chain = 0; chain < numSources; chain++) {
---
>             xp = (randm(&seed)-0.5)*cubeSize;
>             yp = (randm(&seed)-0.5)*cubeSize;
>             zp = (randm(&seed)-0.5)*cubeSize;
>             
>             //modified to produce precipitates in single glide plane (AL)
>             //xp = (randm(&seed)-0.5)*cubeSize*0.8;
>             //yp = (randm(&seed)-0.5)*cubeSize;
>             //zp = -xp;
2357,2361c3582
<             numPoints = 3;
<             lineIndex = chain % 4;
<             burgIndex = lineIndex;
<             gpBurgIndex = lineIndex;
<             gpIndex = (chain / 4) % 6;
---
>            
2364,2365c3585,3586
<  *          Reallocate the node array with sufficient size to add
<  *          all the new nodes defining this chain.
---
>  *          Loop over the points and set up the nodes, link them to 
>  *          the neighbor nodes, etc.
2367,2371c3588,3604
<             newNodeIndex = inData->nodeCount;
<             inData->nodeCount += numPoints;
<             inData->node = (Node_t *)realloc(inData->node, inData->nodeCount
<                                                * sizeof(Node_t));
<             memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numPoints);
---
>             
> 
>                 precipitate = &inprecipitateData->precipitate[ic];
> 
>                 precipitate->x = xp;
>                 precipitate->y = yp;
>                 precipitate->z = zp;
> 
> 				precipitate->forcep = forcep;
>                 precipitate->r = r;
>                 precipitate->constraint=0;
>                 
>                 precipitate->myTag.domainID = 0;
>                 precipitate->myTag.index = ic;
> 				
>             
> 
2374,2375c3607,3611
<  *          Length of the frank read source should be a random length
<  *          between srcLenMin and srcLenMax.
---
>  *          The initial segments created are not necessarily limited to
>  *          param->maxSegLen, so a call to InitRemesh() is needed to
>  *          chop any excessively long segments into proper lengths.
>  *          When we've generated the nodal data for the final chain,
>  *          write the current block of nodal data to the file.
2377,2381c3613
<             if (lenRange > 0) {
<                 srcLen = (srcLenMin + (randm(&seed) * lenRange)) / sqrt(3.0);
<             } else {
<                 srcLen = srcLenMin / sqrt(3.0);
<             }
---
>             
2383,2386c3615,3623
<             for (i = 0; i < 3; i++) {
<                 p0[i] = (randm(&seed)-0.5) * cubeLength * 2.0;
<                 p1[i] = p0[i] - (srcLen * linedir[lineIndex][i]);
<                 p2[i] = p0[i] + (srcLen * linedir[lineIndex][i]);
---
>             lastBlock = (ic == numPrecipitates - 1);
>             if (lastBlock) {
>                 
>                 param->precipitateCount = inprecipitateData->precipitateCount;
>               
>                 WriteInitialPrecipitateData(home, inprecipitateData, lastBlock,nodezero);
>              
>                 FreeInPrecipitateArray(inprecipitateData, inprecipitateData->precipitateCount);
>                 inprecipitateData->precipitateCount = 0;
2389,2393c3626,3629
< /*
<  *          Set up the 3 nodes we're using to define the dislocation line.
<  *          First do point p0.
<  */
<             node = &inData->node[baseNodeID];
---
>            
>         }
>         return;
> }
2395,2399d3630
<             node->myTag.domainID = dislocType;
<             node->myTag.index    = baseNodeID;
<             node->x = p0[0];
<             node->y = p0[1];
<             node->z = p0[2];
2401d3631
<             node->constraint = UNCONSTRAINED;
2403c3633,3693
<             AllocNodeArms(node, 2);
---
> /*---------------------------------------------------------------------------
>  *
>  *      Function:     CreateSingleEdges
>  *      Description:  Creates "edge" dislocations in the [100] [010] and [001]
>  *                    directions.  Each of the three line senses has 4
>  *                    combinations of burgers vector and normals.  With the
>  *                    opposite sign of the line sense vectors considered,
>  *                    the total number of types of dislocations is 24,
>  *                    so number of chains specified should be a multiple
>  *                    of 24 to include all types of.  (Okay, they're
>  *                    not pure edge, but they are not screw either)
>  *
>  *      Arguments:
>  *          cubeLength  Length of cubic problem space in a single
>  *                      dimension (units of b)
>  *          numChains   Number of chains to create
>  *          seed        Seed value for random number generator
>  *
>  *-------------------------------------------------------------------------*/
> void CreateSingleEdges(Home_t *home, InData_t *inData, int cubeLength,
>                  int numChains, int seed, real8 *totDislocLen,
>                  int dislocType,int nodezero)
> {
>         int           ic, ip, np, id0, lastBlock, signFactor;
>         int           newNodeIndex, startRemeshIndex;
>         int           ldIndex, gpIndex, burgIndex, nbr1Index, nbr2Index;
>         real8         posFactor, cubeSize;
>         real8         xp[3], yp[3], zp[3];
>         Param_t       *param;
>         Node_t        *node;
>         
>         //static real8  lineDir[3][3] = {
>                       //{0.0, 0.0, 1.0},
>                       //{0.0, 1.0, 0.0},
>                       //{1.0, 0.0, 0.0}};
>         // ld 1 1 2              
>         static real8  lineDir[3][3] = {
>                       {1.0, 0.0,0.0},
>                       {0.0, 1.0, 0.0},
>                       {1.0, 0.0, 0.0}};              
>                       
>                       
>                       
>         static real8  burg[4][3] = {
>                       { 0.0,  1.0,  0.0},
>                       { 0.0,  0.0, 1.0},
>                       { 0.0, 1.0,  0.0},
>                       {-0.5773503,  0.5773503,  0.5773503}};
>         static real8  glidePlane[12][3] = {
>                       {  0.0, 0.0,  1.0},  /* ldir [001] b [111] 1 */
>                       {  1.0, 0.0,  0.0},  /* ldir [001] b [11-1] 2*/
>                       {  0.0,  0.0,  1.0},  /* ldir [001] b [1-11] 3*/
>                       {  0.7071068,  0.7071068,  0},  /* ldir [001] b [-111] 4*/
>                       {  0.7071068,  0, -0.7071068},  /* ldir [010] b [111]  5*/
>                       {  0.7071068,  0,  0.7071068},  /* ldir [010] b [11-1] 6*/
>                       {  0.7071068,  0, -0.7071068},  /* ldir [010] b [1-11] 7*/
>                       {  0.7071068,  0,  0.7071068},  /* ldir [010] b [-111] 8*/
>                       {  0, -0.7071068,  0.7071068},  /* ldir [100] b [111]  9*/
>                       {  0,  0.7071068,  0.7071068},  /* ldir [100] b [11-1] 10*/
>                       {  0,  0.7071068,  0.7071068},  /* ldir [100] b [1-11] 11*/
>                       {  0, -0.7071068,  0.7071068}}; /* ldir [100] b [-111] 12*/
2405,2408d3694
<             node->nbrTag[0].domainID = dislocType;
<             node->nbrTag[0].index    = baseNodeID + 1;
<             node->nbrTag[1].domainID = dislocType;
<             node->nbrTag[1].index    = baseNodeID + 2;
2410,2412c3696,3699
<             node->burgX[0] = linedir[burgIndex][0];
<             node->burgY[0] = linedir[burgIndex][1];
<             node->burgZ[0] = linedir[burgIndex][2];
---
>         if (numChains <= 0) {
>             Fatal("%s: numChains is %d, but must be > 0.\n",
>                   "CreateEdges", numChains);
>         }
2414,2416c3701,3706
<             node->burgX[1] = -linedir[burgIndex][0];
<             node->burgY[1] = -linedir[burgIndex][1];
<             node->burgZ[1] = -linedir[burgIndex][2];
---
>         param = inData->param;
>         cubeSize = (real8)cubeLength;
>         id0 = 0;
>         inData->nodeCount = 0;
>         startRemeshIndex = 0;
>         posFactor = 0.333 * cubeSize;
2419,2421c3709
<  *          For screw dislocations, use the glide plane from the table.
<  *          For edge dislocations, glide plane is the cross product
<  *          of the burgers vector and line direction.
---
>  *      Create the specified number of chains.
2423,2441c3711,3713
<             if (isScrew) {
<                 node->nx[0] = glidePlane[gpBurgIndex][gpIndex][0];
<                 node->ny[0] = glidePlane[gpBurgIndex][gpIndex][1];
<                 node->nz[0] = glidePlane[gpBurgIndex][gpIndex][2];
<                 node->nx[1] = glidePlane[gpBurgIndex][gpIndex][0];
<                 node->ny[1] = glidePlane[gpBurgIndex][gpIndex][1];
<                 node->nz[1] = glidePlane[gpBurgIndex][gpIndex][2];
<             } else {
<                 cross(linedir[burgIndex], linedir[lineIndex], plane);
<                 plane[0] = (floor(plane[0] * 1.0e+07)) * 1.0e-07;
<                 plane[1] = (floor(plane[1] * 1.0e+07)) * 1.0e-07;
<                 plane[2] = (floor(plane[2] * 1.0e+07)) * 1.0e-07;
<                 node->nx[0] = plane[0];
<                 node->ny[0] = plane[1];
<                 node->nz[0] = plane[2];
<                 node->nx[1] = plane[0];
<                 node->ny[1] = plane[1];
<                 node->nz[1] = plane[2];
<             }
---
>  
> 		for (ic = 0; ic < numChains; ic++) {
>         //for (ic = 0; ic < numChains; ic++) {
2442a3715,3718
>             gpIndex = ic % 2; 
>             burgIndex = gpIndex % 1;
>             ldIndex = gpIndex / 1;
>             
2444c3720,3722
<  *          Now point p1...
---
>  *          First 12 burgers vector/normal sets use positive line
>  *          sense, next set of 12 uses opposite line sense, and
>  *          so on.
2446c3724
<             node = &inData->node[baseNodeID+1];
---
>             signFactor = ((ic / 12) & 0x01) ? -1 : 1;
2448,2449c3726,3728
<             node->myTag.domainID = dislocType;
<             node->myTag.index    = baseNodeID+1;
---
>             np = 3;
>             newNodeIndex = inData->nodeCount;
>             inData->nodeCount += np;
2451,2453c3730,3732
<             node->x = p1[0];
<             node->y = p1[1];
<             node->z = p1[2];
---
>             inData->node = (Node_t *)realloc(inData->node,
>                            inData->nodeCount * sizeof(Node_t));
>             memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);
2455c3734,3746
<             node->constraint = PINNED_NODE;
---
> /*
>  *          Set up 3 initial points for the line.  Point 1 is a base position
>  *          at a random location, point 0 is in the negative direction along
>  *          the line and point 2 is in the positive direction along the line.
>  */
>             //xp[1] = (randm(&seed)-0.5)*cubeSize;
>             //yp[1] = (randm(&seed)-0.5)*cubeSize;
>             //zp[1] = (randm(&seed)-0.5)*cubeSize;
>             
>             //for single glide plane (AL)
>             xp[1] = (randm(&seed)-0.5)*cubeSize;
>             yp[1] = 0.0;
>             zp[1] = -1000 *(ldIndex%2) ;
2457c3748
<             AllocNodeArms(node, 1);
---
> //Teeppä tästä kerrostalo versio
2459,2460d3749
<             node->nbrTag[0].domainID = dislocType;
<             node->nbrTag[0].index    = baseNodeID;
2462,2467c3751,3757
<             node->burgX[0] = -linedir[burgIndex][0];
<             node->burgY[0] = -linedir[burgIndex][1];
<             node->burgZ[0] = -linedir[burgIndex][2];
<             node->nx[0] = glidePlane[gpBurgIndex][gpIndex][0];
<             node->ny[0] = glidePlane[gpBurgIndex][gpIndex][1];
<             node->nz[0] = glidePlane[gpBurgIndex][gpIndex][2];
---
>             xp[0] = xp[1] - (posFactor * signFactor * lineDir[ldIndex][X]);
>             yp[0] = yp[1] - (posFactor * signFactor * lineDir[ldIndex][Y]);
>             zp[0] = zp[1] - (posFactor * signFactor * lineDir[ldIndex][Z]);
> 
>             xp[2] = xp[1] + (posFactor * signFactor * lineDir[ldIndex][X]);
>             yp[2] = yp[1] + (posFactor * signFactor * lineDir[ldIndex][Y]);
>             zp[2] = zp[1] + (posFactor * signFactor * lineDir[ldIndex][Z]);
2470c3760,3761
<  *          Now point p2...
---
>  *          Loop over the points and set up the nodes, link them to 
>  *          the neighbor nodes, etc.
2472c3763
<             node = &inData->node[baseNodeID+2];
---
>             for (ip = 0; ip < np; ip++) {
2474,2475c3765
<             node->myTag.domainID = dislocType;
<             node->myTag.index    = baseNodeID+2;
---
>                 node = &inData->node[ip+id0];
2477,2479c3767,3769
<             node->x = p2[0];
<             node->y = p2[1];
<             node->z = p2[2];
---
>                 node->x = xp[ip];
>                 node->y = yp[ip];
>                 node->z = zp[ip];
2481c3771,3773
<             node->constraint = PINNED_NODE;
---
>                 node->constraint = UNCONSTRAINED;
>                 node->myTag.domainID = dislocType;
>                 node->myTag.index = ip+id0;
2483c3775
<             AllocNodeArms(node, 1);
---
>                 AllocNodeArms(node, 2);
2485,2486c3777,3778
<             node->nbrTag[0].domainID = dislocType;
<             node->nbrTag[0].index    = baseNodeID;
---
>                 if ((nbr1Index = ip + 1) >= np) nbr1Index = 0;
>                 if ((nbr2Index= ip - 1) < 0) nbr2Index = np - 1;
2488,2493c3780,3798
<             node->burgX[0] = linedir[burgIndex][0];
<             node->burgY[0] = linedir[burgIndex][1];
<             node->burgZ[0] = linedir[burgIndex][2];
<             node->nx[0] = glidePlane[gpBurgIndex][gpIndex][0];
<             node->ny[0] = glidePlane[gpBurgIndex][gpIndex][1];
<             node->nz[0] = glidePlane[gpBurgIndex][gpIndex][2];
---
>                 node->nbrTag[0].domainID = dislocType;
>                 node->nbrTag[0].index = id0 + nbr1Index;
>                 node->burgX[0] = burg[burgIndex][0];
>                 node->burgY[0] = burg[burgIndex][1];
>                 node->burgZ[0] = burg[burgIndex][2];
>                 node->nx[0] = glidePlane[gpIndex][X];
>                 node->ny[0] = glidePlane[gpIndex][Y];
>                 node->nz[0] = glidePlane[gpIndex][Z];
>             
>                 node->nbrTag[1].domainID = dislocType;
>                 node->nbrTag[1].index = id0 + nbr2Index;
>                 node->burgX[1] = -burg[burgIndex][0];
>                 node->burgY[1] = -burg[burgIndex][1];
>                 node->burgZ[1] = -burg[burgIndex][2];
>                 node->nx[1] = glidePlane[gpIndex][X];
>                 node->ny[1] = glidePlane[gpIndex][Y];
>                 node->nz[1] = glidePlane[gpIndex][Z];
> 
>             }
2500c3805
<  *          write the block of nodal data to the file.
---
>  *          write the current block of nodal data to the file.
2504,2505c3809
<             lastBlock = (chain == (numSources - 1));
< 
---
>             lastBlock = (ic == numChains - 1);
2507,2511c3811,3815
<                     IncDislocationDensity(inData, totDislocLen);
<                     param->nodeCount = inData->nodeCount;
<                     WriteInitialNodeData(home, inData, lastBlock);
<                     FreeInNodeArray(inData, inData->nodeCount);
<                     inData->nodeCount = 0;
---
>                 IncDislocationDensity(inData, totDislocLen);
>                 param->nodeCount = inData->nodeCount;
>                 WriteInitialNodeData(home, inData, lastBlock,nodezero);
>                 FreeInNodeArray(inData, inData->nodeCount);
>                 inData->nodeCount = 0;
2514,2518c3818,3820
<             baseNodeID = inData->nodeCount;
<             startRemeshIndex = baseNodeID;
< 
<         }  /* for (chain = 0; chain < numChains; ...) */
< 
---
>             id0 = inData->nodeCount;
>             startRemeshIndex = id0;
>         }
2522d3823
< 
2525c3826
<  *      Function:     CreateEdges
---
>  *      Function:     Create2DEdges
2542c3843
< void CreateEdges(Home_t *home, InData_t *inData, int cubeLength,
---
> void CreateD2Edges(Home_t *home, InData_t *inData, int cubeLength,
2544c3845
<                  int dislocType)
---
>                  int dislocType,int nodezero)
2552a3854,3859
>         
>         //static real8  lineDir[3][3] = {
>                       //{0.0, 0.0, 1.0},
>                       //{0.0, 1.0, 0.0},
>                       //{1.0, 0.0, 0.0}};
>         // ld 1 1 2              
2554,2556c3861,3866
<                       {0.0, 0.0, 1.0},
<                       {0.0, 1.0, 0.0},
<                       {1.0, 0.0, 0.0}};
---
>                       {1.0, 0.0,0.0},
>                       {1.0, 0.0, 0.0},
>                       {1.0, 0.0, 0.0}};              
>                       
>                       
>                       
2558,2561c3868,3871
<                       { 0.5773503,  0.5773503,  0.5773503},
<                       { 0.5773503,  0.5773503, -0.5773503},
<                       { 0.5773503, -0.5773503,  0.5773503},
<                       {-0.5773503,  0.5773503,  0.5773503}};
---
>                       { 0.0,  1.0,  0.0},
>                       { 0.0,  -1.0, 0.0},
>                       { 0.0, 1.0,  0.0},
>                       { 0.0, -1.0,  0.0}};
2563,2574c3873,3884
<                       {  0.7071068, -0.7071068,  0},  /* ldir [001] b [111]  */
<                       {  0.7071068, -0.7071068,  0},  /* ldir [001] b [11-1] */
<                       {  0.7071068,  0.7071068,  0},  /* ldir [001] b [1-11] */
<                       {  0.7071068,  0.7071068,  0},  /* ldir [001] b [-111] */
<                       {  0.7071068,  0, -0.7071068},  /* ldir [010] b [111]  */
<                       {  0.7071068,  0,  0.7071068},  /* ldir [010] b [11-1] */
<                       {  0.7071068,  0, -0.7071068},  /* ldir [010] b [1-11] */
<                       {  0.7071068,  0,  0.7071068},  /* ldir [010] b [-111] */
<                       {  0, -0.7071068,  0.7071068},  /* ldir [100] b [111]  */
<                       {  0,  0.7071068,  0.7071068},  /* ldir [100] b [11-1] */
<                       {  0,  0.7071068,  0.7071068},  /* ldir [100] b [1-11] */
<                       {  0, -0.7071068,  0.7071068}}; /* ldir [100] b [-111] */
---
>                       {  0.0, 0.0,  1.0},  /* ldir [001] b [111] 1 */
>                       {  0.0, 0.0,  1.0},  /* ldir [001] b [11-1] 2*/
>                       {  0.0,  0.0,  1.0},  /* ldir [001] b [1-11] 3*/
>                       {  0.0,  0.0,  1.0},  /* ldir [001] b [-111] 4*/
>                       {  0.0,  0.0,  1.0},  /* ldir [010] b [111]  5*/
>                       {  0.0,  0.0,  1.0},  /* ldir [010] b [11-1] 6*/
>                       {  0.0,  0.0,  1.0},  /* ldir [010] b [1-11] 7*/
>                       {  0.0,  0.0,  1.0},  /* ldir [010] b [-111] 8*/
>                       {  0.0,  0.0,  1.0},  /* ldir [100] b [111]  9*/
>                       {  0.0,  0.0,  1.0},  /* ldir [100] b [11-1] 10*/
>                       {  0.0,  0.0,  1.0},  /* ldir [100] b [1-11] 11*/
>                       {  0.0,  0.0,  1.0}}; /* ldir [100] b [-111] 12*/
2592c3902,3904
<         for (ic = 0; ic < numChains; ic++) {
---
>  
> 		for (ic = 0; ic < numChains; ic++) {
>         //for (ic = 0; ic < numChains; ic++) {
2603c3915
<             signFactor = ((ic / 12) & 0x01) ? -1 : 1;
---
>             signFactor = ((ic / 2) & 0x01) ? -1 : 1;
2620a3933,3937
>             
>             
> 
> //Teeppä tästä kerrostalo versio
> 
2684c4001
<                 WriteInitialNodeData(home, inData, lastBlock);
---
>                 WriteInitialNodeData(home, inData, lastBlock,nodezero);
