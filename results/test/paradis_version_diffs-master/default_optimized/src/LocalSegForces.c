2385a2386,2387
> 
> 
2431a2434
>         
2747a2751
>  
2766a2771,2772
>             real8	vx1,vy1,vz1;
>             real8	vx2,vy2,vz2;
2798c2804,2805
< 
---
>             //printf("homedomain,node1,cell %d ( %d %d ) %d %d\n ",home->myDomain,node1->myTag.domainID,node1->myTag.index,node1->cellIdx,node1->numPNbrs);	
> 			//printf("homedomain,node2,cell %d ( %d %d ) %d \n ",home->myDomain,node2->myTag.domainID,node2->myTag.index,node2->cellIdx);
2802c2809,2818
< 
---
>             
> 			vx1=node1->vX;
> 			vy1=node1->vY;
> 			vz1=node1->vZ;
> 			
> 			vx2=node2->vX;
> 			vy2=node2->vY;
> 			vz2=node2->vZ;
> 			//printf("Flag %e",vx2);
> 			
2821a2838,2840
> 				//printf("Nodeforce %d \n",home->myDomain);
> 				//PrintNode(node1);
> 				//PrintNode(node2);
2864a2884,2885
> 		//printf("Flag %e",f2[1]);
> 
2880a2902,2918
>            
>             
> 				//printf("before disorder,node1,node2  ( %d %d ) (%d %d)\n ",node1->myTag.domainID,node1->myTag.index,node2->myTag.domainID,node2->myTag.index);
>                 DisorderForce(home,node1,node2, f1, f2);
> 				//printf("After disorder,node1,node2  ( %d %d ) (%d %d)\n ",node1->myTag.domainID,node1->myTag.index,node2->myTag.domainID,node2->myTag.index);
> 
>                 VECTOR_ADD(node1SegForce, f1);
>                 VECTOR_ADD(node2SegForce, f2);
> 
>                 for (j = 0; j < 3; j++) {
>                     nativeSegList[i].seg->f1[j] += f1[j];
>                     nativeSegList[i].seg->f2[j] += f2[j];
>                 }
>             
>             
>             
>             
2932c2970,2971
< 
---
>  
>   
3010c3049
< 
---
>   
3037c3076
< 
---
> //printf("In localforce  before MPI_COMM_world %d cycle %d \n",home->myDomain,home->cycle+1);
3041a3081
> //printf("In localforce  after MPI_COMM_world %d cycle %d \n",home->myDomain,home->cycle+1);
