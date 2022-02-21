664c664
<  *      Description:  This function does no force calulcations directly
---
>  *      Description:  This function does no force calculations directly
704a705
> 				 
705a707
>                 
853a856
>   
855a859
>         
1187a1192,1193
>             real8	vx1,vy1,vz1;
>             real8	vx2,vy2,vz2;
1286a1293,1300
> 			vx1=node1->vX;
> 			vy1=node1->vY;
> 			vz1=node1->vZ;
> 			
> 			vx2=node2->vX;
> 			vy2=node2->vY;
> 			vz2=node2->vZ;
> 
1291a1306
> 				
1331a1347
> 			
1357a1374,1385
>             
>             
> 				DisorderForce(home,node1,node2, f1, f2);
>                 
> 
>                 VECTOR_ADD(fseg1node1, f1);
>                 VECTOR_ADD(fseg1node2, f2);
>             
>             
>             
>             
>             
