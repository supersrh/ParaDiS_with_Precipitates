503c503
<         int     i, thisDomain, splitOK, splitSeg2;
---
>         int     i, thisDomain, splitOK, splitSeg2,j,out;
508a509
>         real8	dx,dy,dz,x1,y1,z1,x2,y2,z2,prx,pry,pz;
568c569,594
< 
---
> 				
> 				x1=node->x;
> 				y1=node->y;
> 				z1=node->z;
> 				
> 				x2=nbr1->x;
> 				y2=nbr1->y;
> 				z2=nbr1->z;
> 				
> 				//printf("nodeinprecipitate %d %e %e %e \n",NodeinPrecipitate(1,5,-5,10,x1, y1,z1,x2,y2,z2),x1, y1,z1);
> 				
> 				
> 				
> 			//for (j = 0; j < newPrecipitateKeyPtr; j++) {
> 			
> 			//if ((precipitate = home->precipitateKeys[j]) == (Precipitate_t *)NULL) {
>                 //continue;
>                 
>                 
>             //}
>             //prx=precipitate->x;
>             //prx=precipitate->y;
>             //prx=precipitate->z;
>             
> 			//}
> 				
