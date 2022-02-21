41c41
<         int     i, j, nbrs;
---
>         int     i, j, nbrs,newPrecipitateKeyPtr;
48a49
>         real8   xp,yp,zp,r;
63a65
>         Precipitate_t *precipitate;
117a120,146
> 
> /* If node encounters a precipitate (AL)*/
> //newPrecipitateKeyPtr= home->newPrecipitateKeyPtr;
>         
>         //for (i = 0; i < newPrecipitateKeyPtr; i++) {
> 			
> 			//if ((precipitate = home->precipitateKeys[i]) == (Precipitate_t *)NULL) {
>                 //continue;
>             //}
> 
> 		//xp=node->x;
> 		//yp=node->y;
> 		//zp=node->z;
> 		//r=sqrt((xp-precipitate->x)*(xp-precipitate->x)+(yp-precipitate->y)*(yp-precipitate->y)+(zp-precipitate->z)*(zp-precipitate->z));
> 		////printf("f %e   \n",node->fX);
> 		//if (r<=precipitate->r) {
>             //node->vX =0;
>             //node->vY = 0;
>             //node->vZ = 0;
>             //return(0);
>         //}
> 
> //}
> 
> 
> 
> 
297c326
<         nForce[0] = node->fX;
---
> 		nForce[0] = node->fX;
299a329,343
> 
> 
> 
> 		//if (r<=2000.0) {
> 		//nForce[0] = node->fX+1e+13;
>         //nForce[1] = node->fY+1e+13;
>         //nForce[2] = node->fZ;	
> 			
> 		//}
> 		
>      //printf("f %e %e %e   \n",node->fX,node->fY,node->fZ);   
>         
>         
>         
>         
