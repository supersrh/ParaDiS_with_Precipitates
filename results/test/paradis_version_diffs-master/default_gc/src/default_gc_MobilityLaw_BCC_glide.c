46c46
<         int     i, nbrs;
---
>         int     i, nbrs,newPrecipitateKeyPtr;
52a53
>         real8 	xp,yp,zp,r;
72a74
>         Precipitate_t *precipitate;
83a86,112
> 
> 
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
