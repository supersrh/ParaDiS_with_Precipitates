46c46
<         int     color, transmit, radius, newNodeKeyPtr;
---
>         int     color, transmit, radius, newNodeKeyPtr,radius2,newPrecipitateKeyPtr;
52a53
>         Precipitate_t *precipitate;
177c178
<                          
---
>                     
185a187,189
>                     
>                          
>                    
188a193,211
> 
> 
>  newPrecipitateKeyPtr= home->newPrecipitateKeyPtr;
>         
>         for (i = 0; i < newPrecipitateKeyPtr; i++) {
> 			
> 			if ((precipitate = home->precipitateKeys[i]) == (Precipitate_t *)NULL) {
>                 continue;
>             }
> 
>       fprintf(fp, "sphere{ <%e,%e,%e>,%e "
>                             "pigment { color color%02d "
>                             "transmit transmit01 } "
>                             "finish { phong 1 metallic }"
>                             " }\n",
>                             precipitate->x, precipitate->y, precipitate->z,
>                             precipitate->r,4,transmit);
>                             
> 		}
