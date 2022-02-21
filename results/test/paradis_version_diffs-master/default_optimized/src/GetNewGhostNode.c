54a55,61
> 		
> 		newNode->nodedx= 0.0;
>         newNode->nodedy= 0.0;
>         newNode->nodedz= 0.0;
> 
> 
> 
70a78,158
> 
> /***************************************************************************
>  *
>  *  Function    : GetNewGhostPrecipitate
>  *  Description : Get a free precipitate, and assign it to the specified domain 
>  *                and index. Extend the RemoteDomain's precipitateKeys table, if
>  *                necessary. Also, queue the precipitate onto the ghost precipitate queue
>  *
>  ***************************************************************************/
> 
> #include "Home.h"
> #include "QueueOps.h"
> #include "Util.h"
> 
> Precipitate_t *GetNewGhostPrecipitate(Home_t *home, int domain, int index)
> {
>         int            i;
>         Precipitate_t         *newPrecipitate;
>         RemoteDomain_t *remDom;
> 
>         newPrecipitate = PopFreePrecipitateQ (home);
>         remDom  = home->remoteDomainKeys[domain];
> 
>         if (!remDom) {
>             Fatal("GetNewGhostPrecipitate: domain %d is not a neighbor of domain %d",
>                   domain, home->myDomain);
>         }
> 
> /*
>  *      If index is beyond current limit on remDom->precipitateKeys, expand it big
>  *      enough to include index, and initialize new elements to zero
>  */
>         if (index >= remDom->maxprecipitateTagIndex) {
>             remDom->precipitateKeys = (Precipitate_t **)realloc(remDom->precipitateKeys, 
>                                                   (index+1)*sizeof(Precipitate_t *));
> 
>         for (i = remDom->maxprecipitateTagIndex; i < index; i++)
>             remDom->precipitateKeys[i] = 0;
>             remDom->maxprecipitateTagIndex = index + 1;
>         }
> 
>         remDom->precipitateKeys[index] = newPrecipitate;
>         newPrecipitate->myTag.domainID = domain;
>         newPrecipitate->myTag.index    = index;
>         newPrecipitate->cellIdx        = -1;
>         newPrecipitate->cell2Idx       = -1;
>         newPrecipitate->cell2QentIdx   = -1;
> 
>         newPrecipitate->vX = 0.0;
>         newPrecipitate->vY = 0.0;
>         newPrecipitate->vZ = 0.0;
> 
>         newPrecipitate->oldvX = 0.0;
>         newPrecipitate->oldvY = 0.0;
>         newPrecipitate->oldvZ = 0.0;
> 
>         newPrecipitate->flags = 0;
> 
>         PushGhostPrecipitateQ(home, newPrecipitate);
> 
> #ifdef _FEM
>         newPrecipitate->fem_Surface[0] = 0;
>         newPrecipitate->fem_Surface[1] = 0;
> 
>         newPrecipitate->fem_Surface_Norm[0] = 0.0;
>         newPrecipitate->fem_Surface_Norm[1] = 0.0;
>         newPrecipitate->fem_Surface_Norm[2] = 0.0;
> #endif
> 
>         return(newPrecipitate);
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
