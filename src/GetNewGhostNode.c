/***************************************************************************
 *
 *  Function    : GetNewGhostNode
 *  Description : Get a free node, and assign it to the specified domain 
 *                and index. Extend the RemoteDomain's nodeKeys table, if
 *                necessary. Also, queue the node onto the ghost node queue
 *
 ***************************************************************************/

#include "Home.h"
#include "QueueOps.h"
#include "Util.h"

Node_t *GetNewGhostNode(Home_t *home, int domain, int index)
{
        int            i;
        Node_t         *newNode;
        RemoteDomain_t *remDom;

        newNode = PopFreeNodeQ (home);
        remDom  = home->remoteDomainKeys[domain];

        if (!remDom) {
            Fatal("GetNewGhostNode: domain %d is not a neighbor of domain %d",
                  domain, home->myDomain);
        }

/*
 *      If index is beyond current limit on remDom->nodeKeys, expand it big
 *      enough to include index, and initialize new elements to zero
 */
        if (index >= remDom->maxTagIndex) {
            remDom->nodeKeys = (Node_t **)realloc(remDom->nodeKeys, 
                                                  (index+1)*sizeof(Node_t *));

        for (i = remDom->maxTagIndex; i < index; i++)
            remDom->nodeKeys[i] = 0;
            remDom->maxTagIndex = index + 1;
        }

        remDom->nodeKeys[index] = newNode;
        newNode->myTag.domainID = domain;
        newNode->myTag.index    = index;
        newNode->cellIdx        = -1;
        newNode->cell2Idx       = -1;
        newNode->cell2QentIdx   = -1;

        newNode->vX = 0.0;
        newNode->vY = 0.0;
        newNode->vZ = 0.0;

        newNode->oldvX = 0.0;
        newNode->oldvY = 0.0;
        newNode->oldvZ = 0.0;
		newNode->nodedx= 0.0;
        newNode->nodedy= 0.0;
        newNode->nodedz= 0.0;

        newNode->flags = 0;

        PushGhostNodeQ(home, newNode);

#ifdef _FEM
        newNode->fem_Surface[0] = 0;
        newNode->fem_Surface[1] = 0;

        newNode->fem_Surface_Norm[0] = 0.0;
        newNode->fem_Surface_Norm[1] = 0.0;
        newNode->fem_Surface_Norm[2] = 0.0;
#endif

        return(newNode);
}


/***************************************************************************
 *
 *  Function    : GetNewGhostPrecipitate
 *  Description : Get a free precipitate, and assign it to the specified domain 
 *                and index. Extend the RemoteDomain's precipitateKeys table, if
 *                necessary. Also, queue the precipitate onto the ghost precipitate queue
 *
 ***************************************************************************/

#include "Home.h"
#include "QueueOps.h"
#include "Util.h"

Precipitate_t *GetNewGhostPrecipitate(Home_t *home, int domain, int index)
{
        int            i;
        Precipitate_t         *newPrecipitate;
        RemoteDomain_t *remDom;

        newPrecipitate = PopFreePrecipitateQ (home);
        remDom  = home->remoteDomainKeys[domain];

        if (!remDom) {
            Fatal("GetNewGhostPrecipitate: domain %d is not a neighbor of domain %d",
                  domain, home->myDomain);
        }

/*
 *      If index is beyond current limit on remDom->precipitateKeys, expand it big
 *      enough to include index, and initialize new elements to zero
 */
        if (index >= remDom->maxprecipitateTagIndex) {
            remDom->precipitateKeys = (Precipitate_t **)realloc(remDom->precipitateKeys, 
                                                  (index+1)*sizeof(Precipitate_t *));

        for (i = remDom->maxprecipitateTagIndex; i < index; i++)
            remDom->precipitateKeys[i] = 0;
            remDom->maxprecipitateTagIndex = index + 1;
        }

        remDom->precipitateKeys[index] = newPrecipitate;
        newPrecipitate->myTag.domainID = domain;
        newPrecipitate->myTag.index    = index;
        newPrecipitate->cellIdx        = -1;
        newPrecipitate->cell2Idx       = -1;
        newPrecipitate->cell2QentIdx   = -1;

        newPrecipitate->vX = 0.0;
        newPrecipitate->vY = 0.0;
        newPrecipitate->vZ = 0.0;

        newPrecipitate->oldvX = 0.0;
        newPrecipitate->oldvY = 0.0;
        newPrecipitate->oldvZ = 0.0;

        newPrecipitate->flags = 0;

        PushGhostPrecipitateQ(home, newPrecipitate);

#ifdef _FEM
        newPrecipitate->fem_Surface[0] = 0;
        newPrecipitate->fem_Surface[1] = 0;

        newPrecipitate->fem_Surface_Norm[0] = 0.0;
        newPrecipitate->fem_Surface_Norm[1] = 0.0;
        newPrecipitate->fem_Surface_Norm[2] = 0.0;
#endif

        return(newPrecipitate);
}
