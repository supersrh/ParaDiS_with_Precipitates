/*****************************************************************************
 *
 *  Module      : QueueOps
 *  Description : Contains the various queue operations on the native, ghost,
 *                and free queues.
 *
 ****************************************************************************/

#include "Home.h"
#include "ParadisThread.h"
#include "Node.h"
#include "Precipitate.h"
#include "QueueOps.h"

/*---------------------------------------------------------------------------
 *
 * Function    : PopFreeNodeQ
 * Description : dequeue the node at the top of the freeNodeQ and return it
 *               to the caller. If there are no free nodes, allocate a new
 *               block of nodes, queue them all on the freeNodeQ, and then
 *               pop the first one
 *
 *---------------------------------------------------------------------------*/
Node_t *PopFreeNodeQ(Home_t *home)
{
	int		i;
	NodeBlock_t	*nodeBlock;
	Node_t		*currNode, *freeNode;

	if (home->freeNodeQ == 0) {

		nodeBlock = (NodeBlock_t *) malloc(sizeof(NodeBlock_t));
		nodeBlock->next = home->nodeBlockQ;
		home->nodeBlockQ = nodeBlock;

		nodeBlock->nodes = (Node_t *)calloc(1, NODE_BLOCK_COUNT *
						    sizeof(Node_t));

		currNode = nodeBlock->nodes;
		home->freeNodeQ = currNode;

		for (i = 1; i < NODE_BLOCK_COUNT; i++) {
			INIT_LOCK(&currNode->nodeLock);
			currNode->next = currNode + 1;
			currNode++;
		}

		currNode->next = 0;    /* last free node */
		INIT_LOCK(&currNode->nodeLock);
		home->lastFreeNode = currNode;
	}

/*
 *	Dequeue the first free node and return it to the caller
 */
	freeNode = home->freeNodeQ;
	home->freeNodeQ = freeNode->next;

	return(freeNode);
}


/*--------------------------------------------------------------------------
 *
 *	Function:	PushFreeNodeQ
 *	Description:	Push a node onto the top of the free queue
 *
 *-------------------------------------------------------------------------*/
void PushFreeNodeQ(Home_t *home, Node_t *node)
{
	node->next = home->freeNodeQ;
	home->freeNodeQ = node;

	return;
}


/*--------------------------------------------------------------------------
 *
 *	Function:	PushNativeNodeQ
 *	Description:	Push a node onto the top of the native queue
 *
 *-------------------------------------------------------------------------*/
void PushNativeNodeQ(Home_t *home, Node_t *node)
{
	node->next = home->nativeNodeQ;
	home->nativeNodeQ = node;

	return;
}


/*--------------------------------------------------------------------------
 *
 *	Function:	PushGhostNodeQ
 *	Description:	Push a node onto the top of the ghost queue
 *
 *-------------------------------------------------------------------------*/
void PushGhostNodeQ(Home_t *home, Node_t *node)
{
	if (home->ghostNodeQ == 0) home->lastGhostNode = node;
	node->next = home->ghostNodeQ;
	home->ghostNodeQ = node;

	return;
}


/*--------------------------------------------------------------------------
 *
 *	Function:	RecycleGhostNodes
 *	Description:	Put all the nodes in the ghostNodeQ back onto
 *			the freeNodeQ and reset the ghostNodeQ to empty
 *
 *-------------------------------------------------------------------------*/
void RecycleGhostNodes(Home_t *home)
{
	if (home->ghostNodeQ == NULL) return;  /* nothing to do */

	if (home->freeNodeQ == NULL) {
		home->freeNodeQ = home->ghostNodeQ;
	} else {
		home->lastFreeNode->next = home->ghostNodeQ;
	}

	home->lastFreeNode = home->lastGhostNode;
	home->lastGhostNode = NULL;
	home->ghostNodeQ = NULL;

	return;
}


/*--------------------------------------------------------------------------
 *
 *	Function:	RemoveNodeFromCellQ
 *	Description:	Remove the specified node from the queue for
 *                      the cell containing the node.  Typically only
 *                      needed when deleting a node during topological
 *                      changes.
 *
 *-------------------------------------------------------------------------*/
void RemoveNodeFromCellQ(Home_t *home, Node_t *node)
{
        Cell_t *cell;
        Node_t *prevNode, *nextNode;

        if (node == (Node_t *)NULL) return;
        if (node->cellIdx < 0) return;

        cell = home->cellKeys[node->cellIdx];

        if (node == cell->nodeQ) {
            cell->nodeQ = node->nextInCell;
            cell->nodeCount--;
        } else {
            nextNode = cell->nodeQ;
            while (nextNode != (Node_t *)NULL) {
                if (nextNode == node) {
                    prevNode->nextInCell = node->nextInCell;
                    cell->nodeCount--;
                    return;
                }
                prevNode = nextNode;
                nextNode = nextNode->nextInCell;
            }
        }

        return;
}


/*--------------------------------------------------------------------------
 *
 *	Function:	RemoveNodeFromCell2Q
 *	Description:	Remove the specified node from the queue for
 *                      the cell2 containing the node.  Typically only
 *                      needed when deleting a node during topological
 *                      changes.
 *
 *-------------------------------------------------------------------------*/
void RemoveNodeFromCell2Q(Home_t *home, Node_t *node)
{

	if (home->cell2QentArray == (C2Qent_t *)NULL) {
		return;
	}

        if (node == (Node_t *)NULL) {
		return;
	}

        if ((node->cell2Idx < 0) || (node->cell2QentIdx < 0)) {
		return;
	}

        home->cell2QentArray[node->cell2QentIdx].node = (Node_t *)NULL;
        home->cell2QentArray[node->cell2QentIdx].next = -1;

        return;
}
 
 Precipitate_t *PopFreePrecipitateQ(Home_t *home)
 {
 	int		i;
 	PrecipitateBlock_t	*precipitateBlock;
 	Precipitate_t		*currPrecipitate, *freePrecipitate;
 
 	if (home->freePrecipitateQ == 0) {
 
 		precipitateBlock = (PrecipitateBlock_t *) malloc(sizeof(PrecipitateBlock_t));
 		precipitateBlock->next = home->precipitateBlockQ;
 		home->precipitateBlockQ = precipitateBlock;
 
 		precipitateBlock->precipitates = (Precipitate_t *)calloc(1, PRECIPITATE_BLOCK_COUNT *
 						    sizeof(Precipitate_t));
 
 		currPrecipitate = precipitateBlock->precipitates;
 		home->freePrecipitateQ = currPrecipitate;
 
 		for (i = 1; i < PRECIPITATE_BLOCK_COUNT; i++) {
 			INIT_LOCK(&currPrecipitate->precipitateLock);
 			currPrecipitate->next = currPrecipitate + 1;
 			currPrecipitate++;
 		}
 
 		currPrecipitate->next = 0;    /* last free precipitate */
 		INIT_LOCK(&currPrecipitate->precipitateLock);
 		home->lastFreePrecipitate = currPrecipitate;
 	}
 
 /*
  *	Dequeue the first free precipitate and return it to the caller
  */
 	freePrecipitate = home->freePrecipitateQ;
 	home->freePrecipitateQ = freePrecipitate->next;
 
 	return(freePrecipitate);
 }
 
 
 /*--------------------------------------------------------------------------
  *
  *	Function:	PushFreePrecipitateQ
  *	Description:	Push a precipitate onto the top of the free queue
  *
  *-------------------------------------------------------------------------*/
 void PushFreePrecipitateQ(Home_t *home, Precipitate_t *precipitate)
 {
 	precipitate->next = home->freePrecipitateQ;
 	home->freePrecipitateQ = precipitate;
 
 	return;
 }
 
 
 /*--------------------------------------------------------------------------
  *
  *	Function:	PushNativePrecipitateQ
  *	Description:	Push a precipitate onto the top of the native queue
  *
  *-------------------------------------------------------------------------*/
 void PushNativePrecipitateQ(Home_t *home, Precipitate_t *precipitate)
 {
 	precipitate->next = home->nativePrecipitateQ;
 	home->nativePrecipitateQ = precipitate;
 
 	return;
 }
 
 
 /*--------------------------------------------------------------------------
  *
  *	Function:	PushGhostPrecipitateQ
  *	Description:	Push a precipitate onto the top of the ghost queue
  *
  *-------------------------------------------------------------------------*/
 void PushGhostPrecipitateQ(Home_t *home, Precipitate_t *precipitate)
 {
 	if (home->ghostPrecipitateQ == 0) home->lastGhostPrecipitate = precipitate;
 	precipitate->next = home->ghostPrecipitateQ;
 	home->ghostPrecipitateQ = precipitate;
 
 	return;
 }
 
 
 /*--------------------------------------------------------------------------
  *
  *	Function:	RecycleGhostPrecipitates
  *	Description:	Put all the precipitates in the ghostPrecipitateQ back onto
  *			the freePrecipitateQ and reset the ghostPrecipitateQ to empty
  *
  *-------------------------------------------------------------------------*/
 void RecycleGhostPrecipitates(Home_t *home)
 {
 	if (home->ghostPrecipitateQ == NULL) return;  /* nothing to do */
 
 	if (home->freePrecipitateQ == NULL) {
 		home->freePrecipitateQ = home->ghostPrecipitateQ;
 	} else {
 		home->lastFreePrecipitate->next = home->ghostPrecipitateQ;
 	}
 
 	home->lastFreePrecipitate = home->lastGhostPrecipitate;
 	home->lastGhostPrecipitate = NULL;
 	home->ghostPrecipitateQ = NULL;
 
 	return;
 }
 
 
 /*--------------------------------------------------------------------------
  *
  *	Function:	RemovePrecipitateFromCellQ
  *	Description:	Remove the specified precipitate from the queue for
  *                      the cell containing the precipitate.  Typically only
  *                      needed when deleting a precipitate during topological
  *                      changes.
  *
  *-------------------------------------------------------------------------*/
 void RemovePrecipitateFromCellQ(Home_t *home, Precipitate_t *precipitate)
 {
         Cell_t *cell;
         Precipitate_t *prevPrecipitate, *nextPrecipitate;
 
         if (precipitate == (Precipitate_t *)NULL) return;
         if (precipitate->cellIdx < 0) return;
 
         cell = home->cellKeys[precipitate->cellIdx];
 
         if (precipitate == cell->precipitateQ) {
             cell->precipitateQ = precipitate->nextInCell;
             cell->precipitateCount--;
         } else {
             nextPrecipitate = cell->precipitateQ;
             while (nextPrecipitate != (Precipitate_t *)NULL) {
                 if (nextPrecipitate == precipitate) {
                     prevPrecipitate->nextInCell = precipitate->nextInCell;
                     cell->precipitateCount--;
                     return;
                 }
                 prevPrecipitate = nextPrecipitate;
                 nextPrecipitate = nextPrecipitate->nextInCell;
             }
         }
 
         return;
 }
 
 
 /*--------------------------------------------------------------------------
  *
  *	Function:	RemovePrecipitateFromCell2Q
  *	Description:	Remove the specified precipitate from the queue for
  *                      the cell2 containing the precipitate.  Typically only
  *                      needed when deleting a precipitate during topological
  *                      changes.
  *
  *-------------------------------------------------------------------------*/
 void RemovePrecipitateFromCell2Q(Home_t *home, Precipitate_t *precipitate)
 {
 
 	if (home->cell2QentArray == (C2Qent_t *)NULL) {
 		return;
 	}
 
         if (precipitate == (Precipitate_t *)NULL) {
 		return;
 	}
 
         if ((precipitate->cell2Idx < 0) || (precipitate->cell2QentIdx < 0)) {
 		return;
 	}
 
         home->cell2QentArray[precipitate->cell2QentIdx].precipitate = (Precipitate_t *)NULL;
         home->cell2QentArray[precipitate->cell2QentIdx].next = -1;
 
         return;
 }