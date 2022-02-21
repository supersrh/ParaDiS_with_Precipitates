/*****************************************************************************
 *
 *  QueueOps.h   Define the Queue operation prototypes
 *
 ****************************************************************************/

#ifndef _QueueOps_h
#define _QueueOps_h

Node_t *PopFreeNodeQ(Home_t *home);
void PushFreeNodeQ(Home_t *home, Node_t *node);
void PushNativeNodeQ(Home_t *home, Node_t *node);
void PushGhostNodeQ(Home_t *home, Node_t *node);
void RecycleGhostNodes(Home_t *home);
void RemoveNodeFromCellQ(Home_t *home, Node_t *node);
void RemoveNodeFromCell2Q(Home_t *home, Node_t *node);
Precipitate_t *PopFreePrecipitateQ(Home_t *home);
void PushFreePrecipitateQ(Home_t *home, Precipitate_t *precipitate);
void PushNativePrecipitateQ(Home_t *home, Precipitate_t *precipitate);
void PushGhostPrecipitateQ(Home_t *home, Precipitate_t *precipitate);
void RecycleGhostPrecipitates(Home_t *home);
void RemovePrecipitateFromCellQ(Home_t *home, Precipitate_t *precipitate);
void RemovePrecipitateFromCell2Q(Home_t *home, Precipitate_t *precipitate);

#endif
