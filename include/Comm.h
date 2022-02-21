/*************************************************************************
 *
 *  Comm.h - define various constants used during MPI communications
 *
 *************************************************************************/

#ifndef _Comm_h
#define _Comm_h

#include "Home.h"

#ifdef _FEM
#define FLTS_PER_GHOST_NODE      25
#define FLTS_PER_GHOST2_NODE     27
#else
#define FLTS_PER_GHOST_NODE      20
#define FLTS_PER_GHOST2_NODE     22
#define FLTS_PER_GHOST_PRECIPITATE      18
#define FLTS_PER_GHOST2_PRECIPITATE     20
#endif

#define FLTS_PER_GHOST_PNRB       2
#define FLTS_PER_GHOST_ARM       11
#define FLTS_PER_GHOST_CELL       3
#define EXTRA_GHOST_FLTS          4

#define INIT_PARAMS         1
#define INIT_VALS_PER_NODE  11
#define INIT_VALS_PER_PRECIPITATE  8  
#define INIT_VALS_PER_ARM   8 
#define INIT_VALS_PER_PNRB   2 

#define INTS_PER_TAG        2
/*
 *	For downloading entire problem space to domain zero.
 *	Used for various plotting functions and writing text
 * 	restart files.
 */
#define INTS_PER_MIRROR_NODE 5
#define INTS_PER_MIRROR_ARM  2
#define FLTS_PER_MIRROR_NODE 6
#define FLTS_PER_MIRROR_ARM  9
#define EXTRA_MIRROR_INTS    6

/*
 *	message tags
 */
#define MSG_INIT_LENS     1000
#define MSG_INIT_NODES    1002
#define MSG_INIT_PLENS     1003
#define MSG_INIT_PRECIPITATES    1004
#define MSG_GHOST_LEN     1020
#define MSG_GHOST         1030
#define MSG_GHOST_PR      1035
#define MSG_GHOST_PR_LEN  1037
#define MSG_MIG_LENGTH 3000
#define MSG_MIGP_LENGTH 3001
#define MSG_MIG_NODES  3002
#define MSG_MIG_PRECIPITATES  3003
#define MSG_OLDNEW_LEN    1060
#define MSG_OLDNEW        1070
#define MSG_MIRRORS       1090
#define	MSG_SEND_MIRROR   1091
#define MSG_GHOST2_REQ_LEN      1100
#define MSG_GHOST2_REQ          1101
#define MSG_GHOST2_RESPONSE_LEN 1102
#define MSG_GHOST2_RESPONSE     1103
#define MSG_VELOCITY      2000
#define MSG_VELOCITY_LEN  2001
#define MSG_OPLIST_LEN    2010
#define MSG_OPLIST        2020
#define MSG_TOKEN_RING    2030
#define MSG_REMESH_LEN    2040
#define MSG_REMESH        2050
#define	MSG_TAGREMAP_LEN  2060
#define	MSG_TAGREMAP	  2070
#define MSG_SEGDATA_LEN   2080
#define MSG_SEGDATA       2090
#define MSG_VISIT_COUNTS  2100
#define	MSG_PR_TAGREMAP_LEN  3060
#define	MSG_PR_TAGREMAP	  3070
#define MSG_GHOST_PR      1035
#define MSG_GHOST_PR_LEN  1037

/*
 *	Prototypes
 */
void CommSendGhosts(Home_t *home);
void CommSendPrecipitateGhosts(Home_t *home);
void CommSendGhostPlanes(Home_t *home);
void CommSendRemesh(Home_t *home);
void CommSendSecondaryGhosts(Home_t *home);
void CommSendSegments(Home_t *home, int numRecvBufs, int numSendBufs,
         int *sendDomList, Segment_t **cellSegLists, int *cellSegCnts);
void CommSendVelocity(Home_t *home);
void CommSendSecondaryPrecipitateGhosts(Home_t *home);

#endif
