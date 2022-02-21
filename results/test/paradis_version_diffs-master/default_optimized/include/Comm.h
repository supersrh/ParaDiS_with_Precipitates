13,14c13,14
< #define FLTS_PER_GHOST_NODE      21
< #define FLTS_PER_GHOST2_NODE     23
---
> #define FLTS_PER_GHOST_NODE      25
> #define FLTS_PER_GHOST2_NODE     27
16,17c16,19
< #define FLTS_PER_GHOST_NODE      16
< #define FLTS_PER_GHOST2_NODE     18
---
> #define FLTS_PER_GHOST_NODE      20
> #define FLTS_PER_GHOST2_NODE     22
> #define FLTS_PER_GHOST_PRECIPITATE      18
> #define FLTS_PER_GHOST2_PRECIPITATE     20
19a22
> #define FLTS_PER_GHOST_PNRB       2
23a27
> 
25c29,30
< #define INIT_VALS_PER_NODE  7 
---
> #define INIT_VALS_PER_NODE  11
> #define INIT_VALS_PER_PRECIPITATE  8  
26a32
> #define INIT_VALS_PER_PNRB   2 
34c40
< #define INTS_PER_MIRROR_NODE 4
---
> #define INTS_PER_MIRROR_NODE 5
44a51,52
> #define MSG_INIT_PLENS     1003
> #define MSG_INIT_PRECIPITATES    1004
47,48c55,60
< #define MSG_MIG_LEN       1040
< #define MSG_MIG           1050
---
> #define MSG_GHOST_PR      1035
> #define MSG_GHOST_PR_LEN  1037
> #define MSG_MIG_LENGTH 3000
> #define MSG_MIGP_LENGTH 3001
> #define MSG_MIG_NODES  3002
> #define MSG_MIG_PRECIPITATES  3003
68a81,82
> #define	MSG_PR_TAGREMAP_LEN  3060
> #define	MSG_PR_TAGREMAP	  3070
73a88
> void CommSendPrecipitateGhosts(Home_t *home);
76a92
> void CommSendSecondaryPrecipitateGhosts(Home_t *home);
