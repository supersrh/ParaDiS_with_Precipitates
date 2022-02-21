76a77,78
> void   AllocNodePrecipitates(Node_t *node, int n);
> void   AllocPrecipitateArms(Precipitate_t *precipitate, int n);
77a80
> void   FreePrecipitate(Home_t *home, int index);
78a82
> void   FreeNodePNbrs(Node_t *node);
81a86,88
> void PrintAllprecipitateNBR(Home_t *home);
> Precipitate_t *GetPrecipitateFromTag(Home_t *home, Tag_t tag);
> Precipitate_t *GetPrecipitateFromIndex(Home_t *home, int domID, int index);
86a94,96
> void   PrintNodelist(Home_t *home);
> int NodeinPrecipitateforNbrList(real8 R,real8 xp,real8 yp,real8 zp, real8 x1,real8 y1,real8 z1);
> int   NodeinPrecipitate(real8 R,real8 xp,real8 yp,real8 zp, real8 x1,real8 y1,real8 z1,real8 x2,real8 y2,real8 z2);
100a111,115
> int    GetFreePrecipitateTag(Home_t *home);
> int    GetRecycledPrecipitateTag(Home_t *home);
> void   RecyclePrecipitateTag(Home_t *home, int tagIdx);
> 
> 
145a161
> void   SortNativePrecipitates(Home_t *home);
