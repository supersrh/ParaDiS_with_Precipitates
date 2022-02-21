10c10,16
<           int numChains, int seed, real8 *totDislocLen, int dislocType);
---
>           int numChains, int seed, real8 *totDislocLen, int dislocType,int nodezero);
> void  CreateSingleEdges(Home_t *home, InData_t *inData, int cubeLength,
>           int numChains, int seed, real8 *totDislocLen, int dislocType,int nodezero);
> void  CreateD2Edges(Home_t *home, InData_t *inData, int cubeLength,
>           int numChains, int seed, real8 *totDislocLen, int dislocType,int nodezero);                    
> void  CreatePrecipitates(Home_t *home, InPrecipitateData_t *inprecipitateData, int cubeLength,
>                  int numPrecipitates, int seed, real8 r,real8 forcep,int nodezero);          
12c18
<           int numChains, int seed, real8 *totDislocLen, int dislocType);
---
>           int numChains, int seed, real8 *totDislocLen, int dislocType,int nodezero);
14c20
<           int numChains, int seed, real8 *totDislocLen, int dislocType);
---
>           int numChains, int seed, real8 *totDislocLen, int dislocType,int nodezero);
17c23
<           int seed, real8 *totDislocLen, int dislocType);
---
>           int seed, real8 *totDislocLen, int dislocType,int nodezero);
20c26
<           real8 *totDislocLen, int dislocType);
---
>           real8 *totDislocLen, int dislocType,int nodezero);
22c28,32
<           int numChains, int seed, real8 *totDislocLen, int dislocType);
---
>           int numChains, int seed, real8 *totDislocLen, int dislocType,int nodezero);
>           
>  void  CreateFCC2DConfig(Home_t *home, InData_t *inData, int cubeLength,
>           int numChains, int seed, real8 *totDislocLen, int dislocType,int nodezero);         
>           
25c35
<           real8 *totDislocLen, int dislocType);
---
>           real8 *totDislocLen, int dislocType,int nodezero);
27c37
<           int numChains, int seed, real8 *totDislocLen, int dislocType);
---
>           int numChains, int seed, real8 *totDislocLen, int dislocType,int nodezero);
30c40
< void  WriteInitialNodeData(Home_t *home, InData_t *inData, int lastBlock);
---
> void  WriteInitialNodeData(Home_t *home, InData_t *inData, int lastBlock,int nodezero);
31a42
> void WriteInitialPrecipitateData(Home_t *home, InPrecipitateData_t *inprecipitateData, int lastBlock, int precipitatezero);
45c56,61
< #define FTYPE_MAX		8
---
> #define FTYPE_PRECIPITATES		8
> #define FTYPE_FCC2D		9
> #define FTYPE_SEDGE		10
> #define FTYPE_D2EDGE		11
> #define FTYPE_MAX		12
> 
48a65,66
> #define FNAME_SEDGE		"sedge"
> #define FNAME_D2EDGE		"d2edge"
54a73,74
> #define FNAME_PRECIPITATES	"precipitates"
> #define FNAME_FCC2D		"fcc2d"
94c114,120
< #define OPT_MAX		18
---
> #define OPT_FORCEP	18
> #define OPT_NPRECIPITATES 19
> #define OPT_NODEZERO 20
> #define OPT_MAX		21
> 
> 
> 
135a162,165
>         real8	forcep;
>         int		numPrecipitates; 
>         int		nodezero; 
>               
