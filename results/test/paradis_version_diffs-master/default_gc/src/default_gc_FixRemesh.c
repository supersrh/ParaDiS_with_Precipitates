25c25
<         real8           cellSize[3], cellIndex[3];
---
>     real8           cellSize[3], cellIndex[3];
110c110
< 				
---
> 				FreeNodePNbrs(node2);
120c120,123
< 
---
> 				
> 				node2->nodedx = op->nodedx;
> 				node2->nodedy = op->nodedy;
> 				node2->nodedz = op->nodedz;
228a232,237
> 				
> 				node1->nodedx = op->nodedx;
> 				node1->nodedy = op->nodedy;
> 				node1->nodedz = op->nodedz;
> 
> 
