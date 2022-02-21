342a343
> 		intBuf[iIdx++] = node->numPNbrs;
377a379,381
> 		fltBuf[fIdx++] = node->nodedx;
> 		fltBuf[fIdx++] = node->nodedy;
> 		fltBuf[fIdx++] = node->nodedz;
423c427
< 	int		intCount, fltCount, i, inode, numNbrs, armCount;
---
> 	int		intCount, fltCount, i, inode, numNbrs,numPNbrs, armCount;
485a490
> 		numPNbrs                  = intBuf[iIdx++];
487a493
> 		AllocNodePrecipitates(node,numPNbrs );
