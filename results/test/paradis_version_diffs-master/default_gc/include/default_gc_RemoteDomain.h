23,24c23,29
< 				/* to Node_t                           */
< 
---
> 				/* to Node_t  
> 				 *
> 				                          */
> 	int	maxprecipitateTagIndex;			                          
> 	Precipitate_t	**precipitateKeys;
> 	/* indexed by precipitates's tag.index, points */
> 				/* to Precipitate_t */
25a31
> 	int	inPBufLen;
26a33,34
> 	char	*inPBuf;
> 	
27a36
> 	int	outPBufLen;
28a38
> 	char	*outPBuf;
