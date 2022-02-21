65c65
<         int     i, newNodeKeyPtr;
---
>         int     i, newNodeKeyPtr,newPrecipitateKeyPtr;
68,69c68,69
<         char    fileName[256], ctrlFile[256], dataFile[256];
<         FILE    *fp, *fpCtrl;
---
>         char    fileName[256], ctrlFile[256], dataFile[256], precipitatedataFile[256];
>         FILE    *fp, *fpCtrl,*fprecip;
70a71
>         Precipitate_t  *precipitate;
114a116,119
>                      
>              snprintf(precipitatedataFile, sizeof(precipitatedataFile), "%s/%s%s",
>                      DIR_RESTART, fileName, PRECIPITATEDATA_FILE_SUFFIX);        
>                      
117a123,126
>                      
>             snprintf(precipitatedataFile, sizeof(precipitatedataFile), "%s/%s%s",
>                      DIR_RESTART, fileName, PRECIPITATEDATA_FILE_SUFFIX,ioGroup);          
>                      
145a155,160
>             
>             if ((fprecip = fopen(precipitatedataFile, "w")) == (FILE *)NULL) {
>                 Fatal("WriteCN: Open error %d on %s\n", errno, dataFile);
>             }
>             
>             
187a203,227
>                              
>                              
>                        
>          /*Write the data file parameters
>  */
>                 WriteParam(home->precipitateParamList, -1, fprecip);
> 
> /*
>  *              Write the domain decomposition into nodal data file
>  *              and then some comment lines describing the nodal
>  *              data that will follow.
>  */
>                 fprintf(fprecip, "\n#\n#  END OF DATA FILE PARAMETERS\n#\n\n");
>                 fprintf(fprecip, "domainDecomposition = \n");
>                 WriteDecompBounds(home, fprecip);
> 
>                 fprintf(fprecip, "precipitateData = \n");
>                 fprintf(fprecip, "#  Primary lines: precipitate_tag, x, y, z, "
>                         "forcep, r\n");
>                 fprintf(fprecip, "#  Secondary lines: arm_tag, burgx, burgy, "
>                              "burgz, nx, ny, nz\n");
>                     
>                              
>                          
>                              
197a238,244
>             
>             if ((fprecip = fopen(precipitatedataFile, "a")) == (FILE *)NULL) {
>                 Fatal("WriteCN: Open error %d on %s\n", errno, dataFile);
>             }
>             
>             
>             
201c248,249
<  *      Now dump the data for all nodes in this block.
---
>  *      Now dump the data for all nodes
>  *  in this block.
243a292,331
>        
>         
>           
>         
>         
>         
>          newPrecipitateKeyPtr = home->newPrecipitateKeyPtr;
>         
>         for (i = 0; i < newPrecipitateKeyPtr; i++) {
> 			
>             if ((precipitate = home->precipitateKeys[i]) == (Precipitate_t *)NULL) {
> 				
>                 continue;
>                  
>             }
> 
> /*
>  *          For now, add a temporary sanity check. If we find any 
>  *          unconstrained node (i.e. not a pinned node, surface
>  *          node, etc) that is singly-linked we've got a problem,
>  *          so abort without writing the restart file.
>  *
>  *          Once we identify the problem in the code that is leading to
>  *          the singly-linked nodes, we get get rid of this check.
>  */
>             
> 			//printf("Precipitate %d,%d %.8f %.8f %.8f %.8f  %.8f", precipitate->myTag.domainID, precipitate->myTag.index,
>                     //precipitate->x, precipitate->y, precipitate->z, precipitate->forcep,
>                     //precipitate->r);
>             fprintf(fprecip,
>                     " %d,%d %.8f %.8f %.8f %.8f  %.8f \n",
>                     precipitate->myTag.domainID, precipitate->myTag.index,
>                     precipitate->x, precipitate->y, precipitate->z, precipitate->forcep,
>                     precipitate->r);
>        
>     
>            
>         }
>         
>         fclose(fprecip);
244a333,335
>         
>         
>         
