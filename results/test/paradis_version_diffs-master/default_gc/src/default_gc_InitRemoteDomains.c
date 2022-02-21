90a91,94
>                     remDom->maxprecipitateTagIndex = 0;
>                     remDom->inBufLen=0;
>                     remDom->inPBufLen=0;
>                     
214a219,231
>         home->inPRequests  = (MPI_Request *) malloc(home->numDomains *
>                                                   sizeof(MPI_Request));
>         home->outPRequests = (MPI_Request *) malloc(home->numDomains *
>                                                   sizeof(MPI_Request));
>         home->inPStatus    = (MPI_Status *) malloc(home->numDomains *
>                                                   sizeof(MPI_Status));
>         home->outPStatus   = (MPI_Status *) malloc(home->numDomains *
>                                                   sizeof(MPI_Status));                                          
>                                                   
>                                                   
>                                                   
>                                                   
>                                                   
