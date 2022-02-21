#include "Home.h"

/*---------------------------------------------------------------------------
 *
 *      Function:     DisorderForce by AL
 *      Description:  This routine calculates the gaussian  from  precipitate
 * 						 
 *                    
 *
 *                    
 *
 *      Arguments:
 *          x1,y1,z1  Coordinates of segment's first endpoint
 *          x2,y2,z2  Coordinates of segment's second endpoint
 *          bx,by,bz  Components of burgers vector from first endpoint
 *                    to second endpoint
 *          f1        Array in which to return osmotic force at pos1
 *          f2        Array in which to return osmotic force at pos2
 *
 *      Author:       Arttu Lehtinen
 *
 *-------------------------------------------------------------------------*/
void DisorderForce(Home_t *home,Node_t *node1, Node_t *node2, 
                  real8 f1[3], real8 f2[3])
{
		real8	 x1,y1,z1,x2,y2,z2,bx,by,bz;	
        real8    eps, gaussForce1,gaussForce2, bDotrt,r,forcep;
        real8   xCellSize, yCellSize, zCellSize, minCellSize;
        real8    rtNorm,rn1Norm,rn2Norm,r1Norm,r2Norm,rpNorm,checkScrewNorm, bEdgeNorm,r1pangle,r2pangle,r12angle;
        real8    temperature, atomicVol, cv, cvFloor, cvEquilibrium,vproj;
        real8    burg[3],xs,ys,zs,a,b,c,t1,t2;
        real8    r0[3],rt[3],rn1[3],rn2[3],r1[3],r2[3],rp[3],rtUnit[3], burgUnit[3], checkScrew[3],velocity1[3],velocity2[3];
        real8	 rt1proj,rt2proj;
        real8    bEdge[3], tmp3[3],bEdgeunit[3],rm[3];
        int 	 newPrecipitateKeyPtr,i,in,nop,nop1,nop2,noptot,j,truth;	
        Precipitate_t *precipitate;
        Param_t *param;
        Tag_t	tag,tagtemp,*tags;
		in=0;
        eps = 1.0e-12;
		t1=-1.0;
		t2=-1.0;
        param = home->param;
        cv = param->vacancyConc;
        cvEquilibrium = param->vacancyConcEquilibrium;
        temperature = param->TempK;
        atomicVol = param->burgMag * param->burgMag * param->burgMag;
		
        
        
        xCellSize = param->Lx / param->nXcells;
        yCellSize = param->Ly / param->nYcells;
        zCellSize = param->Lz / param->nZcells;

        minCellSize = MIN(xCellSize, yCellSize);
        minCellSize = MIN(minCellSize, zCellSize);
        
        f1[0] = 0.0;
	f1[1] = 0.0;
        f1[2] = 0.0;

        f2[0] = 0.0;
        f2[1] = 0.0;
        f2[2] = 0.0;
        
        
        
        
        x1=node1->x;
        y1=node1->y;
        z1=node1->z;
        
        x2=node2->x;
        y2=node2->y;
        z2=node2->z;
        
        noptot=node1->numPNbrs+node2->numPNbrs;
        nop1=node1->numPNbrs;
        nop2=node2->numPNbrs;
        
        if(noptot==0){
		f1[0] =f1[0]+0.0; 
        f1[1] =f1[1]+0.0; 
        f1[2] =f1[2]+0.0; 
        
        f2[0] =f2[0]+0.0; 
        f2[1] =f2[1]+0.0; 
        f2[2] =f2[2]+0.0;	
			
			
        return;   
	    }
        
        //for (j =0; j < nop; j++) {
        
	//}	
	//allocating a joint neigborlist
	tags = (Tag_t *)calloc(noptot,sizeof(Tag_t) );
	
	
	
	if((nop1>0) && (nop2 == 0)){
		for (i = 0; i < nop1; i++) {
			tags[i]=node1->PnbrTag[i];		
		}
		nop=nop1;	
	}

	if((nop1 == 0) && (nop2 > 0)){
		for (i = 0; i < nop2; i++) {
			tags[i]=node2->PnbrTag[i];		
		}
		nop=nop2;	
	}
	
	
	
	if(nop1 > 0 && nop2 > 0){
		
		for (i = 0; i < nop1; i++) {
			tags[i]=node1->PnbrTag[i];		
		}
		
	//Finding the diffrence between neighborlists
		nop=nop1;
		for (i = 0; i < nop2; i++) {
			truth=-1;
			for (j = 0; j < nop1; j++) {
				if(OrderTags(&node1->PnbrTag[j], &node2->PnbrTag[i])==0){
					truth=1;
				//printf("Samat \n");
				}
			
			}
			if(truth<0){
			tags[nop]=node2->PnbrTag[i];
			nop=nop+1;
			}
		}	
	}	
		
		//for (i = 0; i < nop; i++) {
			//printf("Disorder jointPNBRlist (%d %d) \n",tags[i].domainID,tags[i].index);
	
		//}
	
        for (i = 0; i < nop; i++) {
			
			//if ((precipitate = home->precipitateKeys[i]) == (Precipitate_t *)NULL) {
                //continue;
            //}

		tag=tags[i];
		
		
        precipitate=GetPrecipitateFromTag(home,tag);
       
       if(precipitate->constraint>0){
		  continue; 
		   }
		
		r=precipitate->r;
		r0[0]=precipitate->x;
		r0[1]=precipitate->y;
		r0[2]=precipitate->z;

		
		r1[0] =x1-r0[0];
        r1[1] = y1-r0[1];
        r1[2] = z1-r0[2];
        r1Norm=Normal(r1);
        
        r2[0] =x2-r0[0];
        r2[1] = y2-r0[1];
        r2[2] = z2-r0[2];
        r2Norm=Normal(r2);
        
        //korjaa sellaiseksi että katsotaan segmentin keskelle!
     
		 /* periodic boundaries*/
		// printf("Distances node1 %e  node 2 %e \n",r1Norm,r2Norm);
		 if((abs(r1[0])>minCellSize) || (abs(r1[1])>minCellSize) || (abs(r1[2])>minCellSize)){ 
			ZImage(param, &r1[0], &r1[1], &r1[2]);
			 r1Norm=Normal(r1);						
		 }
		 
		 
		 if((abs(r2[0])>minCellSize) || (abs(r2[1])>minCellSize) || (abs(r2[2])>minCellSize)){	
			ZImage(param, &r2[0], &r2[1], &r2[2]);
			 r2Norm=Normal(r2);	
                	
		 }	 
		 
		rt[0] = x2 - x1;
        rt[1] = y2 - y1;
        rt[2] = z2 - z1;
        /*PBC*/
        if( (rt[0]>minCellSize) || (rt[1]>minCellSize) || (rt[2]>minCellSize)){
					ZImage(param, &rt[0], &rt[1], &rt[2]);
					
						}	 
	rm[0]=x1+0.5*rt[0];
	rm[1]=y1+0.5*rt[1];
	rm[2]=z1+0.5*rt[2];	 
		 	 

in=1;		
if(in>0){

			
	/*
 *      Get directional vector of the dislocation and the corresponding
 *      unit vector
 */		
		
		       
		rtNorm = Normal(rt);
		        
        rtUnit[0] = rt[0] / rtNorm;
        rtUnit[1] = rt[1] / rtNorm;
        rtUnit[2] = rt[2] / rtNorm;
        
        
     
       
 	/*We only need the component which is normal to the dislocation segment (AL) */	      
       
        rt1proj=DotProduct(r1, rt);
        rn1[0]=r1[0]-rt1proj*rt[0]/(rtNorm * rtNorm);
		rn1[1]=r1[1]-rt1proj*rt[1]/(rtNorm * rtNorm);
		rn1[2]=r1[2]-rt1proj*rt[2]/(rtNorm * rtNorm);
		rn1Norm = Normal(rn1);
       
        rt2proj=DotProduct(r2, rt);
        rn2[0]=r2[0]-rt2proj*rt[0]/(rtNorm * rtNorm);
		rn2[1]=r2[1]-rt2proj*rt[1]/(rtNorm * rtNorm);
		rn2[2]=r2[2]-rt2proj*rt[2]/(rtNorm * rtNorm);
		rn2Norm = Normal(rn2);
       
            forcep=precipitate->forcep;        
            gaussForce1=-2.0*(forcep)*r1Norm*exp(-r1Norm*r1Norm/(r*r))/(r*r)*rtNorm ;
			gaussForce2=-2.0*(forcep)*r2Norm*exp(-r2Norm*r2Norm/(r*r))/(r*r)*rtNorm ;


/*     Force per unit length of the segment at each endpoint...
 */
 
		f1[0] =f1[0] +gaussForce1 * -rn1[0]/rn1Norm; 
        f1[1] =f1[1] + gaussForce1 * -rn1[1]/rn1Norm; 
        f1[2] =f1[2] + gaussForce1 * -rn1[2]/rn1Norm; 
        
        f2[0] =f2[0]+gaussForce2 * -rn2[0]/rn2Norm; 
        f2[1] =f2[1]+ gaussForce2 * -rn2[1]/rn2Norm; 
        f2[2] =f2[2]+gaussForce2 * -rn2[2]/rn2Norm;
                
 //printf( "%d %d \n",precipitate->myTag.index,NodeinPrecipitate(r,r0[0],r0[1],r0[2], x1, y1,z1,x2,y2,z2));       
         
        
       //printf("f1 %e %e %e %e %d \n",f1[0],f1[1],f1[2],rtNorm,in);
       //printf("r1 %e %e %e  \n",r1[0],r1[1],r1[2]); 
        
        
        
        
	}
        
	
	
	
	 
	else{
		f1[0] =f1[0]+0.0; 
        f1[1] =f1[1]+0.0; 
        f1[2] =f1[2]+0.0; 
        
        f2[0] =f2[0]+0.0; 
        f2[1] =f2[1]+0.0; 
        f2[2] =f2[2]+0.0;
        	
	 }
		
			
		
	}	
	
	free(tags);
        return;
}



		




		
