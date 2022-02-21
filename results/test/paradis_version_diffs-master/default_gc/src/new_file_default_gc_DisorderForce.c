#include "Home.h"

/*---------------------------------------------------------------------------
 *
 *      Function:     DisorderForce by AL
 *      Description:  This routine calculates the force  from generic precipitate
 * 						 
 *                    
 *
 *                    The units of the parameters of the force are tied to
 *                    the other units of the simulations but typical units
 *                    are:
 *                        Boltmann's constant -> Joules/Kelvin
 *                        Atomic Volume - > in cubic burgers vectors
 *                        cvEquilibrium - > number density of vacancies in
 *                                          thermal equilibrium per atomic
 *                                          volume
 *                        cv - > number density of vacancies per atomic volume
 *
 *                    When this routine is used, there is a balance between
 *                    the climb of dislocations changes and the vacancy
 *                    concentation in the simulation volume.
 *
 *      Arguments:
 *          x1,y1,z1  Coordinates of segment's first endpoint
 *          x2,y2,z2  Coordinates of segment's second endpoint
 *          bx,by,bz  Components of burgers vector from first endpoint
 *                    to second endpoint
 *          f1        Array in which to return osmotic force at pos1
 *          f2        Array in which to return osmotic force at pos2
 *
 *      Author:       Tuan Hoang
 *
 *-------------------------------------------------------------------------*/
void DisorderForce(Home_t *home,Node_t *node1, Node_t *node2, 
                  real8 f1[3], real8 f2[3])
{
		real8	 x1,y1,z1,x2,y2,z2,bx,by,bz;	
        real8    eps, osForce1,osForce2, bDotrt,r,forcep;
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
		//return;
		
        //PrintNode(node1);
        //PrintNode(node2);
        
        
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
		 	 
//in=NodeinPrecipitateforNbrList(10.0*r, r0[0],r0[1],r0[2],rm[0],rm[1],rm[2]);		 
//in=NodeinPrecipitate(10.0*r,r0[0],r0[1],r0[2], x1, y1,z1,x2,y2,z2);
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
        
        
     
       
 	/*We only need the komponent which is normal to the dislocation segment (AL) */	      
       
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
            osForce1=-2.0*(forcep)*r1Norm*exp(-r1Norm*r1Norm/(r*r))/(r*r)*rtNorm ;
			osForce2=-2.0*(forcep)*r2Norm*exp(-r2Norm*r2Norm/(r*r))/(r*r)*rtNorm ;


/*     Osmotic force per unit length of the segment at each endpoint...
 */
 
		f1[0] =f1[0] +osForce1 * -rn1[0]/rn1Norm; 
        f1[1] =f1[1] + osForce1 * -rn1[1]/rn1Norm; 
        f1[2] =f1[2] + osForce1 * -rn1[2]/rn1Norm; 
        
        f2[0] =f2[0]+osForce2 * -rn2[0]/rn2Norm; 
        f2[1] =f2[1]+ osForce2 * -rn2[1]/rn2Norm; 
        f2[2] =f2[2]+osForce2 * -rn2[2]/rn2Norm;
                
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


/*raytracin hard sphere version */
//void DisorderForce(Home_t *home, real8 x1, real8 y1, real8 z1,
                  //real8 x2, real8 y2, real8 z2,real8 vx1,real8 vy1,real8 vz1,real8 vx2,real8 vy2, real8 vz2, real8 bx, real8 by, real8 bz,
                  //real8 f1[3], real8 f2[3])
//{
        //real8    eps, osForce, bDotrt,r;
        //real8    rtNorm,rn1Norm,rn2Norm,r1Norm,r2Norm,rpNorm,checkScrewNorm, bEdgeNorm,r1pangle,r2pangle,r12angle;
        //real8    temperature, atomicVol, cv, cvFloor, cvEquilibrium,vproj;
        //real8    burg[3],xs,ys,zs,a,b,c,t1,t2;
        //real8    r0[3],rt[3],rn1[3],rn2[3],r1[3],r2[3],rp[3],rtUnit[3], burgUnit[3], checkScrew[3],velocity1[3],velocity2[3];
        //real8	 rt1proj,rt2proj;
        //real8    bEdge[3], tmp3[3],bEdgeunit[3];
        //Param_t *param;
		//int 	 newPrecipitateKeyPtr,i;	
        //Precipitate_t *precipitate;
        //eps = 1.0e-12;
		//t1=-1.0;
		//t2=-1.0;
        //param = home->param;
        //cv = param->vacancyConc;
        //cvEquilibrium = param->vacancyConcEquilibrium;
        //temperature = param->TempK;
        //atomicVol = param->burgMag * param->burgMag * param->burgMag;

//newPrecipitateKeyPtr= home->newPrecipitateKeyPtr;
        
        //for (i = 0; i < newPrecipitateKeyPtr; i++) {
			
			//if ((precipitate = home->precipitateKeys[i]) == (Precipitate_t *)NULL) {
                //continue;
            //}



        //burg[0] = bx;
        //burg[1] = by;
        //burg[2] = bz;


		//r=precipitate->r;
		//xs=precipitate->x;
		//ys=precipitate->y;
		//zs=precipitate->z;
		
		//a=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
		//b=2.0*((x1-xs)*(x2-x1)+(y1-ys)*(y2-y1)+(z1-zs)*(z2-z1));
		//c=(x1-xs)*(x1-xs)+(y1-ys)*(y1-ys)+(z1-zs)*(z1-zs)-r*r;
		
		//if((b*b-4.0*a*c)>=0.0){
		
		//t1=(-b+sqrt(b*b-4.0*a*c))/(2.0*a);
		//t2=(-b-sqrt(b*b-4.0*a*c))/(2.0*a);
		
	
	
	//if((t1>=0.0) && (t1<=1.0)){
	//printf("t1 %e  \n",t1);
		//rp[0] =x1+t1*(x2-x1);
        //rp[1] =y1+t1*(y2-y1);
        //rp[2] =z1+t1*(z2-z1);
		//rpNorm = Normal(rp);
	
	
//}

//if((t2>=0.0) && (t2<=1.0)){
	//printf("t2 %e  \n",t2);
	
		//rp[0] =x1+t2*(x2-x1);
        //rp[1] =y1+t2*(y2-y1);
        //rp[2] =z1+t2*(z2-z1);
		//rpNorm = Normal(rp);
	
//}	

			
	//*
 //*      Get directional vector of the dislocation and the corresponding
 //*      unit vector
 //*/		
		
		//r0[0]=precipitate->x;
		//r0[1]=precipitate->y;
		//r0[2]=precipitate->z;
 
        //rt[0] = x2 - x1;
        //rt[1] = y2 - y1;
        //rt[2] = z2 - z1;
		//rtNorm = Normal(rt);
        
        //rtUnit[0] = rt[0] / rtNorm;
        //rtUnit[1] = rt[1] / rtNorm;
        //rtUnit[2] = rt[2] / rtNorm;
        
        
        //r1[0] =rp[0]-r0[0];
        //r1[1] = rp[1]-r0[1];
        //r1[2] = rp[2]-r0[2];
        //r1Norm=Normal(r1);
       
 	/*We only need the komponent which is normal to the dislocation segment (AL) */	      
       
        //rt1proj=DotProduct(r1, rt);
        //rn1[0]=r1[0]-rt1proj*rt[0]/(rtNorm * rtNorm);
		//rn1[1]=r1[1]-rt1proj*rt[1]/(rtNorm * rtNorm);
		//rn1[2]=r1[2]-rt1proj*rt[2]/(rtNorm * rtNorm);
		//rn1Norm = Normal(rn1);
       
			////printf("forcep %e ",precipitate->forcep);				
            //osForce=precipitate->forcep;/*r*exp(r*r);//(r*r+0.01);//*rtNorm*0.5;
            ////;          
      


 /*      Osmotic force per unit length of the segment at each endpoint...
 //*/
 //if(((t1>=0.0) && (t1<=1.0))|| ( (t2>=0.0) && (t2<=1.0)   ) ){
		//f1[0] =f1[0] +osForce * rn1[0]/rn1Norm; 
        //f1[1] = f1[1]+osForce * rn1[1]/rn1Norm; 
        //f1[2] =f1[2] +osForce * rn1[2]/rn1Norm; 
        
        //f2[0] =f2[0] +osForce * rn1[0]/rn1Norm; 
        //f2[1] = f2[1]+osForce * rn1[1]/rn1Norm; 
        //f2[2] =f2[2] +osForce * rn1[2]/rn1Norm; 
        
        /*printf("rp %e %e %e  \n",rp[0],rp[1],rp[2]);*/
	//}
        
	//} 
	//else{
		//f1[0] =f1[0]+0.0; 
        //f1[1] =f1[1]+0.0; 
        //f1[2] =f1[2]+0.0; 
        
        //f2[0] =f2[0]+0.0; 
        //f2[1] =f2[1]+0.0; 
        //f2[2] =f2[2]+0.0;
	 //}
	//}	        
        //return;
//} 
		




		
