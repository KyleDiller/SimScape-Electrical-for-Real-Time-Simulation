#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "MatrixFun.h"
#define MAX_NB_SWITCH 20
int MateDiscretize(double A[], double B1[],double B2[], double C[], double D[], double h, int opt, int holdtype, int NSTA, int NIN, int NOUT)
{

	double a1,a2,a3,a4;
	double *I;
	int not_alloc=-9;
	double *tmp,*tmp1,*denx,*Ad,*Bd1,*Bd2,*Ah,*Ctmp,*Dtmp,*Bf=NULL,*Bz=NULL;
	double fa=(2.0-sqrt(2.0))/2.0;
	double *tmp0_r, *tmp0_i, *tmp1_r, *tmp1_i, *den1x_r, *den1x_i, *den2x, *num1_r, *num1_i;

	int jz,jf,j,i,k;
	denx = 0;
	tmp = 0;
	Ad=0;
	Bd1=0;


	if ((I=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
	ident(I,NSTA);
	if ((Ah=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }


	


	switch (opt) {

	case 5:			/* Radau IIA(5) */
		if ((den1x_r=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((den1x_i=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((den2x=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((num1_r=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((num1_i=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp0_r=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp0_i=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp1_r=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp1_i=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
		if ((Bd1=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
		}
		if ((Bd2=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
		}


        a1=-2.681082873627752133895790743;
        a2=-3.0504301992474105694263776247;
        a3=- 3.6378342527444957322084185135777;

        opmx(Ah,A,&h,NSTA,NSTA,2);          //A*h

        opmx(tmp0_r,I,&a1,NSTA,NSTA,2);        //
        opmx(den1x_r,Ah,tmp0_r,NSTA,NSTA,0);   // 
        opmx(den1x_i,I,&a2,NSTA,NSTA,2);    // 
        matinvcpx(den1x_r,den1x_i,NSTA);    //

        opmx(tmp0_r,I,&a3,NSTA,NSTA,2);        // 
        opmx(den2x,Ah,tmp0_r,NSTA,NSTA,0);     //
        matinv2(den2x,NSTA);                // 

        a1= 4.0;
        a2= 2.0;
        a3=-3.0;

        opmx(tmp0_r,I,&a1,NSTA,NSTA,2);
        opmx(num1_r,Ah,tmp0_r,NSTA,NSTA,0);    // A*h +4*I

        opmx(num1_i,I,&a2,NSTA,NSTA,2);     // 2*i*I

        mulmxcpx(tmp0_r,tmp0_i,den1x_r,den1x_i,0,num1_r,num1_i,1,NSTA,NSTA,NSTA);
        mulmxcpx(tmp1_r,tmp1_i,tmp0_r,tmp0_i,0,den1x_r,den1x_i,1,NSTA,NSTA,NSTA);
        mulmxcpx(tmp0_r,tmp0_i,tmp1_r,tmp1_i,0,num1_r,num1_i,0,NSTA,NSTA,NSTA);
        mulmx(tmp1_r,tmp0_r,den2x,NSTA,NSTA,NSTA);
        opmx(A,tmp1_r,&a3,NSTA,NSTA,2);


		if (holdtype==1){



			a1=-4.0;
			a2=-3.741657386773941385583748732316549301756;
			a3=-1.0*h;

			opmx(tmp0_r,I,&a1,NSTA,NSTA,2);
			opmx(num1_r,Ah,tmp0_r,NSTA,NSTA,0);    // A*h - 4*I

			opmx(num1_i,I,&a2,NSTA,NSTA,2);     // 

			mulmxcpx(tmp0_r,tmp0_i,den1x_r,den1x_i,0,num1_r,num1_i,0,NSTA,NSTA,NSTA);
			mulmxcpx(tmp1_r,tmp1_i,tmp0_r,tmp0_i,0,den1x_r,den1x_i,1,NSTA,NSTA,NSTA);
			mulmxcpx(tmp0_r,tmp0_i,tmp1_r,tmp1_i,0,num1_r,num1_i,1,NSTA,NSTA,NSTA);
			mulmx(tmp1_r,tmp0_r,den2x,NSTA,NSTA,NSTA);
			opmx(tmp1_r,tmp1_r,&a3,NSTA,NSTA,2);


			mulmx(Bd2,tmp1_r,B1,NSTA,NIN,NSTA);	/* Bd2 not Bd1!   */

			// 1/2*I + 1/30*h*A =
			// 0.03333333333333*(A*h + 15*I)

			a1=15.0;
			a2=-2.0*h;

			opmx(tmp0_r,I,&a1,NSTA,NSTA,2);
			opmx(num1_r,Ah,tmp0_r,NSTA,NSTA,0);    // A*h +15*I
			opmx(num1_i,I,I,NSTA,NSTA,1);    //  0*i*I

			mulmxcpx(tmp0_r,tmp0_i,den1x_r,den1x_i,0,num1_r,num1_i,0,NSTA,NSTA,NSTA);
			mulmxcpx(tmp1_r,tmp1_i,tmp0_r,tmp0_i,0,den1x_r,den1x_i,1,NSTA,NSTA,NSTA);
			mulmx(tmp0_r,tmp1_r,den2x,NSTA,NSTA,NSTA);
			opmx(tmp1_r,tmp0_r,&a2,NSTA,NSTA,2);

			mulmx(Bd1,tmp1_r,B1,NSTA,NIN,NSTA);	/* Bd1 , not Bd2  ! */



			
			cpmx(B1,Bd1,NSTA*NIN);
			cpmx(B2,Bd2,NSTA*NIN);
		
		}
		
		free(tmp0_r);
		free(tmp0_i);
		free(tmp1_r);
		free(tmp1_i);
		free(den1x_r);
		free(den1x_i);
		free(den2x);
		free(num1_r);
		free(num1_i);
		free(Bd1);
		free(Bd2);
		break;


	case 31:
		fa=(2.0+sqrt(2.0))/2.0;
	case 3:
		if ((denx=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp1=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
		if ((Bd1=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
		}
		if ((Bd2=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
		a4=0.0-fa;
		opmx(Ah,A,&h,NSTA,NSTA,2);          //
		opmx(denx,Ah,&a4,NSTA,NSTA,2);     // 
		opmx(denx,I,denx,NSTA,NSTA,0);	// 
		matinv2(denx,NSTA);


		a1=1.0-2.0*fa;
		opmx(tmp,Ah,&a1,NSTA,NSTA,2);       //  
		opmx(tmp,I,tmp,NSTA,NSTA,0);	    // 

		mulmx(tmp1,denx,tmp,NSTA,NSTA,NSTA);// 
		mulmx(A,tmp1,denx,NSTA,NSTA,NSTA);// 

		if (holdtype==1){
			a1=2.0*fa-fa*fa;
			a2=0.0-fa*fa;
			opmx(tmp,Ah,&a2,NSTA,NSTA,2);  // 
			opmx(tmp1,I,&a1,NSTA,NSTA,2);  //   
			opmx(tmp,tmp,tmp1,NSTA,NSTA,0);	//  

			mulmx(tmp1,denx,tmp,NSTA,NSTA,NSTA);
			mulmx(tmp,tmp1,denx,NSTA,NSTA,NSTA);
			opmx(tmp1,tmp,&h,NSTA,NSTA,2);
			//mulmx(Bd1,tmp1,B1,NSTA,NIN,NSTA);	

		    mulmx(Bd1,tmp1,B1,NSTA,NIN,NSTA);	/* Bd1   */

			a3=(1.0-2.0*fa+fa*fa)*h;
			mulmx(tmp,denx,denx,NSTA,NSTA,NSTA);
			opmx(tmp,tmp,&a3,NSTA,NSTA,2);

			//mulmx(Bd2,tmp,B1,NSTA,NIN,NSTA);	

			mulmx(Bd2,tmp,B1,NSTA,NIN,NSTA);	/* Bd2   */

			
			
			cpmx(B1,Bd1,NSTA*NIN);
			cpmx(B2,Bd2,NSTA*NIN);
			
			//cpmx(B1,Bd1,NSTA*NIN);
			//cpmx(B2,Bd2,NSTA*NIN);
		}
		

		free(denx);
		free(tmp);
		free(tmp1);
		free(Bd1);
		free(Bd2);

		break;

	case 2:

		if ((denx=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
		if ((Bd1=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
		}
		if ((Bd2=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
		a1=0.0-0.5;
		a2=0.5;
		opmx(Ah,A,&h,NSTA,NSTA,2);  /* A*h */
		opmx(denx,Ah,&a1,NSTA,NSTA,2); /* -0.5*h*A */
		opmx(denx,I,denx,NSTA,NSTA,0);	/* I-0.5*h*A */
		matinv2(denx,NSTA);
		opmx(tmp,Ah,&a2,NSTA,NSTA,2); /* 0.5*h*A */
		opmx(tmp,I,tmp,NSTA,NSTA,0);	/* I+0.5*h*A */
		mulmx(A,denx,tmp,NSTA,NSTA,NSTA);  /* 	Ad=denx*(I+h*A/2) */


		opmx(tmp,denx,&h,NSTA,NSTA,2); /* denx*h */
		mulmx(Bd1,tmp,B1,NSTA,NIN,NSTA);  /* denx*h*B	 */
		opmx(Bd2,Bd1,&a2,NSTA,NIN,2); /* denx*h*B/2 */
		cpmx(B1,Bd2,NSTA*NIN);
	    cpmx(B2,Bd2,NSTA*NIN); 
	    
		free(denx);
		free(tmp);
		free(Bd1);
		free(Bd2);


		break;
		
	case 1:

		if ((denx=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
		if ((Bd1=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
		}
		if ((Bd2=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
		a1=1.0;
		a2=0.0;
		opmx(Ah,A,&h,NSTA,NSTA,2);  /* A*h */
		opmx(denx,I,Ah,NSTA,NSTA,1);	/* I-h*A */
		matinv2(denx,NSTA); // denx=inv(I-hA)
		cpmx(A,denx,NSTA*NSTA);
	
		opmx(tmp,A,&h,NSTA,NSTA,2); /* denx*h */
		mulmx(Bd2,tmp,B1,NSTA,NIN,NSTA);  /* denx*h*B	 */
		opmx(Bd1,Bd2,&a2,NSTA,NIN,2); /* Bd1=0 */
		cpmx(B1,Bd1,NSTA*NIN);
	    cpmx(B2,Bd2,NSTA*NIN); 
	    
		free(denx);
		free(tmp);
		free(Bd1);
		free(Bd2);


		break;
		

	}



	free(Ah);
	free(I);
	return(0);
}

int MatePermutize(double Amm[], double Bmm[], double Cmm[], double Dmm[], double Rswitch[],int permutation, int NSTATE, int NINPUT, int NOUTPUT, int NSWITCH)
{


// IA=eye(length(Am));
double *IA, *I, Yswitch[MAX_NB_SWITCH],YswitchDiagonal[MAX_NB_SWITCH];
int ETATSW[MAX_NB_SWITCH], NIV_SWITCH[MAX_NB_SWITCH];
int not_alloc=-9;

	

if ((IA=(double*)malloc((NSTATE*NSTATE+1)*sizeof(double)))  == NULL) {
	return(not_alloc);
}
ident(IA,NSTATE);





for (int i=0; i<NSWITCH;i++) ETATSW[i]=-1;
    
// NIV_SWITCH=1:NSWITCH;

for (int i=0; i<NSWITCH;i++) NIV_SWITCH[i]=i;


// Yswitch=1./Rswitch;

for (int i=0; i<NSWITCH;i++) Yswitch[i]=1.0/Rswitch[i];
    
// I=diag(ones(1,NOUTPUT));
if ((I=(double*)malloc((NOUTPUT*NOUTPUT+1)*sizeof(double)))  == NULL) {
	return(not_alloc);
}
ident(I,NOUTPUT);




int idx;

int i;

i=permutation;
    for (int j=0; j<NSWITCH; j++) ETATSW[j]=0;
    int tmp=0;
    int div=1;
    for (int j=0; j<NSWITCH-1; j++) div=div*2;
    idx=i;
    for (int j=0; j<NSWITCH; j++) {
        if ( (idx-div)>=0) {
            ETATSW[j]=1;
            idx=idx-div;
        }
            div=div/2;
    }


    for (int i=0; i<NSWITCH; i++){ 
        if (ETATSW[i]==1) YswitchDiagonal[i]=Yswitch[i];
        else YswitchDiagonal[i]=0;
    }


    double temp, *DxCol, *BdCol, *tmp1, *tmp2;
    if ((DxCol=(double*)malloc((NOUTPUT+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
    if ((BdCol=(double*)malloc((NSTATE+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
    if ((tmp1=(double*)malloc((NSTATE+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
    if ((tmp2=(double*)malloc((NINPUT+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
  
    for(int k=0;k<NSWITCH;k++){


        temp=1.0/(1.0-Dmm[k*NINPUT+k]*YswitchDiagonal[k]);
//  
        for(int m=0;m<NOUTPUT;m++) DxCol[m]=Dmm[m*NINPUT+k]*YswitchDiagonal[k]*temp;
//     
        DxCol[k]=temp;
 



        for(int m=0;m<NSTATE;m++) BdCol[m]=Bmm[m*NINPUT+k]*YswitchDiagonal[k];


        for(int m=0;m<NSTATE;m++) {
            tmp1[m]=Cmm[k*NSTATE+m];
            Cmm[k*NSTATE+m]=0;
        }
 


        for(int m=0;m<NINPUT;m++) {
            tmp2[m]=Dmm[k*NINPUT+m];
            Dmm[k*NINPUT+m]=0;
        }



        for(int m=0;m<NOUTPUT;m++) {
            for(int n=0;n<NSTATE;n++) Cmm[m*NSTATE+n]=Cmm[m*NSTATE+n]+DxCol[m]*tmp1[n];
            for(int n=0;n<NINPUT;n++) Dmm[m*NINPUT+n]=Dmm[m*NINPUT+n]+DxCol[m]*tmp2[n];
        }
 


        for(int m=0;m<NSTATE;m++) {
            for(int n=0;n<NSTATE;n++) Amm[m*NSTATE+n]=Amm[m*NSTATE+n]+BdCol[m]*Cmm[k*NSTATE+n];
            for(int n=0;n<NINPUT;n++) Bmm[m*NINPUT+n]=Bmm[m*NINPUT+n]+BdCol[m]*Dmm[k*NINPUT+n];
        }




    }

    free(tmp1);
    free(tmp2);
    free(BdCol);
    free(DxCol);
    
    free(I);
    free(IA);
    return(0);
}