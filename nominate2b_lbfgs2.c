/*
gcc -I"c:/Program Files/R/R-3.0.0/include" -L"C:/Program Files/R/R-3.0.0/bin/i386" -Wall %1.c -o %1.exe -lRlapack -lRBLAS

 nominate2b sen85kh.ord
                       
*/
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "lapacke.h"
#include "blas.h"
#include "lbfgs.h"

//#include <c:/Program Files/R/R-2.4.0/include/R_ext/Lapack.h>

typedef struct {int icount;
                double data;
} kpsorter;

void kp_clean();
void kp_reader_counts(char *);
void kp_reader(char *);
void kp_double_center(double *, double *);
void xsvd(int kpnp, int kpnq, double *, double *, double *, double *);
int real2_cmp(const void *s1, const void *s2);
void mainlbfgs(int kpnp, int kpnq, double *, double *);
double runif(void);

double runif(void){
	return(rand() / ((double)RAND_MAX + 1));
}

void timestamp ( void );

FILE *fp;
FILE *jp;
FILE *kp;
FILE *lp;
/* Remember to Dimension this 1 greater than the Number
    of Characters in the Name of the Matrix!!*/
/*char s[22]="sen90kh.ord";*/
int nq,np;
int nrowX,ncolX,NS,N;
double BETA;
static int minlegvotes=20;
static double minrcmargin=0.025;
static int ns=2;
static int jne=36;
char *kname;
int *legscalable;
int *rcscalable;
int *rollcalls;
double *X, *CONSTRAINTS;

int main(int argc, char *argv[])
{
	double time1, timedif;
	double *agree_scores, *y, *xdata;
	double *ZCOORDS3, *ZCOORDS4;
	double xsum,ysum;
	double *u, *lambda, *vt, *u_sort;
	int i,j,jj,k,kpnp,kpnq;
	int numberyeas, numbernays;
	int *L_sort;
/* clock() is part of time.h -- returns the implementation's
 * best approximationto the processor time elapsed since the
 * program was invoke, divide by CLOCKS_PER_SEC to get the time
 * in seconds           */
	time1 = (double) clock();            /* get initial time */
	time1 = time1 / CLOCKS_PER_SEC;      /*    in seconds    */

      
	printf("Welcome to NOMINATE\n");
	lp=fopen("coords.dat","w");
	jp=fopen("diagnostics.dat","w");
/*
 Read data matrix to get row and column counts
*/
	kp_reader_counts(argv[1]);
	kname = calloc(np*jne,sizeof(char));
	legscalable = calloc(np,sizeof(int));
	rcscalable = calloc(nq,sizeof(int));
	rollcalls = calloc(np*nq,sizeof(int));
/**/
/*
 Read in the data matrix
*/
	kp_reader(argv[1]);
/*
 Throw out lopsided roll calls and members who did not vote enough
*/
	kp_clean();
/*
 Dynamically allocate memory for SVD to save space
*/
	u      = calloc(nq*nq+np*np,sizeof(double));
	lambda = calloc(nq*nq+np*np,sizeof(double));
	vt     = calloc(nq*nq+np*np,sizeof(double));
	agree_scores=calloc(np*np,sizeof(double));
	y = calloc(np*np,sizeof(double));
	xdata=calloc(np*ns,sizeof(double));
/*
  Create Agreement Score Matrix and Double-Center it after converting
          it to squared distances
*/
	kp_double_center(agree_scores,y);
/*
  Decompose the Double-Centered matrix of squared distances to get
          starting values for the legislators
*/	kpnp=np;
	kpnq=np;
	xsvd(kpnp,kpnq,y,u,lambda,vt);
/*
  Compute Starting Coordinates Using Torgerson Method

 Dynamically allocate memory for Sort Routine to save space
*/
        u_sort = calloc(np,sizeof(double));
        L_sort = calloc(np,sizeof(int));

	ysum=-99999.0;
//        k_edge=0;
	for (i=0;i<np;i++)
	{
		xsum=0.0;
		fprintf(kp,"lambda[%i]=%12.7f\n",i,lambda[i]);
		for(j=0;j<jne;j=j+1)
		{
			fprintf(lp,"%c",kname[(i)*jne+j]);
		}
		for(k=0;k<ns;k++)
		{
//  NOTE THAT XDATA IS ARRANGED FORTRAN STYLE -- COLUMN STACKED
			xdata[i+(k*np)]=u[i+(k*np)]*sqrt(lambda[k]);
			xsum+=xdata[i+(k*np)]*xdata[i+(k*np)];
		}
		u_sort[i]=xsum;
		L_sort[i]=i;
		fprintf(lp,"%12.7f %12.7f\n",xdata[i],xdata[i+np]);
		if(xsum > ysum)
		{
			ysum=xsum;
//			k_edge=i;
		}
	}
//	n_real_sort=np;
/*	kp_real_sort(n_real_sort,u_sort,L_sort);
 */
	int recordcount=np;
	kpsorter recordset[np];
	for(i=0;i<np;i++)
	{
		recordset[i].icount=i;
		recordset[i].data=u_sort[i];
	}
	/*	int recordcount=sizeof(recordset)/sizeof(kpsorter);*/
	qsort(recordset,recordcount,sizeof(kpsorter),real2_cmp);
	for(i=0;i<np;i++)
	{
		fprintf(lp,"%5d %5d %10.6f\n",i,recordset[i].icount,recordset[i].data);
	}
	for(i=0;i<np;i++)
	{
		for(j=0;j<jne;j=j+1)
		{
//			fprintf(lp,"%c",kname[(i)*jne+j]);
			fprintf(lp,"%c",kname[(recordset[i].icount)*jne+j]);
		}
		fprintf(lp,"%12.7f %12.7f %12.7f %d\n",xdata[recordset[i].icount],xdata[recordset[i].icount+np],recordset[i].data,recordset[i].icount);
	}
	for(i=0;i<np;i++)
	{
		for(k=0;k<ns;k++)
		{
			xdata[i+(k*np)]=xdata[i+(k*np)]*(1.0/sqrt(ysum));
		}
		for(j=0;j<jne;j=j+1)
		{
			fprintf(lp,"%c",kname[(i)*jne+j]);
		}
		fprintf(lp,"%12.7f %12.7f\n",xdata[i],xdata[i+np]);
	}
	nrowX=np;
	ncolX=nq;
	NS=ns;
	N=(NS*(nrowX+2*ncolX))-((NS*(NS+1))/2);
	BETA=1.2;
// SET UP FOR LBFGS CALL -- X is the Roll Call Matrix	
	X = (double *) malloc (nrowX*ncolX*sizeof(double));
	CONSTRAINTS = (double *) malloc (((NS*(nrowX+2*ncolX))+1)*sizeof(double));
	ZCOORDS3 = (double *) malloc (((NS*(nrowX+2*ncolX))+1)*sizeof(double));
	ZCOORDS4 = (double *) malloc (((NS*(nrowX+2*ncolX))+1)*sizeof(double));
	for(i=0;i<(NS*(2*ncolX));i++)
	{
		CONSTRAINTS[i]=1.0;
	}
	for(i=0;i<(NS*nrowX);i++)
	{
		CONSTRAINTS[i+(NS*(2*ncolX))]=0.0;
		if(i<((NS*nrowX)-(NS*(NS+1))/2))
		{
			CONSTRAINTS[i+(NS*(2*ncolX))]=1.0;
		}
	}
// Pass Roll Call Matrix into X
	numberyeas=0;
	numbernays=0;
	for(i=0;i<nrowX;i++)
	{
		for(j=0;j<ncolX;j++)
		{
			X[i*ncolX+j]=0.0;
			if(rollcalls[i*ncolX+j]==1){
				X[i*ncolX+j]=1.0;
				numberyeas=numberyeas+1;
			}
			if(rollcalls[i*ncolX+j]==6){
				X[i*ncolX+j]=6.0;
				numbernays=numbernays+1;
			}
		}
	}
	for(i=0;i<nrowX;i++)
	{
		for(j=0;j<jne;j=j+1)
		{
			fprintf(lp,"%c",kname[(i)*jne+j]);
		}
		for(j=0;j<ncolX;j++)
		{
			fprintf(lp,"%2.0f",X[i*ncolX+j]);
		}
		fprintf(lp,"\n");
	}
//	
	printf("Parameters to use in LBFGS %8d %8d %8d %8d %12.7f %10d %10d\n",nrowX,ncolX,NS,N,BETA,numberyeas,numbernays);
	fprintf(kp,"Parameters to use in LBFGS %8d %8d %8d %8d %12.7f %10d %10d\n",nrowX,ncolX,NS,N,BETA,numberyeas,numbernays);
//
//  CALL Limited-Memory Broyden-Fletcher-Goldfarb-Shanno
//     ZCOORDS3 are the starting coordinates
//     Solution Passed back in ZCOORDS4
//
	for(j=0;j<ncolX;j++)
	{
//  Each Roll Call has 2*s parameters -- Yea and Nay outcomes
		for(jj=0;jj<(2*NS);jj++)
		{
//			ZCOORDS3[j*NS+jj]=ZCOORDS2[j*NS+jj];
			ZCOORDS3[j*(2*NS)+jj]=(0.5-runif())*2.0;
		}
	}
	for(j=0;j<nrowX;j++)
	{
//  LEGISLATORS MUST BE STACKED BY ROW -- C STYLE
		for(jj=0;jj<NS;jj++)
		{
			ZCOORDS3[j*NS+jj+ncolX*NS*2]=xdata[j + jj*nrowX];
		}
	}
//	ZCOORDS3[((nrowX+2*ncolX)*NS)]=BETA;
//  
	mainlbfgs(nrowX,ncolX,ZCOORDS4,ZCOORDS3);
//
	timedif = ( ((double) clock()) / CLOCKS_PER_SEC) - time1;
	printf("\nThe total elapsed time of the program is %12.3f seconds\n", timedif);
	fprintf(kp,"\nThe total elapsed time of the program is %12.3f seconds\n", timedif);
	timestamp ( );
//
	free(u);
	free(lambda);
	free(vt);
	free(agree_scores);
	free(y);
	free(xdata);
	free(kname);
	free(legscalable);
	free(rcscalable);
	free(rollcalls);
	free(X);
	free(CONSTRAINTS);
	free(ZCOORDS3);
	free(ZCOORDS4);
/**/
	printf("This is the End of the Program Beavis\n");
	fclose(kp);
	fclose(lp);
	return(0);
}
/*


Subroutine kp_reader_counts -- Reads Roll Call Matrix and does an intial
                       counts of rows and columns
*/
void kp_reader_counts(char *argv)
{
	int i,j,imax,jmax,kk;
	char str0[1830];
	char *dummyread;
	jmax = 0;
	imax = 0;
	kp=fopen("hello2.dat","w");
	if((fp=fopen(argv,"r"))==NULL)
		printf("The file was not opened\n");
	else
	{
		printf("The file was opened successfully Beavis\n");
/* Get the record length of the Roll Call File */
		fgets(str0,1820,fp);
		j=strlen(str0);
		jmax=j;
		j=j-(jne+1);
		printf("There are %d Roll Calls\n",j);
		fprintf(kp,"There are %d Roll Calls\n",j);
		nq=j;
		fclose(fp);
	}
	dummyread = calloc(jmax,sizeof(char));
/* Reopen the Roll Call File and read in the Names and
     Find the number of Records */
	if((fp=fopen(argv,"r"))==NULL)
		printf("The file was not opened\n");
	else
	{
		printf("The file was opened successfully for the 2nd Time Beavis\n");
		for(i=0;;i=i+1)
		{
			for(j=0;j<jmax;j=j+1)
			{
				fscanf(fp,"%c",&dummyread[j]);
			}
			kk=feof(fp);
			if(kk == 0) imax=i;
			if(kk != 0) break;
		}
		printf("Number of Rows and Columns of Roll Call Matrix %d %d\n",imax+1,jmax);
		fprintf(kp,"Number of Rows and Columns of Roll Call Matrix %d %d\n",imax+1,jmax);
		np=imax+1;
		fprintf(kp,"Number of Legislators and Number of Roll Calls %d %d\n",np,nq);
		printf("Number of Legislators and Number of Roll Calls %d %d\n",np,nq);
	}
	fclose(fp);
	free(dummyread);
}
/*


Subroutine kp_reader -- Reads Roll Call Matrix and does an intial
                       recoding of the roll calls to 1=Yea and 6=Nay

*/
	void kp_reader(char *argv)
	{
	int i,j,imax,jmax,kk;
	char *ldata;
	jmax = nq+jne+1;
	imax = 0;
	ldata = calloc(nq,sizeof(char));

/* Reopen the Roll Call File and read in the Names and
     Find the number of Records */
	if((fp=fopen(argv,"r"))==NULL)
/*	if((fp=fopen(s,"r"))==NULL) */
		printf("The file was not opened\n");
	else
	{
		printf("The file was opened successfully for the 2nd Time Beavis\n");
		for(i=0;;i=i+1)
/*			for(i=0;i<700;i=i+1)  */
		{
			for(j=0;j<jne;j=j+1)
			{
				fscanf(fp,"%c",&kname[(i)*jne+j]);
			}
			for(j=jne;j<jmax;j=j+1)
			{
				fscanf(fp,"%c",&ldata[j-jne]);
			} 
			kk=feof(fp);
			if(kk != 0) fprintf(kp,"end of file marker %d %d\n",i,kk);
			if(kk==0){
			for(j=0;j<jne;j=j+1)
			{
				fprintf(kp,"%c",kname[(i)*jne+j]);
			}
			fprintf(kp,"\n");
			for(j=0;j<nq;j=j+1)
			{
/*				if(ldata[j-jne] == '0') rollcalls[i][j-jne]=0;
*/
				if(ldata[j] == '0') rollcalls[i*nq+j]=0;  
				if(ldata[j] == '1') rollcalls[i*nq+j]=1;
				if(ldata[j] == '2') rollcalls[i*nq+j]=1;
				if(ldata[j] == '3') rollcalls[i*nq+j]=1;
				if(ldata[j] == '4') rollcalls[i*nq+j]=6;
				if(ldata[j] == '5') rollcalls[i*nq+j]=6;
				if(ldata[j] == '6') rollcalls[i*nq+j]=6;
				if(ldata[j] == '7') rollcalls[i*nq+j]=0;
				if(ldata[j] == '8') rollcalls[i*nq+j]=0;
				if(ldata[j] == '9') rollcalls[i*nq+j]=0;
			}
			}
			if(kk == 0) imax=i;
			if(kk != 0) break;
		}
		printf("Number of Rows and Columns of Roll Call Matrix %d %d\n",imax+1,jmax-1);
		fprintf(kp,"Number of Rows and Columns of Roll Call Matrix %d %d\n",imax+1,jmax-1);
		np=imax+1;
		fprintf(kp,"Number of Legislators and Number of Roll Calls %d %d\n",np,nq);
		printf("Number of Legislators and Number of Roll Calls %d %d\n",np,nq);
	}
	fclose(fp);
	free(ldata);
}
/*
 
Subroutine kp_clean -- Drops lopsided roll calls and Legislators who
                       do not vote enough times
*/
void kp_clean()
{
	int *rollcalls2;
	double minkyeskno,xmargin;
	int i,j,kk,kvote,kyes,kno,npscale,nqscale;
	int ktrack, kstack;
	char *kkname;
	rollcalls2=calloc(np*nq,sizeof(int));
	kkname=calloc(np*jne,sizeof(char));
	npscale=0;
    for(i=0;i<np;i++)
	{
		kvote=0;
		legscalable[i]=0;
        for(j=0;j<nq;j++) {
		if(rollcalls[i*nq+j]!=0)kvote++;
	}
	if(kvote<minlegvotes)legscalable[i]=1;
	if(kvote>minlegvotes)npscale++;
		for(j=0;j<jne;j=j+1)
		{
			fprintf(kp,"%c",kname[(i)*jne+j]);
		}
		fprintf(kp,"%d,%d\n",kvote,legscalable[i]);
	}
    nqscale=0;
    for(j=0;j<nq;j++){
	    kyes=0;
	    kno=0;
	    for(i=0;i<np;i++){
		    if(rollcalls[i*nq+j]==1)kyes++;
		    if(rollcalls[i*nq+j]==6)kno++;
	    }
	    xmargin=0.0;
	    rcscalable[j]=0;
	    if(kyes<=kno)minkyeskno=kyes;
	    if(kyes>kno)minkyeskno=kno;
	    if((kyes+kno)>0)xmargin=minkyeskno/(kyes+kno);
	    if(xmargin<minrcmargin)rcscalable[j]=1;
	    if(xmargin>=minrcmargin)nqscale++;
	    fprintf(kp,"%5i %5i %5i %10f %5i\n",j+1,kyes,kno,xmargin,rcscalable[j]);
    }
    fprintf(kp,"Number of Scalable Legislators and Roll Calls %5d %5d\n",npscale,nqscale);
    printf("Number of Scalable Legislators and Roll Calls %5d %5d\n",npscale,nqscale);
/*
 Delete the unscalable legislators from the roll call matrix
 */
    kk=0;
    ktrack=0;
    for(i=0;i<np;i++)
    {
	    if(legscalable[i]==0)
	    {
   	       for(j=0;j<nq;j++) {
		    rollcalls2[j+(kk*nq)]=rollcalls[i*nq+j];
	       }
	       for(j=0;j<jne;j++) {
		       kkname[j+(kk*jne)]=kname[(i)*jne+j];
	       }
	       kk=kk+1;
	    }
	    if(legscalable[i]==1)
	    {
		    if(ktrack==0)
		    {
			    fprintf(kp,"Unscalable Legislators\n");
		    }
		    for(j=0;j<jne;j=j+1)
		    {
			    fprintf(kp,"%c",kname[(i)*jne+j]);
		    }
		    fprintf(kp,"\n");
		    ktrack=1;
	    }
    }
//    nqscale=nq;
    for(i=0;i<kk;i++)
    {
	    kstack=0;
		    for(j=0;j<nq;j++) {
//  STORE ROLLCALLS BY ROWS STACKED IN C STYLE
			    if(rcscalable[j]==0)
			    {
				    rollcalls[kstack+nqscale*i]=rollcalls2[j+(i*nq)];
				    kstack=kstack+1;
			    }
		    }
		    for(j=0;j<jne;j++) {
			    kname[(i)*jne+j]=kkname[j+(i*jne)];
		    }
    }
    fprintf(kp,"Number of Scalable Legislators and Roll Calls %5d %5d\n",kk,nqscale);
    printf("Number of Scalable Legislators and Roll Calls %5d %5d\n",kk,nqscale);
/*
Reset Number of Legislators
*/
    np=npscale;
    nq=nqscale;
    free(rollcalls2);
    free(kkname);
}
/*

Subroutine kp_double_center -- Computes Agreement Score Matrix and
                                Double-Centers it
*/
void kp_double_center(double *agree_scores, double *y)
{
	int i,ii,j,kboth,kbothyea,kbothnay;
	double xboth,xagree,matrix_mean,dist_agree;
	double *row_means, *column_means;
	row_means=calloc(np,sizeof(double));
	column_means=calloc(np,sizeof(double));
/*
Compute Agreement Score Matrix
*/
	for(i=0;i<np;i++)
	{
		row_means[i]=0.0;
		column_means[i]=0.0;
	}
	agree_scores[(np*np)-1]=1.0;
	for(i=0;i<np-1;i++)
	{
		agree_scores[i*(np+1)]=1.0;
		for(ii=i+1;ii<np;ii++)
		{
			kboth=0;
			kbothyea=0;
			kbothnay=0;
			for(j=0;j<nq;j++)
			{
				if(rcscalable[j]==0)
				{
					if((rollcalls[i*nq+j]!=0) && (rollcalls[ii*nq+j]!=0))
					{
						kboth++;
					}
					if((rollcalls[i*nq+j]==1) && (rollcalls[ii*nq+j]==1))
					{
						kbothyea++;
					}
					if((rollcalls[i*nq+j]==6) && (rollcalls[ii*nq+j]==6))
					{
						kbothnay++;
					}
				}
			}
			if(kboth > 0)
			{
				xboth=kboth;
				xagree=kbothyea+kbothnay;
				agree_scores[ii+i*np]=xagree/xboth;
				agree_scores[i+ii*np]=xagree/xboth;
			}
			else
			{
				agree_scores[ii+i*np]=0.5;
				agree_scores[i+ii*np]=0.5;
			}
		}
	}
	for(ii=0;ii<np;ii++)
	{
        	for(j=0;j<jne;j=j+1)
	        {
		        fprintf(kp,"%c",kname[(ii)*jne+j]);
        	}
	        for(i=0;i<np;i++)
        	{
	        	fprintf(kp,"%8.4f",agree_scores[i+ii*np]);
        	}
	        fprintf(kp,"\n");
	}
/*
  Double-Center the Agreement Score Matrix
*/
	matrix_mean=0.0;
	for(ii=0;ii<np;ii++)
	{
		for(i=0;i<np;i++)
		{
			dist_agree=(1.0-agree_scores[i+ii*np])*(1.0-agree_scores[i+ii*np]);
			row_means[ii]+=dist_agree;
			matrix_mean+=dist_agree;
		}
		row_means[ii]=row_means[ii]/np;
	}
	matrix_mean=matrix_mean/(np*np);
	for(ii=0;ii<np;ii++)
	{
		for(i=0;i<np;i++)
		{
			dist_agree=(1.0-agree_scores[i+ii*np])*(1.0-agree_scores[i+ii*np]);
			y[i+ii*np]=(dist_agree-row_means[ii]-row_means[i]+matrix_mean)/(-2.0);
		}
	}
	for(ii=0;ii<np;ii++)
	{
		for(j=0;j<jne;j=j+1)
		{
			fprintf(kp,"%c",kname[(ii)*jne+j]);
		}
		for(i=0;i<np;i++)
		{
			fprintf(kp,"%10.4f",y[i+ii*np]);
		}
		fprintf(kp,"\n");
	}
	free(row_means);
	free(column_means);
}
/*

Singular Value Decomposition Subroutine

*/
void xsvd(int kpnp, int kpnq, double *y, double *u, double *lambda, double *vt) {
/*   double a[25000], work[25000];
*/
	
   double *a, *work;
   double sumulv, svd_error_sum, svd_error_sum_2;
   int i, j, jj;
   int  info = 12;
   int  lwork= kpnp*kpnp+kpnq*kpnq;
   int  lda,ldu,ldvt;

   a     = calloc( kpnp*kpnq, sizeof(double));
   work  = calloc( lwork, sizeof(double));

   fprintf(kp,"entering svd...\n");

   lda  = kpnp;
   ldu  = kpnp;
   ldvt = kpnq;
   for (i=0;i<lwork;i++){
     work[i] = 0;
   }
   
   fprintf(kp,"lwork=%i\n",lwork);

   for (j=0;j<kpnq;j++) {
      for (i=0;i<kpnp;i++) {
	      a[(j*kpnp)+i] = y[(j*kpnp)+i];
      }
   }

   dgesvd_("A","A", &kpnp, &kpnq, a, &lda, lambda,
	    u, &ldu, vt, &ldvt, work, &lwork, &info);

   fprintf(kp,"Info = %i\n",info);
/*
  Do simple check of SVD

*/
   svd_error_sum=0.0;
   svd_error_sum_2=0.0;
   for (i=0;i<kpnp;i++)
   {
	   for (jj=0;jj<kpnq;jj++)
	   {
		   sumulv=0.0;
		   for (j=0;j<kpnq;j++)
		   {
			   sumulv+=u[(j*kpnp)+i]*lambda[j]*vt[j+(jj*kpnq)];
		   }
		   svd_error_sum+=(y[i+(jj*kpnp)]-sumulv)*(y[i+(jj*kpnp)]-sumulv);
		   svd_error_sum_2+=fabs(y[i+(jj*kpnp)]-sumulv);
/*		   fprintf(kp,"%12.7g %12.7g\n",y[i+(jj*kpnp)],sumulv);
 */
		   }
   }
   fprintf(kp,"%12.7g %12.7g\n",svd_error_sum,svd_error_sum_2);
   fprintf(kp,"Leaving svd...\n");	 
   free(work);
   free(a);
}

int real2_cmp(const void *s1, const void *s2)
{
	kpsorter *p1=(kpsorter *)s1;
	kpsorter *p2=(kpsorter *)s2;

	if(p1->data < p2->data)
		return -1;
	else if (p1->data == p2->data)
		return 0;
	else
		return 1;
}



/*
 *
Listing of LBFGS.H file

*/
/*
 *      C library of Limited memory BFGS (L-BFGS).
 *
 * Copyright (c) 1990, Jorge Nocedal
 * Copyright (c) 2007-2010 Naoaki Okazaki
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* $Id$ */

#ifndef __LBFGS_H__
#define __LBFGS_H__

#ifdef  __cplusplus
extern "C" {
#endif/*__cplusplus*/

/*
 * The default precision of floating point values is 64bit (double).
 */
#ifndef LBFGS_FLOAT
#define LBFGS_FLOAT     64
#endif/*LBFGS_FLOAT*/

/*
 * Activate optimization routines for IEEE754 floating point values.
 */
#ifndef LBFGS_IEEE_FLOAT
#define LBFGS_IEEE_FLOAT    1
#endif/*LBFGS_IEEE_FLOAT*/

#if     LBFGS_FLOAT == 32
typedef float lbfgsfloatval_t;

#elif   LBFGS_FLOAT == 64
typedef double lbfgsfloatval_t;

#else
#error "libLBFGS supports single (float; LBFGS_FLOAT = 32) or double (double; LBFGS_FLOAT=64) precision only."

#endif


/** 
 * \addtogroup liblbfgs_api libLBFGS API
 * @{
 *
 *  The libLBFGS API.
 */

/**
 * Return values of lbfgs().
 * 
 *  Roughly speaking, a negative value indicates an error.
 */
enum {
    /** L-BFGS reaches convergence. */
    LBFGS_SUCCESS = 0,
    LBFGS_CONVERGENCE = 0,
    LBFGS_STOP,
    /** The initial variables already minimize the objective function. */
    LBFGS_ALREADY_MINIMIZED,

    /** Unknown error. */
    LBFGSERR_UNKNOWNERROR = -1024,
    /** Logic error. */
    LBFGSERR_LOGICERROR,
    /** Insufficient memory. */
    LBFGSERR_OUTOFMEMORY,
    /** The minimization process has been canceled. */
    LBFGSERR_CANCELED,
    /** Invalid number of variables specified. */
    LBFGSERR_INVALID_N,
    /** Invalid number of variables (for SSE) specified. */
    LBFGSERR_INVALID_N_SSE,
    /** The array x must be aligned to 16 (for SSE). */
    LBFGSERR_INVALID_X_SSE,
    /** Invalid parameter lbfgs_parameter_t::epsilon specified. */
    LBFGSERR_INVALID_EPSILON,
    /** Invalid parameter lbfgs_parameter_t::past specified. */
    LBFGSERR_INVALID_TESTPERIOD,
    /** Invalid parameter lbfgs_parameter_t::delta specified. */
    LBFGSERR_INVALID_DELTA,
    /** Invalid parameter lbfgs_parameter_t::linesearch specified. */
    LBFGSERR_INVALID_LINESEARCH,
    /** Invalid parameter lbfgs_parameter_t::max_step specified. */
    LBFGSERR_INVALID_MINSTEP,
    /** Invalid parameter lbfgs_parameter_t::max_step specified. */
    LBFGSERR_INVALID_MAXSTEP,
    /** Invalid parameter lbfgs_parameter_t::ftol specified. */
    LBFGSERR_INVALID_FTOL,
    /** Invalid parameter lbfgs_parameter_t::wolfe specified. */
    LBFGSERR_INVALID_WOLFE,
    /** Invalid parameter lbfgs_parameter_t::gtol specified. */
    LBFGSERR_INVALID_GTOL,
    /** Invalid parameter lbfgs_parameter_t::xtol specified. */
    LBFGSERR_INVALID_XTOL,
    /** Invalid parameter lbfgs_parameter_t::max_linesearch specified. */
    LBFGSERR_INVALID_MAXLINESEARCH,
    /** Invalid parameter lbfgs_parameter_t::orthantwise_c specified. */
    LBFGSERR_INVALID_ORTHANTWISE,
    /** Invalid parameter lbfgs_parameter_t::orthantwise_start specified. */
    LBFGSERR_INVALID_ORTHANTWISE_START,
    /** Invalid parameter lbfgs_parameter_t::orthantwise_end specified. */
    LBFGSERR_INVALID_ORTHANTWISE_END,
    /** The line-search step went out of the interval of uncertainty. */
    LBFGSERR_OUTOFINTERVAL,
    /** A logic error occurred; alternatively, the interval of uncertainty
        became too small. */
    LBFGSERR_INCORRECT_TMINMAX,
    /** A rounding error occurred; alternatively, no line-search step
        satisfies the sufficient decrease and curvature conditions. */
    LBFGSERR_ROUNDING_ERROR,
    /** The line-search step became smaller than lbfgs_parameter_t::min_step. */
    LBFGSERR_MINIMUMSTEP,
    /** The line-search step became larger than lbfgs_parameter_t::max_step. */
    LBFGSERR_MAXIMUMSTEP,
    /** The line-search routine reaches the maximum number of evaluations. */
    LBFGSERR_MAXIMUMLINESEARCH,
    /** The algorithm routine reaches the maximum number of iterations. */
    LBFGSERR_MAXIMUMITERATION,
    /** Relative width of the interval of uncertainty is at most
        lbfgs_parameter_t::xtol. */
    LBFGSERR_WIDTHTOOSMALL,
    /** A logic error (negative line-search step) occurred. */
    LBFGSERR_INVALIDPARAMETERS,
    /** The current search direction increases the objective function value. */
    LBFGSERR_INCREASEGRADIENT,
};

/**
 * Line search algorithms.
 */
enum {
    /** The default algorithm (MoreThuente method). */
    LBFGS_LINESEARCH_DEFAULT = 0,
    /** MoreThuente method proposd by More and Thuente. */
    LBFGS_LINESEARCH_MORETHUENTE = 0,
    /**
     * Backtracking method with the Armijo condition.
     *  The backtracking method finds the step length such that it satisfies
     *  the sufficient decrease (Armijo) condition,
     *    - f(x + a * d) <= f(x) + lbfgs_parameter_t::ftol * a * g(x)^T d,
     *
     *  where x is the current point, d is the current search direction, and
     *  a is the step length.
     */
    LBFGS_LINESEARCH_BACKTRACKING_ARMIJO = 1,
    /** The backtracking method with the defualt (regular Wolfe) condition. */
    LBFGS_LINESEARCH_BACKTRACKING = 2,
    /**
     * Backtracking method with regular Wolfe condition.
     *  The backtracking method finds the step length such that it satisfies
     *  both the Armijo condition (LBFGS_LINESEARCH_BACKTRACKING_ARMIJO)
     *  and the curvature condition,
     *    - g(x + a * d)^T d >= lbfgs_parameter_t::wolfe * g(x)^T d,
     *
     *  where x is the current point, d is the current search direction, and
     *  a is the step length.
     */
    LBFGS_LINESEARCH_BACKTRACKING_WOLFE = 2,
    /**
     * Backtracking method with strong Wolfe condition.
     *  The backtracking method finds the step length such that it satisfies
     *  both the Armijo condition (LBFGS_LINESEARCH_BACKTRACKING_ARMIJO)
     *  and the following condition,
     *    - |g(x + a * d)^T d| <= lbfgs_parameter_t::wolfe * |g(x)^T d|,
     *
     *  where x is the current point, d is the current search direction, and
     *  a is the step length.
     */
    LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE = 3,
};

/**
 * L-BFGS optimization parameters.
 *  Call lbfgs_parameter_init() function to initialize parameters to the
 *  default values.
 */
typedef struct {
    /**
     * The number of corrections to approximate the inverse hessian matrix.
     *  The L-BFGS routine stores the computation results of previous \ref m
     *  iterations to approximate the inverse hessian matrix of the current
     *  iteration. This parameter controls the size of the limited memories
     *  (corrections). The default value is \c 6. Values less than \c 3 are
     *  not recommended. Large values will result in excessive computing time.
     */
    int             m;

    /**
     * Epsilon for convergence test.
     *  This parameter determines the accuracy with which the solution is to
     *  be found. A minimization terminates when
     *      ||g|| < \ref epsilon * max(1, ||x||),
     *  where ||.|| denotes the Euclidean (L2) norm. The default value is
     *  \c 1e-5.
     */
    lbfgsfloatval_t epsilon;

    /**
     * Distance for delta-based convergence test.
     *  This parameter determines the distance, in iterations, to compute
     *  the rate of decrease of the objective function. If the value of this
     *  parameter is zero, the library does not perform the delta-based
     *  convergence test. The default value is \c 0.
     */
    int             past;

    /**
     * Delta for convergence test.
     *  This parameter determines the minimum rate of decrease of the
     *  objective function. The library stops iterations when the
     *  following condition is met:
     *      (f' - f) / f < \ref delta,
     *  where f' is the objective value of \ref past iterations ago, and f is
     *  the objective value of the current iteration.
     *  The default value is \c 0.
     */
    lbfgsfloatval_t delta;

    /**
     * The maximum number of iterations.
     *  The lbfgs() function terminates an optimization process with
     *  ::LBFGSERR_MAXIMUMITERATION status code when the iteration count
     *  exceedes this parameter. Setting this parameter to zero continues an
     *  optimization process until a convergence or error. The default value
     *  is \c 0.
     */
    int             max_iterations;

    /**
     * The line search algorithm.
     *  This parameter specifies a line search algorithm to be used by the
     *  L-BFGS routine.
     */
    int             linesearch;

    /**
     * The maximum number of trials for the line search.
     *  This parameter controls the number of function and gradients evaluations
     *  per iteration for the line search routine. The default value is \c 20.
     */
    int             max_linesearch;

    /**
     * The minimum step of the line search routine.
     *  The default value is \c 1e-20. This value need not be modified unless
     *  the exponents are too large for the machine being used, or unless the
     *  problem is extremely badly scaled (in which case the exponents should
     *  be increased).
     */
    lbfgsfloatval_t min_step;

    /**
     * The maximum step of the line search.
     *  The default value is \c 1e+20. This value need not be modified unless
     *  the exponents are too large for the machine being used, or unless the
     *  problem is extremely badly scaled (in which case the exponents should
     *  be increased).
     */
    lbfgsfloatval_t max_step;

    /**
     * A parameter to control the accuracy of the line search routine.
     *  The default value is \c 1e-4. This parameter should be greater
     *  than zero and smaller than \c 0.5.
     */
    lbfgsfloatval_t ftol;

    /**
     * A coefficient for the Wolfe condition.
     *  This parameter is valid only when the backtracking line-search
     *  algorithm is used with the Wolfe condition,
     *  ::LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE or
     *  ::LBFGS_LINESEARCH_BACKTRACKING_WOLFE .
     *  The default value is \c 0.9. This parameter should be greater
     *  the \ref ftol parameter and smaller than \c 1.0.
     */
    lbfgsfloatval_t wolfe;

    /**
     * A parameter to control the accuracy of the line search routine.
     *  The default value is \c 0.9. If the function and gradient
     *  evaluations are inexpensive with respect to the cost of the
     *  iteration (which is sometimes the case when solving very large
     *  problems) it may be advantageous to set this parameter to a small
     *  value. A typical small value is \c 0.1. This parameter shuold be
     *  greater than the \ref ftol parameter (\c 1e-4) and smaller than
     *  \c 1.0.
     */
    lbfgsfloatval_t gtol;

    /**
     * The machine precision for floating-point values.
     *  This parameter must be a positive value set by a client program to
     *  estimate the machine precision. The line search routine will terminate
     *  with the status code (::LBFGSERR_ROUNDING_ERROR) if the relative width
     *  of the interval of uncertainty is less than this parameter.
     */
    lbfgsfloatval_t xtol;

    /**
     * Coeefficient for the L1 norm of variables.
     *  This parameter should be set to zero for standard minimization
     *  problems. Setting this parameter to a positive value activates
     *  Orthant-Wise Limited-memory Quasi-Newton (OWL-QN) method, which
     *  minimizes the objective function F(x) combined with the L1 norm |x|
     *  of the variables, {F(x) + C |x|}. This parameter is the coeefficient
     *  for the |x|, i.e., C. As the L1 norm |x| is not differentiable at
     *  zero, the library modifies function and gradient evaluations from
     *  a client program suitably; a client program thus have only to return
     *  the function value F(x) and gradients G(x) as usual. The default value
     *  is zero.
     */
    lbfgsfloatval_t orthantwise_c;

    /**
     * Start index for computing L1 norm of the variables.
     *  This parameter is valid only for OWL-QN method
     *  (i.e., \ref orthantwise_c != 0). This parameter b (0 <= b < N)
     *  specifies the index number from which the library computes the
     *  L1 norm of the variables x,
     *      |x| := |x_{b}| + |x_{b+1}| + ... + |x_{N}| .
     *  In other words, variables x_1, ..., x_{b-1} are not used for
     *  computing the L1 norm. Setting b (0 < b < N), one can protect
     *  variables, x_1, ..., x_{b-1} (e.g., a bias term of logistic
     *  regression) from being regularized. The default value is zero.
     */
    int             orthantwise_start;

    /**
     * End index for computing L1 norm of the variables.
     *  This parameter is valid only for OWL-QN method
     *  (i.e., \ref orthantwise_c != 0). This parameter e (0 < e <= N)
     *  specifies the index number at which the library stops computing the
     *  L1 norm of the variables x,
     */
    int             orthantwise_end;
} lbfgs_parameter_t;


/**
 * Callback interface to provide objective function and gradient evaluations.
 *
 *  The lbfgs() function call this function to obtain the values of objective
 *  function and its gradients when needed. A client program must implement
 *  this function to evaluate the values of the objective function and its
 *  gradients, given current values of variables.
 *  
 *  @param  instance    The user data sent for lbfgs() function by the client.
 *  @param  x           The current values of variables.
 *  @param  g           The gradient vector. The callback function must compute
 *                      the gradient values for the current variables.
 *  @param  n           The number of variables.
 *  @param  step        The current step of the line search routine.
 *  @retval lbfgsfloatval_t The value of the objective function for the current
 *                          variables.
 */
typedef lbfgsfloatval_t (*lbfgs_evaluate_t)(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    );

/**
 * Callback interface to receive the progress of the optimization process.
 *
 *  The lbfgs() function call this function for each iteration. Implementing
 *  this function, a client program can store or display the current progress
 *  of the optimization process.
 *
 *  @param  instance    The user data sent for lbfgs() function by the client.
 *  @param  x           The current values of variables.
 *  @param  g           The current gradient values of variables.
 *  @param  fx          The current value of the objective function.
 *  @param  xnorm       The Euclidean norm of the variables.
 *  @param  gnorm       The Euclidean norm of the gradients.
 *  @param  step        The line-search step used for this iteration.
 *  @param  n           The number of variables.
 *  @param  k           The iteration count.
 *  @param  ls          The number of evaluations called for this iteration.
 *  @retval int         Zero to continue the optimization process. Returning a
 *                      non-zero value will cancel the optimization process.
 */
typedef int (*lbfgs_progress_t)(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    );

/*
A user must implement a function compatible with ::lbfgs_evaluate_t (evaluation
callback) and pass the pointer to the callback function to lbfgs() arguments.
Similarly, a user can implement a function compatible with ::lbfgs_progress_t
(progress callback) to obtain the current progress (e.g., variables, function
value, ||G||, etc) and to cancel the iteration process if necessary.
Implementation of a progress callback is optional: a user can pass \c NULL if
progress notification is not necessary.

In addition, a user must preserve two requirements:
    - The number of variables must be multiples of 16 (this is not 4).
    - The memory block of variable array ::x must be aligned to 16.

This algorithm terminates an optimization
when:

    ||G|| < \epsilon \cdot \max(1, ||x||) .

In this formula, ||.|| denotes the Euclidean norm.
*/

/**
 * Start a L-BFGS optimization.
 *
 *  @param  n           The number of variables.
 *  @param  x           The array of variables. A client program can set
 *                      default values for the optimization and receive the
 *                      optimization result through this array. This array
 *                      must be allocated by ::lbfgs_malloc function
 *                      for libLBFGS built with SSE/SSE2 optimization routine
 *                      enabled. The library built without SSE/SSE2
 *                      optimization does not have such a requirement.
 *  @param  ptr_fx      The pointer to the variable that receives the final
 *                      value of the objective function for the variables.
 *                      This argument can be set to \c NULL if the final
 *                      value of the objective function is unnecessary.
 *  @param  proc_evaluate   The callback function to provide function and
 *                          gradient evaluations given a current values of
 *                          variables. A client program must implement a
 *                          callback function compatible with \ref
 *                          lbfgs_evaluate_t and pass the pointer to the
 *                          callback function.
 *  @param  proc_progress   The callback function to receive the progress
 *                          (the number of iterations, the current value of
 *                          the objective function) of the minimization
 *                          process. This argument can be set to \c NULL if
 *                          a progress report is unnecessary.
 *  @param  instance    A user data for the client program. The callback
 *                      functions will receive the value of this argument.
 *  @param  param       The pointer to a structure representing parameters for
 *                      L-BFGS optimization. A client program can set this
 *                      parameter to \c NULL to use the default parameters.
 *                      Call lbfgs_parameter_init() function to fill a
 *                      structure with the default values.
 *  @retval int         The status code. This function returns zero if the
 *                      minimization process terminates without an error. A
 *                      non-zero value indicates an error.
 */
int lbfgs(
    int n,
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *ptr_fx,
    lbfgs_evaluate_t proc_evaluate,
    lbfgs_progress_t proc_progress,
    void *instance,
    lbfgs_parameter_t *param
    );

/**
 * Initialize L-BFGS parameters to the default values.
 *
 *  Call this function to fill a parameter structure with the default values
 *  and overwrite parameter values if necessary.
 *
 *  @param  param       The pointer to the parameter structure.
 */
void lbfgs_parameter_init(lbfgs_parameter_t *param);

/**
 * Allocate an array for variables.
 *
 *  This function allocates an array of variables for the convenience of
 *  ::lbfgs function; the function has a requreiemt for a variable array
 *  when libLBFGS is built with SSE/SSE2 optimization routines. A user does
 *  not have to use this function for libLBFGS built without SSE/SSE2
 *  optimization.
 *  
 *  @param  n           The number of variables.
 */
lbfgsfloatval_t* lbfgs_malloc(int n);

/**
 * Free an array of variables.
 *  
 *  @param  x           The array of variables allocated by ::lbfgs_malloc
 *                      function.
 */
void lbfgs_free(lbfgsfloatval_t *x);

/** @} */

#ifdef  __cplusplus
}
#endif/*__cplusplus*/



/**
@mainpage libLBFGS: a library of Limited-memory Broyden-Fletcher-Goldfarb-Shanno (L-BFGS)

@section intro Introduction

This library is a C port of the implementation of Limited-memory
Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) method written by Jorge Nocedal.
The original FORTRAN source code is available at:
http://www.ece.northwestern.edu/~nocedal/lbfgs.html

The L-BFGS method solves the unconstrainted minimization problem,

<pre>
    minimize F(x), x = (x1, x2, ..., xN),
</pre>

only if the objective function F(x) and its gradient G(x) are computable. The
well-known Newton's method requires computation of the inverse of the hessian
matrix of the objective function. However, the computational cost for the
inverse hessian matrix is expensive especially when the objective function
takes a large number of variables. The L-BFGS method iteratively finds a
minimizer by approximating the inverse hessian matrix by information from last
m iterations. This innovation saves the memory storage and computational time
drastically for large-scaled problems.

Among the various ports of L-BFGS, this library provides several features:
- <b>Optimization with L1-norm (Orthant-Wise Limited-memory Quasi-Newton
  (OWL-QN) method)</b>:
  In addition to standard minimization problems, the library can minimize
  a function F(x) combined with L1-norm |x| of the variables,
  {F(x) + C |x|}, where C is a constant scalar parameter. This feature is
  useful for estimating parameters of sparse log-linear models (e.g.,
  logistic regression and maximum entropy) with L1-regularization (or
  Laplacian prior).
- <b>Clean C code</b>:
  Unlike C codes generated automatically by f2c (Fortran 77 into C converter),
  this port includes changes based on my interpretations, improvements,
  optimizations, and clean-ups so that the ported code would be well-suited
  for a C code. In addition to comments inherited from the original code,
  a number of comments were added through my interpretations.
- <b>Callback interface</b>:
  The library receives function and gradient values via a callback interface.
  The library also notifies the progress of the optimization by invoking a
  callback function. In the original implementation, a user had to set
  function and gradient values every time the function returns for obtaining
  updated values.
- <b>Thread safe</b>:
  The library is thread-safe, which is the secondary gain from the callback
  interface.
- <b>Cross platform.</b> The source code can be compiled on Microsoft Visual
  Studio 2010, GNU C Compiler (gcc), etc.
- <b>Configurable precision</b>: A user can choose single-precision (float)
  or double-precision (double) accuracy by changing ::LBFGS_FLOAT macro.
- <b>SSE/SSE2 optimization</b>:
  This library includes SSE/SSE2 optimization (written in compiler intrinsics)
  for vector arithmetic operations on Intel/AMD processors. The library uses
  SSE for float values and SSE2 for double values. The SSE/SSE2 optimization
  routine is disabled by default.

This library is used by:
- <a href="http://www.chokkan.org/software/crfsuite/">CRFsuite: A fast implementation of Conditional Random Fields (CRFs)</a>
- <a href="http://www.chokkan.org/software/classias/">Classias: A collection of machine-learning algorithms for classification</a>
- <a href="http://www.public.iastate.edu/~gdancik/mlegp/">mlegp: an R package for maximum likelihood estimates for Gaussian processes</a>
- <a href="http://infmath.uibk.ac.at/~matthiasf/imaging2/">imaging2: the imaging2 class library</a>
- <a href="http://search.cpan.org/~laye/Algorithm-LBFGS-0.16/">Algorithm::LBFGS - Perl extension for L-BFGS</a>
- <a href="http://www.cs.kuleuven.be/~bernd/yap-lbfgs/">YAP-LBFGS (an interface to call libLBFGS from YAP Prolog)</a>

@section download Download

- <a href="https://github.com/downloads/chokkan/liblbfgs/liblbfgs-1.10.tar.gz">Source code</a>
- <a href="https://github.com/chokkan/liblbfgs">GitHub repository</a>

libLBFGS is distributed under the term of the
<a href="http://opensource.org/licenses/mit-license.php">MIT license</a>.

@section changelog History
- Version 1.10 (2010-12-22):
    - Fixed compiling errors on Mac OS X; this patch was kindly submitted by
      Nic Schraudolph.
    - Reduced compiling warnings on Mac OS X; this patch was kindly submitted
      by Tamas Nepusz.
    - Replaced memalign() with posix_memalign().
    - Updated solution and project files for Microsoft Visual Studio 2010.
- Version 1.9 (2010-01-29):
    - Fixed a mistake in checking the validity of the parameters "ftol" and
      "wolfe"; this was discovered by Kevin S. Van Horn.
- Version 1.8 (2009-07-13):
    - Accepted the patch submitted by Takashi Imamichi;
      the backtracking method now has three criteria for choosing the step
      length:
        - ::LBFGS_LINESEARCH_BACKTRACKING_ARMIJO: sufficient decrease (Armijo)
          condition only
        - ::LBFGS_LINESEARCH_BACKTRACKING_WOLFE: regular Wolfe condition
          (sufficient decrease condition + curvature condition)
        - ::LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE: strong Wolfe condition
    - Updated the documentation to explain the above three criteria.
- Version 1.7 (2009-02-28):
    - Improved OWL-QN routines for stability.
    - Removed the support of OWL-QN method in MoreThuente algorithm because
      it accidentally fails in early stages of iterations for some objectives.
      Because of this change, <b>the OW-LQN method must be used with the
      backtracking algorithm (::LBFGS_LINESEARCH_BACKTRACKING)</b>, or the
      library returns ::LBFGSERR_INVALID_LINESEARCH.
    - Renamed line search algorithms as follows:
        - ::LBFGS_LINESEARCH_BACKTRACKING: regular Wolfe condition.
        - ::LBFGS_LINESEARCH_BACKTRACKING_LOOSE: regular Wolfe condition.
        - ::LBFGS_LINESEARCH_BACKTRACKING_STRONG: strong Wolfe condition.
    - Source code clean-up.
- Version 1.6 (2008-11-02):
    - Improved line-search algorithm with strong Wolfe condition, which was
      contributed by Takashi Imamichi. This routine is now default for
      ::LBFGS_LINESEARCH_BACKTRACKING. The previous line search algorithm
      with regular Wolfe condition is still available as
      ::LBFGS_LINESEARCH_BACKTRACKING_LOOSE.
    - Configurable stop index for L1-norm computation. A member variable
      ::lbfgs_parameter_t::orthantwise_end was added to specify the index
      number at which the library stops computing the L1 norm of the
      variables. This is useful to prevent some variables from being
      regularized by the OW-LQN method.
    - A sample program written in C++ (sample/sample.cpp).
- Version 1.5 (2008-07-10):
    - Configurable starting index for L1-norm computation. A member variable
      ::lbfgs_parameter_t::orthantwise_start was added to specify the index
      number from which the library computes the L1 norm of the variables.
      This is useful to prevent some variables from being regularized by the
      OWL-QN method.
    - Fixed a zero-division error when the initial variables have already
      been a minimizer (reported by Takashi Imamichi). In this case, the
      library returns ::LBFGS_ALREADY_MINIMIZED status code.
    - Defined ::LBFGS_SUCCESS status code as zero; removed unused constants,
      LBFGSFALSE and LBFGSTRUE.
    - Fixed a compile error in an implicit down-cast.
- Version 1.4 (2008-04-25):
    - Configurable line search algorithms. A member variable
      ::lbfgs_parameter_t::linesearch was added to choose either MoreThuente
      method (::LBFGS_LINESEARCH_MORETHUENTE) or backtracking algorithm
      (::LBFGS_LINESEARCH_BACKTRACKING).
    - Fixed a bug: the previous version did not compute psuedo-gradients
      properly in the line search routines for OWL-QN. This bug might quit
      an iteration process too early when the OWL-QN routine was activated
      (0 < ::lbfgs_parameter_t::orthantwise_c).
    - Configure script for POSIX environments.
    - SSE/SSE2 optimizations with GCC.
    - New functions ::lbfgs_malloc and ::lbfgs_free to use SSE/SSE2 routines
      transparently. It is uncessary to use these functions for libLBFGS built
      without SSE/SSE2 routines; you can still use any memory allocators if
      SSE/SSE2 routines are disabled in libLBFGS.
- Version 1.3 (2007-12-16):
    - An API change. An argument was added to lbfgs() function to receive the
      final value of the objective function. This argument can be set to
      \c NULL if the final value is unnecessary.
    - Fixed a null-pointer bug in the sample code (reported by Takashi Imamichi).
    - Added build scripts for Microsoft Visual Studio 2005 and GCC.
    - Added README file.
- Version 1.2 (2007-12-13):
    - Fixed a serious bug in orthant-wise L-BFGS.
      An important variable was used without initialization.
- Version 1.1 (2007-12-01):
    - Implemented orthant-wise L-BFGS.
    - Implemented lbfgs_parameter_init() function.
    - Fixed several bugs.
    - API documentation.
- Version 1.0 (2007-09-20):
    - Initial release.

@section api Documentation

- @ref liblbfgs_api "libLBFGS API"

@section sample Sample code

@include sample.c

@section ack Acknowledgements

The L-BFGS algorithm is described in:
    - Jorge Nocedal.
      Updating Quasi-Newton Matrices with Limited Storage.
      <i>Mathematics of Computation</i>, Vol. 35, No. 151, pp. 773--782, 1980.
    - Dong C. Liu and Jorge Nocedal.
      On the limited memory BFGS method for large scale optimization.
      <i>Mathematical Programming</i> B, Vol. 45, No. 3, pp. 503-528, 1989.

The line search algorithms used in this implementation are described in:
    - John E. Dennis and Robert B. Schnabel.
      <i>Numerical Methods for Unconstrained Optimization and Nonlinear
      Equations</i>, Englewood Cliffs, 1983.
    - Jorge J. More and David J. Thuente.
      Line search algorithm with guaranteed sufficient decrease.
      <i>ACM Transactions on Mathematical Software (TOMS)</i>, Vol. 20, No. 3,
      pp. 286-307, 1994.

This library also implements Orthant-Wise Limited-memory Quasi-Newton (OWL-QN)
method presented in:
    - Galen Andrew and Jianfeng Gao.
      Scalable training of L1-regularized log-linear models.
      In <i>Proceedings of the 24th International Conference on Machine
      Learning (ICML 2007)</i>, pp. 33-40, 2007.

Special thanks go to:
    - Yoshimasa Tsuruoka and Daisuke Okanohara for technical information about
      OWL-QN
    - Takashi Imamichi for the useful enhancements of the backtracking method
    - Kevin S. Van Horn, Nic Schraudolph, and Tamas Nepusz for bug fixes

Finally I would like to thank the original author, Jorge Nocedal, who has been
distributing the effieicnt and explanatory implementation in an open source
licence.

@section reference Reference

- <a href="http://www.ece.northwestern.edu/~nocedal/lbfgs.html">L-BFGS</a> by Jorge Nocedal.
- <a href="http://research.microsoft.com/en-us/downloads/b1eb1016-1738-4bd5-83a9-370c9d498a03/default.aspx">Orthant-Wise Limited-memory Quasi-Newton Optimizer for L1-regularized Objectives</a> by Galen Andrew.
- <a href="http://chasen.org/~taku/software/misc/lbfgs/">C port (via f2c)</a> by Taku Kudo.
- <a href="http://www.alglib.net/optimization/lbfgs.php">C#/C++/Delphi/VisualBasic6 port</a> in ALGLIB.
- <a href="http://cctbx.sourceforge.net/">Computational Crystallography Toolbox</a> includes
  <a href="http://cctbx.sourceforge.net/current_cvs/c_plus_plus/namespacescitbx_1_1lbfgs.html">scitbx::lbfgs</a>.
*/

#endif/*__LBFGS_H__*/

/*

END OF LBFGS.H Listing

*/
/*

THIS IS THE BEGINNING OF THE LBFGS CODE BLOCK -- USER MUST PROGRAM THE
DERIVATIVES AND THE LOSS FUNCTION BELOW

*/
static lbfgsfloatval_t evaluate(
				void *instance,
				const lbfgsfloatval_t *x,
				lbfgsfloatval_t *g,
				const int n,
				const lbfgsfloatval_t step
			       )
{
	int i, j, k, jj, kk, istop, idebug;
	lbfgsfloatval_t fx = 0.0;
	double sumsquared=0;
	double yeanay, smallx, value, xsum, ysum;
	double distyea, distnay, psiyea, psinay, phiyea, phinay, PHIYEA, PHINAY;
	double INVMILLSYEA, INVMILLSNAY;
	double PI=3.141592653589793;
	double *ZCOORDSD, *XCOORDSD, *GG;
	ZCOORDSD     = calloc( (((nrowX+2*ncolX)*NS)+1), sizeof(double));
	XCOORDSD     = calloc( (((nrowX+2*ncolX)*NS)+1), sizeof(double));
	GG     = calloc((((nrowX+2*ncolX)*NS)+1), sizeof(double));
/*
 */
/*
*
*/
//  TRANSFER ROLL CALL OUTCOME POINTS FROM GRADIENT VECTOR INTO COORDINATE VECTOR
	kk=0;
	for(j=0;j<(NS*2*ncolX);j++)
	{
		GG[j]=0.0;
		if(CONSTRAINTS[j]<1.0){
			ZCOORDSD[j]=0.0;
		}
		if(CONSTRAINTS[j]>0.0){
			ZCOORDSD[j]=x[kk];
			kk=kk+1;
		}
	}
//  TRANSFER LEGISLATOR IDEAL POINTS FROM GRADIENT VECTOR INTO COORDINATE VECTOR
	for(j=0;j<(NS*nrowX);j++)
	{
		GG[j+(NS*2*ncolX)]=0.0;
		if(CONSTRAINTS[j+(NS*2*ncolX)]<1.0){
			XCOORDSD[j]=0.0;
		}
		if(CONSTRAINTS[j+(NS*2*ncolX)]>0.0){
			XCOORDSD[j]=x[kk];
			kk=kk+1;
		}
	}
//	kk=kk+1;
//	BETA=x[kk];
//	GG[((nrowX+2*ncolX)*NS)+1]=0.0;
	ysum=-99999.0;
	for (i=0;i<nrowX;i++)
	{
		xsum=0.0;
		for(k=0;k<NS;k++)
		{
//  NORMALIZE TO UNIT HYPERSPHERE
			xsum=xsum+XCOORDSD[NS*i+k]*XCOORDSD[NS*i+k];
		}
		if(xsum > ysum)
		{
			ysum=xsum;
		}
	}
//
	for(i=0;i<nrowX;i++)
	{
		for(j=0;j<jne;j++)
		{
			fprintf(jp,"%c",kname[(i)*jne+j]);
		}
		for(j=0;j<NS;j++)
		{
//			XCOORDSD[i*NS+j]=XCOORDSD[i*NS+j]*(1.0/sqrt(ysum));
			fprintf(jp,"%10.6f",XCOORDSD[i*NS+j]);
		}
		fprintf(jp,"\n");
	}
	for(i=0;i<ncolX;i++)
	{
		for(j=0;j<NS;j++)
		{
//			ZCOORDSD[2*i*NS+j]=ZCOORDSD[2*i*NS+j]*(1.0/sqrt(ysum));
			if(i < 10)fprintf(jp,"%5d %10.6f",i,ZCOORDSD[2*i*NS+j]);
//			ZCOORDSD[2*i*NS+j+NS]=ZCOORDSD[2*i*NS+j+NS]*(1.0/sqrt(ysum));
			if(i < 10)fprintf(jp,"%5d %10.6f",i,ZCOORDSD[2*i*NS+j+NS]);
		}
		if(i < 10)fprintf(jp,"\n");
	}
//
	for(i=0;i<N;i++)
	{
//		x[i]=x[i]*(1.0/sqrt(ysum));
	}
//
	for(i=0;i<nrowX;i++)
	{
		for(j=0;j<ncolX;j++)
		{
			yeanay = X[i*ncolX+j];
/*
CATCH MISSING DATA HERE
*/
			if(yeanay > 0.0)
			{
				// VOTE YEA
				if(yeanay < 6.0){
				}
				// VOTE NAY
				if(yeanay > 1.0){
				}
				distyea=0.0;
				distnay=0.0;
					for(jj=0;jj<NS;jj++)
					{
						distyea = distyea+pow((XCOORDSD[NS*i+jj]-ZCOORDSD[2*NS*j+jj]),2.0);
						distnay = distnay+pow((XCOORDSD[NS*i+jj]-ZCOORDSD[2*NS*j+jj+NS]),2.0);
					}
					psiyea=exp((-.5)*distyea) - exp((-.5)*distnay);
					psinay=exp((-.5)*distnay) - exp((-.5)*distyea);
					smallx = (BETA*psiyea)/sqrt(2.0);
					value = erf(smallx);       
					PHIYEA = (value/2.0) + 0.5;
					smallx = (BETA*psinay)/sqrt(2.0);
					value = erf(smallx);       
					PHINAY = (value/2.0) + 0.5;
					phiyea = (1/sqrt(2*PI))*exp((-.5)*BETA*BETA*psiyea*psiyea);
					phinay = (1/sqrt(2*PI))*exp((-.5)*BETA*BETA*psinay*psinay);
					if(PHIYEA < .0001){
						PHIYEA=.001;
						PHINAY=.999;
					}
					INVMILLSYEA = phiyea/PHIYEA;
					if(PHINAY < .0001){
						PHINAY=.001;
						PHIYEA=.999;
					}
					INVMILLSNAY = phinay/PHINAY;
					fprintf(jp,"L-BFGS %5d %5d %2.0f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
						  i,j,yeanay,distyea,distnay,psiyea,psinay,phiyea,phinay,PHIYEA,PHINAY,INVMILLSYEA,INVMILLSNAY,BETA);
//  CALCULATE LOG OF LIKELIHOOD FUNCTION
				// VOTE YEA
					if(yeanay < 6.0){
						sumsquared=sumsquared+log(PHIYEA);
					}
				// VOTE NAY
					if(yeanay > 1.0){
						sumsquared=sumsquared+log(PHINAY);
					}

/* DERIVATIVES FOR LEGISLATORS
 * */
					for(jj=0;jj<NS;jj++)
					{
				// VOTE YEA
						if(yeanay < 6.0){
							GG[NS*i+jj+2*NS*ncolX]=GG[NS*i+jj+2*NS*ncolX]-INVMILLSYEA*
							    (((-1.0)*(exp((-.5)*distyea))*(XCOORDSD[NS*i+jj]-ZCOORDSD[2*NS*j+jj]))+
								    ((exp((-.5)*distnay))*(XCOORDSD[NS*i+jj]-ZCOORDSD[2*NS*j+jj+NS])));
						}
				// VOTE NAY
						if(yeanay > 1.0){
							GG[NS*i+jj+2*NS*ncolX]=GG[NS*i+jj+2*NS*ncolX]-INVMILLSNAY*
							    (((-1.0)*(exp((-.5)*distnay))*(XCOORDSD[NS*i+jj]-ZCOORDSD[2*NS*j+jj+NS]))+
								    ((exp((-.5)*distyea))*(XCOORDSD[NS*i+jj]-ZCOORDSD[2*NS*j+jj])));
								
						}
					}
/* DERIVATIVES FOR ROLL CALLS
 * */
					for(jj=0;jj<NS;jj++)
					{
				// VOTE YEA
						if(yeanay < 6.0){
							GG[2*NS*j+jj]=GG[2*NS*j+jj]-INVMILLSYEA*
								((exp((-.5)*distyea))*(XCOORDSD[NS*i+jj]-ZCOORDSD[2*NS*j+jj]));
							GG[2*NS*j+jj+NS]=GG[2*NS*j+jj+NS]-INVMILLSNAY*
								((-1.0)*(exp((-.5)*distnay))*(XCOORDSD[NS*i+jj]-ZCOORDSD[2*NS*j+jj+NS]));								
						}
				// VOTE NAY
						if(yeanay > 1.0){
							GG[2*NS*j+jj]=GG[2*NS*j+jj]-INVMILLSYEA*
								((-1.0)*(exp((-.5)*distyea))*(XCOORDSD[NS*i+jj]-ZCOORDSD[2*NS*j+jj]));
							GG[2*NS*j+jj+NS]=GG[2*NS*j+jj+NS]-INVMILLSNAY*
								((exp((-.5)*distnay))*(XCOORDSD[NS*i+jj]-ZCOORDSD[2*NS*j+jj+NS]));								
						}
					}
			}
		}
	}
	fx=-sumsquared;
	fprintf(kp,"LOSS FUNCTION L-BFGS %5d %5d %5d %5d %10.6f %20.6f\n",n,N,nrowX,ncolX,BETA,fx);
        printf("LOSS FUNCTION L-BFGS %5d %5d %5d %5d %10.6f %20.6f\n",n,N,nrowX,ncolX,BETA,fx);
//
//  TRANSFER INTO GRADIENT VECTOR
//
	kk=0;
	for(j=0;j<(ncolX*2*NS);j++)
	{
		if(CONSTRAINTS[j]>0.0){
			g[kk]=GG[j];
			kk=kk+1;
		}
	}
	for(j=0;j<(nrowX*NS);j++)
	{
		if(CONSTRAINTS[j+(ncolX*NS)]>0.0){
			g[kk]=GG[j+(ncolX*2*NS)];
			kk=kk+1;
		}
	}
//	g[kk]=GG[((nrowX+2*ncolX)*NS)];
//	kk=kk+1;
//	GG[((nrowX+2*ncolX)*NS)+1]
	for (i = 0;i < kk;i++) {
//
//	    fprintf(kp,"INITIAL VALUES GRADIENT %5d %12.6f\n",i,g[i]);
//	    printf("INITIAL VALUES GRADIENT %5d %12.6f\n",i,g[i]);
	}

//
//    SIMPLE METHOD OF STOPPING EXECUTION
//
//	istop=100;
//	if(istop==100)idebug=99998;
//	if(idebug==99998)exit(EXIT_FAILURE);
//

        free(ZCOORDSD);
        free(XCOORDSD);
	free(GG);
	return fx;
}

static int progress(
		    void *instance,
		    const lbfgsfloatval_t *x,
		    const lbfgsfloatval_t *g,
		    const lbfgsfloatval_t fx,
		    const lbfgsfloatval_t xnorm,
		    const lbfgsfloatval_t gnorm,
		    const lbfgsfloatval_t step,
		    int n,
		    int k,
		    int ls
		   )
{
	printf("Iteration %d:\n", k);
	printf("  fx = %f, x[0] = %f,  x[1] = %f,  x[2] = %f, x[3] = %f\n", fx, x[0], x[1], x[2], x[3]);
	printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
	printf("\n");
//
	fprintf(kp,"Iteration %d:\n", k);
	fprintf(kp,"  fx = %f, x[0] = %f,  x[1] = %f,  x[2] = %f, x[3] = %f\n", fx, x[0], x[1], x[2], x[3]);
	fprintf(kp,"  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
	fprintf(kp,"\n");
	return 0;
}

/*
 *
 *  THIS IS THE TOP-LEVEL ROUTINE -- the "MAIN" Code for LBFGS -- AGAIN
 *  THE USER MUST PROGRAM THE INITIAL VALUES WHICH ARE FED TO THE
 *  ROUTINES ABOVE -- "evaluate" and "progress"
 *
 *kpnp = number of rows of y
 *kpnq = number of columns of y
rmatrix holds the starting coordinates from the double-centering
yrotate returns the results from L-BFGS (Limited-Memory
       Broyden-Fletcher-Goldfarb-Shanno algorthm
y is the symmetric matrix of dissimilarities (distances)
*/

void mainlbfgs(int kpnp, int kpnq, double *yrotate, double *rmatrix)
{
    double *a, *b;
    int i, j, istop, idebug, kk, ret = 0;
    lbfgsfloatval_t fx;
    lbfgsfloatval_t *x = lbfgs_malloc(N);
    lbfgs_parameter_t param;
    a     = calloc( (((kpnp+2*kpnq)*NS)+1), sizeof(double));
    b     = calloc( (((kpnp+2*kpnq)*NS)+1), sizeof(double));
// Transfer coordinates
    for(i=0;i<((2*kpnq+kpnp)*NS);i++)
    {
	    a[i]=rmatrix[i];
//	    printf("ZCOORDS3 VALUES %5d %12.6f\n",i,a[i]);
    }
    if (x == NULL) {
        printf("ERROR: Failed to allocate a memory block for variables.\n");
//	return 1;
	return;
    }
    kk=0;
    for(j=0;j<(kpnq*NS*2);j++)
    {
	    if(CONSTRAINTS[j]<1.0){
//		    XCOORDS[j]=0.0;
	    }
	    if(CONSTRAINTS[j]>0.0){
		    x[kk]=a[j];
		    kk=kk+1;
	    }
    }
    for(j=0;j<(kpnp*NS);j++)
    {
	    if(CONSTRAINTS[j+(kpnq*NS*2)]>0.0){
		    x[kk]=a[j+(kpnq*NS*2)];
		    kk=kk+1;
	    }
    }
//    x[kk]=a[((kpnp+2*kpnq)*NS)];
//    kk=kk+1;
//
    /* Initialize the variables. */
    for (i = 0;i < kk;i++) {
//
	    fprintf(kp,"INITIAL VALUES %5d %12.6f\n",i,x[i]);
	    printf("INITIAL VALUES %5d %12.6f\n",i,x[i]);
    }

    /* Initialize the parameters for the L-BFGS optimization. */
    lbfgs_parameter_init(&param);
    /*param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;*/

    /*
        Start the L-BFGS optimization; this will invoke the callback functions
        evaluate() and progress() when necessary.
     */

//    ret = lbfgs(*N, x, &fx, evaluate, progress, NULL, &param);
    ret = lbfgs(kk, x, &fx, evaluate, progress, NULL, &param);

    /* Report the result. */
    printf("L-BFGS optimization terminated with status code = %d\n", ret);
    fprintf(kp,"L-BFGS optimization terminated with status code = %d\n", ret);
//    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  fx = %f, x[0] = %f,  x[1] = %f,  x[2] = %f, x[3] = %f\n", fx, x[0], x[1], x[2], x[3]);
    fprintf(kp,"  fx = %f, x[0] = %f,  x[1] = %f,  x[2] = %f, x[3] = %f\n", fx, x[0], x[1], x[2], x[3]);

// PASS BACK THE SOLUTION
//    
    kk=0;
    for(j=0;j<(kpnq*NS);j++)
    {
	    if(CONSTRAINTS[j]<1.0){
		    yrotate[j]=0.0;
	    }
	    if(CONSTRAINTS[j]>0.0){
		    yrotate[j]=x[kk];
		    kk=kk+1;
	    }
    }
    for(j=0;j<(kpnp*NS);j++)
    {
//	    yrotate[j+(kpnq*NS)]=a[j+(kpnq*NS)];
//	    if(a[j+(kpnq*NS)]> -99.0){
	    if(CONSTRAINTS[j+(kpnq*NS)]<1.0){
		    yrotate[j]=0.0;
	    }
		    if(CONSTRAINTS[j+(kpnq*NS)]>0.0){
		    yrotate[j+(kpnq*NS)]=x[kk];
		    kk=kk+1;
	    }
    }
//
    lbfgs_free(x);
    free(a);
    free(b);
    return;
}
/*

BELOW ARE THE ROUTINES THAT DO THE OPTIMIZATION USING THE FUNCTION
EVALUATION ROUTINES ABOVE

*/
/*
 *      Limited memory BFGS (L-BFGS).
 *
 * Copyright (c) 1990, Jorge Nocedal
 * Copyright (c) 2007-2010 Naoaki Okazaki
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* $Id$ */

/*
This library is a C port of the FORTRAN implementation of Limited-memory
Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) method written by Jorge Nocedal.
The original FORTRAN source code is available at:
http://www.ece.northwestern.edu/~nocedal/lbfgs.html

The L-BFGS algorithm is described in:
    - Jorge Nocedal.
      Updating Quasi-Newton Matrices with Limited Storage.
      <i>Mathematics of Computation</i>, Vol. 35, No. 151, pp. 773--782, 1980.
    - Dong C. Liu and Jorge Nocedal.
      On the limited memory BFGS method for large scale optimization.
      <i>Mathematical Programming</i> B, Vol. 45, No. 3, pp. 503-528, 1989.

The line search algorithms used in this implementation are described in:
    - John E. Dennis and Robert B. Schnabel.
      <i>Numerical Methods for Unconstrained Optimization and Nonlinear
      Equations</i>, Englewood Cliffs, 1983.
    - Jorge J. More and David J. Thuente.
      Line search algorithm with guaranteed sufficient decrease.
      <i>ACM Transactions on Mathematical Software (TOMS)</i>, Vol. 20, No. 3,
      pp. 286-307, 1994.

This library also implements Orthant-Wise Limited-memory Quasi-Newton (OWL-QN)
method presented in:
    - Galen Andrew and Jianfeng Gao.
      Scalable training of L1-regularized log-linear models.
      In <i>Proceedings of the 24th International Conference on Machine
      Learning (ICML 2007)</i>, pp. 33-40, 2007.

I would like to thank the original author, Jorge Nocedal, who has been
distributing the effieicnt and explanatory implementation in an open source
licence.
*/

#ifdef  HAVE_CONFIG_H
#include <config.h>
#endif/*HAVE_CONFIG_H*/

//#include <stdint.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>

//#include <lbfgs.h>

#ifdef  _MSC_VER
#define inline  __inline
#endif/*_MSC_VER*/

#if     defined(USE_SSE) && defined(__SSE2__) && LBFGS_FLOAT == 64
/* Use SSE2 optimization for 64bit double precision. */
#include "arithmetic_sse_double.h"

#elseif   defined(USE_SSE) && defined(__SSE__) && LBFGS_FLOAT == 32
/* Use SSE optimization for 32bit float precision. */
#include "arithmetic_sse_float.h"

#else
/* No CPU specific optimization. */
//#include "arithmetic_ansi.h"

#endif

/*

END OF IN-LINE LBFGS ROUTINES (i.e., not in the Header Files)

*/
/*
 *   ARITHMETIC_ANSI.H listing Below
 *
*/
/*
 *      ANSI C implementation of vector operations.
 *
 * Copyright (c) 2007-2010 Naoaki Okazaki
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* $Id$ */

//#include <stdlib.h>
//#include <memory.h>

#if     LBFGS_FLOAT == 32 && LBFGS_IEEE_FLOAT
#define fsigndiff(x, y) (((*(uint32_t*)(x)) ^ (*(uint32_t*)(y))) & 0x80000000U)
#else
#define fsigndiff(x, y) (*(x) * (*(y) / fabs(*(y))) < 0.)
#endif/*LBFGS_IEEE_FLOAT*/

inline static void* vecalloc(size_t size)
{
	void *memblock = malloc(size);
	if (memblock) {
		memset(memblock, 0, size);
	}
	return memblock;
}

inline static void vecfree(void *memblock)
{
	free(memblock);
}

inline static void vecset(lbfgsfloatval_t *x, const lbfgsfloatval_t c, const int n)
{
	int i;

	for (i = 0;i < n;++i) {
		x[i] = c;
	}
}

inline static void veccpy(lbfgsfloatval_t *y, const lbfgsfloatval_t *x, const int n)
{
	int i;

	for (i = 0;i < n;++i) {
		y[i] = x[i];
	}
}

inline static void vecncpy(lbfgsfloatval_t *y, const lbfgsfloatval_t *x, const int n)
{
	int i;

	for (i = 0;i < n;++i) {
		y[i] = -x[i];
	}
}

inline static void vecadd(lbfgsfloatval_t *y, const lbfgsfloatval_t *x, const lbfgsfloatval_t c, const int n)
{
	int i;

	for (i = 0;i < n;++i) {
		y[i] += c * x[i];
	}
}

inline static void vecdiff(lbfgsfloatval_t *z, const lbfgsfloatval_t *x, const lbfgsfloatval_t *y, const int n)
{
	int i;

	for (i = 0;i < n;++i) {
		z[i] = x[i] - y[i];
	}
}

inline static void vecscale(lbfgsfloatval_t *y, const lbfgsfloatval_t c, const int n)
{
	int i;

	for (i = 0;i < n;++i) {
		y[i] *= c;
	}
}

inline static void vecmul(lbfgsfloatval_t *y, const lbfgsfloatval_t *x, const int n)
{
	int i;

	for (i = 0;i < n;++i) {
		y[i] *= x[i];
	}
}

inline static void vecdot(lbfgsfloatval_t* s, const lbfgsfloatval_t *x, const lbfgsfloatval_t *y, const int n)
{
	int i;
	*s = 0.;
	for (i = 0;i < n;++i) {
		*s += x[i] * y[i];
	}
}

inline static void vec2norm(lbfgsfloatval_t* s, const lbfgsfloatval_t *x, const int n)
{
	vecdot(s, x, x, n);
	*s = (lbfgsfloatval_t)sqrt(*s);
}

inline static void vec2norminv(lbfgsfloatval_t* s, const lbfgsfloatval_t *x, const int n)
{
	vec2norm(s, x, n);
	*s = (lbfgsfloatval_t)(1.0 / *s);
}
/*

END OF ARITHMETIC_ANSI.H CODE BLOCK

*/
//#endif
/*

MORE IN-LINE LBFGS CODE (Not in the Header Files)

*/

#define min2(a, b)      ((a) <= (b) ? (a) : (b))
#define max2(a, b)      ((a) >= (b) ? (a) : (b))
#define max3(a, b, c)   max2(max2((a), (b)), (c));

struct tag_callback_data {
    int n;
    void *instance;
    lbfgs_evaluate_t proc_evaluate;
    lbfgs_progress_t proc_progress;
};
typedef struct tag_callback_data callback_data_t;

struct tag_iteration_data {
    lbfgsfloatval_t alpha;
    lbfgsfloatval_t *s;     /* [n] */
    lbfgsfloatval_t *y;     /* [n] */
    lbfgsfloatval_t ys;     /* vecdot(y, s) */
};
typedef struct tag_iteration_data iteration_data_t;
// STOPPING CRITERIA -- EPSILON IS THE SECOND NUMBER --
// GNORM/XNORM < epsilon IS
// THE STOPPING CRITERIA
/*
Variables: m=6, number of corrections to approximate the inverse
                hessian matrix
           epsilon = .002
           past = 0, distance for delta-based convergence test
           delta=1e-5, the minimum rate of decrease of the objective function.
           max_iterations=0, keeps iterating until convergence
           LBFGS_LINESEARCH_DEFAULT = line search algorithm
           max_linesearch=40, number of trials for the line search
           min_step=1e-20, minimum step of the line search routine
           max_step=1e+20, maximum number trials for the line search
           ftol=1e-4, accuracy of line search routine
           gtol=.9, accuracy of line search routine (2nd parameter)
           wolfe=.9, only used if backtracking is the Wolfe condition
           xtol=1.0e-16, the machine precision for floating point values
           orthanwise_c=0.0, coefficient for the L1 norm of variables
           orthantwise_start=0, Start index for computing L1 norm
           orthantwise_end=-1, End index for computing L1 norm
           
*/
static const lbfgs_parameter_t _defparam = {
//	6, 1e-1, 0, 1e-5,
	6, 0.5, 0, 0.001,
    0, LBFGS_LINESEARCH_DEFAULT, 40,
    1e-20, 1e20, 1e-4, 0.9, 0.9, 1.0e-16,
    0.0, 0, -1,
};

/* Forward function declarations. */

typedef int (*line_search_proc)(
    int n,
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *f,
    lbfgsfloatval_t *g,
    lbfgsfloatval_t *s,
    lbfgsfloatval_t *stp,
    const lbfgsfloatval_t* xp,
    const lbfgsfloatval_t* gp,
    lbfgsfloatval_t *wa,
    callback_data_t *cd,
    const lbfgs_parameter_t *param
    );
    
static int line_search_backtracking(
    int n,
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *f,
    lbfgsfloatval_t *g,
    lbfgsfloatval_t *s,
    lbfgsfloatval_t *stp,
    const lbfgsfloatval_t* xp,
    const lbfgsfloatval_t* gp,
    lbfgsfloatval_t *wa,
    callback_data_t *cd,
    const lbfgs_parameter_t *param
    );

static int line_search_backtracking_owlqn(
    int n,
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *f,
    lbfgsfloatval_t *g,
    lbfgsfloatval_t *s,
    lbfgsfloatval_t *stp,
    const lbfgsfloatval_t* xp,
    const lbfgsfloatval_t* gp,
    lbfgsfloatval_t *wp,
    callback_data_t *cd,
    const lbfgs_parameter_t *param
    );

static int line_search_morethuente(
    int n,
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *f,
    lbfgsfloatval_t *g,
    lbfgsfloatval_t *s,
    lbfgsfloatval_t *stp,
    const lbfgsfloatval_t* xp,
    const lbfgsfloatval_t* gp,
    lbfgsfloatval_t *wa,
    callback_data_t *cd,
    const lbfgs_parameter_t *param
    );

static int update_trial_interval(
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *fx,
    lbfgsfloatval_t *dx,
    lbfgsfloatval_t *y,
    lbfgsfloatval_t *fy,
    lbfgsfloatval_t *dy,
    lbfgsfloatval_t *t,
    lbfgsfloatval_t *ft,
    lbfgsfloatval_t *dt,
    const lbfgsfloatval_t tmin,
    const lbfgsfloatval_t tmax,
    int *brackt
    );

static lbfgsfloatval_t owlqn_x1norm(
    const lbfgsfloatval_t* x,
    const int start,
    const int n
    );

static void owlqn_pseudo_gradient(
    lbfgsfloatval_t* pg,
    const lbfgsfloatval_t* x,
    const lbfgsfloatval_t* g,
    const int n,
    const lbfgsfloatval_t c,
    const int start,
    const int end
    );

static void owlqn_project(
    lbfgsfloatval_t* d,
    const lbfgsfloatval_t* sign,
    const int start,
    const int end
    );


#if     defined(USE_SSE) && (defined(__SSE__) || defined(__SSE2__))
static int round_out_variables(int n)
{
    n += 7;
    n /= 8;
    n *= 8;
    return n;
}
#endif/*defined(USE_SSE)*/

lbfgsfloatval_t* lbfgs_malloc(int n)
{
#if     defined(USE_SSE) && (defined(__SSE__) || defined(__SSE2__))
    n = round_out_variables(n);
#endif/*defined(USE_SSE)*/
    return (lbfgsfloatval_t*)vecalloc(sizeof(lbfgsfloatval_t) * n);
}

void lbfgs_free(lbfgsfloatval_t *x)
{
    vecfree(x);
}

void lbfgs_parameter_init(lbfgs_parameter_t *param)
{
    memcpy(param, &_defparam, sizeof(*param));
}

int lbfgs(
    int n,
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *ptr_fx,
    lbfgs_evaluate_t proc_evaluate,
    lbfgs_progress_t proc_progress,
    void *instance,
    lbfgs_parameter_t *_param
    )
{
    int ret;
    int i, j, k, ls, end, bound;
    lbfgsfloatval_t step;

    /* Constant parameters and their default values. */
    lbfgs_parameter_t param = (_param != NULL) ? (*_param) : _defparam;
    const int m = param.m;

    lbfgsfloatval_t *xp = NULL;
    lbfgsfloatval_t *g = NULL, *gp = NULL, *pg = NULL;
    lbfgsfloatval_t *d = NULL, *w = NULL, *pf = NULL;
    iteration_data_t *lm = NULL, *it = NULL;
    lbfgsfloatval_t ys, yy;
    lbfgsfloatval_t xnorm, gnorm, beta;
    lbfgsfloatval_t fx = 0.;
    lbfgsfloatval_t rate = 0.;
    line_search_proc linesearch = line_search_morethuente;

    /* Construct a callback data. */
    callback_data_t cd;
    cd.n = n;
    cd.instance = instance;
    cd.proc_evaluate = proc_evaluate;
    cd.proc_progress = proc_progress;

#if     defined(USE_SSE) && (defined(__SSE__) || defined(__SSE2__))
    /* Round out the number of variables. */
    n = round_out_variables(n);
#endif/*defined(USE_SSE)*/

    /* Check the input parameters for errors. */
    if (n <= 0) {
        return LBFGSERR_INVALID_N;
    }
#if     defined(USE_SSE) && (defined(__SSE__) || defined(__SSE2__))
    if (n % 8 != 0) {
        return LBFGSERR_INVALID_N_SSE;
    }
    if ((uintptr_t)(const void*)x % 16 != 0) {
        return LBFGSERR_INVALID_X_SSE;
    }
#endif/*defined(USE_SSE)*/
    if (param.epsilon < 0.) {
        return LBFGSERR_INVALID_EPSILON;
    }
    if (param.past < 0) {
        return LBFGSERR_INVALID_TESTPERIOD;
    }
    if (param.delta < 0.) {
        return LBFGSERR_INVALID_DELTA;
    }
    if (param.min_step < 0.) {
        return LBFGSERR_INVALID_MINSTEP;
    }
    if (param.max_step < param.min_step) {
        return LBFGSERR_INVALID_MAXSTEP;
    }
    if (param.ftol < 0.) {
        return LBFGSERR_INVALID_FTOL;
    }
    if (param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_WOLFE ||
        param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE) {
        if (param.wolfe <= param.ftol || 1. <= param.wolfe) {
            return LBFGSERR_INVALID_WOLFE;
        }
    }
    if (param.gtol < 0.) {
        return LBFGSERR_INVALID_GTOL;
    }
    if (param.xtol < 0.) {
        return LBFGSERR_INVALID_XTOL;
    }
    if (param.max_linesearch <= 0) {
        return LBFGSERR_INVALID_MAXLINESEARCH;
    }
    if (param.orthantwise_c < 0.) {
        return LBFGSERR_INVALID_ORTHANTWISE;
    }
    if (param.orthantwise_start < 0 || n < param.orthantwise_start) {
        return LBFGSERR_INVALID_ORTHANTWISE_START;
    }
    if (param.orthantwise_end < 0) {
        param.orthantwise_end = n;
    }
    if (n < param.orthantwise_end) {
        return LBFGSERR_INVALID_ORTHANTWISE_END;
    }
    if (param.orthantwise_c != 0.) {
        switch (param.linesearch) {
        case LBFGS_LINESEARCH_BACKTRACKING:
            linesearch = line_search_backtracking_owlqn;
            break;
        default:
            /* Only the backtracking method is available. */
            return LBFGSERR_INVALID_LINESEARCH;
        }
    } else {
        switch (param.linesearch) {
        case LBFGS_LINESEARCH_MORETHUENTE:
            linesearch = line_search_morethuente;
            break;
        case LBFGS_LINESEARCH_BACKTRACKING_ARMIJO:
        case LBFGS_LINESEARCH_BACKTRACKING_WOLFE:
        case LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE:
            linesearch = line_search_backtracking;
            break;
        default:
            return LBFGSERR_INVALID_LINESEARCH;
        }
    }

    /* Allocate working space. */
    xp = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
    g = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
    gp = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
    d = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
    w = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
    if (xp == NULL || g == NULL || gp == NULL || d == NULL || w == NULL) {
        ret = LBFGSERR_OUTOFMEMORY;
        goto lbfgs_exit;
    }

    if (param.orthantwise_c != 0.) {
        /* Allocate working space for OW-LQN. */
        pg = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
        if (pg == NULL) {
            ret = LBFGSERR_OUTOFMEMORY;
            goto lbfgs_exit;
        }
    }

    /* Allocate limited memory storage. */
    lm = (iteration_data_t*)vecalloc(m * sizeof(iteration_data_t));
    if (lm == NULL) {
        ret = LBFGSERR_OUTOFMEMORY;
        goto lbfgs_exit;
    }

    /* Initialize the limited memory. */
    for (i = 0;i < m;++i) {
        it = &lm[i];
        it->alpha = 0;
        it->ys = 0;
        it->s = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
        it->y = (lbfgsfloatval_t*)vecalloc(n * sizeof(lbfgsfloatval_t));
        if (it->s == NULL || it->y == NULL) {
            ret = LBFGSERR_OUTOFMEMORY;
            goto lbfgs_exit;
        }
    }

    /* Allocate an array for storing previous values of the objective function. */
    if (0 < param.past) {
        pf = (lbfgsfloatval_t*)vecalloc(param.past * sizeof(lbfgsfloatval_t));
    }

    /* Evaluate the function value and its gradient. */
    fx = cd.proc_evaluate(cd.instance, x, g, cd.n, 0);
    if (0. != param.orthantwise_c) {
        /* Compute the L1 norm of the variable and add it to the object value. */
        xnorm = owlqn_x1norm(x, param.orthantwise_start, param.orthantwise_end);
        fx += xnorm * param.orthantwise_c;
        owlqn_pseudo_gradient(
            pg, x, g, n,
            param.orthantwise_c, param.orthantwise_start, param.orthantwise_end
            );
    }

    /* Store the initial value of the objective function. */
    if (pf != NULL) {
        pf[0] = fx;
    }

    /*
        Compute the direction;
        we assume the initial hessian matrix H_0 as the identity matrix.
     */
    if (param.orthantwise_c == 0.) {
        vecncpy(d, g, n);
    } else {
        vecncpy(d, pg, n);
    }

    /*
       Make sure that the initial variables are not a minimizer.
     */
    vec2norm(&xnorm, x, n);
    if (param.orthantwise_c == 0.) {
        vec2norm(&gnorm, g, n);
    } else {
        vec2norm(&gnorm, pg, n);
    }
    if (xnorm < 1.0) xnorm = 1.0;
    if (gnorm / xnorm <= param.epsilon) {
        ret = LBFGS_ALREADY_MINIMIZED;
        goto lbfgs_exit;
    }

    /* Compute the initial step:
        step = 1.0 / sqrt(vecdot(d, d, n))
     */
    vec2norminv(&step, d, n);

    k = 1;
    end = 0;
    for (;;) {
        /* Store the current position and gradient vectors. */
        veccpy(xp, x, n);
        veccpy(gp, g, n);

        /* Search for an optimal step. */
        if (param.orthantwise_c == 0.) {
            ls = linesearch(n, x, &fx, g, d, &step, xp, gp, w, &cd, &param);
        } else {
            ls = linesearch(n, x, &fx, g, d, &step, xp, pg, w, &cd, &param);
            owlqn_pseudo_gradient(
                pg, x, g, n,
                param.orthantwise_c, param.orthantwise_start, param.orthantwise_end
                );
        }
        if (ls < 0) {
            /* Revert to the previous point. */
            veccpy(x, xp, n);
            veccpy(g, gp, n);
            ret = ls;
            goto lbfgs_exit;
        }

        /* Compute x and g norms. */
        vec2norm(&xnorm, x, n);
        if (param.orthantwise_c == 0.) {
            vec2norm(&gnorm, g, n);
        } else {
            vec2norm(&gnorm, pg, n);
        }

        /* Report the progress. */
        if (cd.proc_progress) {
            if ((ret = cd.proc_progress(cd.instance, x, g, fx, xnorm, gnorm, step, cd.n, k, ls))) {
                goto lbfgs_exit;
            }
        }

        /*
            Convergence test.
            The criterion is given by the following formula:
                |g(x)| / \max(1, |x|) < \epsilon
         */
        if (xnorm < 1.0) xnorm = 1.0;
        if (gnorm / xnorm <= param.epsilon) {
            /* Convergence. */
            ret = LBFGS_SUCCESS;
            break;
        }

        /*
            Test for stopping criterion.
            The criterion is given by the following formula:
                (f(past_x) - f(x)) / f(x) < \delta
         */
        if (pf != NULL) {
            /* We don't test the stopping criterion while k < past. */
            if (param.past <= k) {
                /* Compute the relative improvement from the past. */
                rate = (pf[k % param.past] - fx) / fx;

                /* The stopping criterion. */
                if (rate < param.delta) {
                    ret = LBFGS_STOP;
                    break;
                }
            }

            /* Store the current value of the objective function. */
            pf[k % param.past] = fx;
        }

        if (param.max_iterations != 0 && param.max_iterations < k+1) {
            /* Maximum number of iterations. */
            ret = LBFGSERR_MAXIMUMITERATION;
            break;
        }

        /*
            Update vectors s and y:
                s_{k+1} = x_{k+1} - x_{k} = \step * d_{k}.
                y_{k+1} = g_{k+1} - g_{k}.
         */
        it = &lm[end];
        vecdiff(it->s, x, xp, n);
        vecdiff(it->y, g, gp, n);

        /*
            Compute scalars ys and yy:
                ys = y^t \cdot s = 1 / \rho.
                yy = y^t \cdot y.
            Notice that yy is used for scaling the hessian matrix H_0 (Cholesky factor).
         */
        vecdot(&ys, it->y, it->s, n);
        vecdot(&yy, it->y, it->y, n);
        it->ys = ys;

        /*
            Recursive formula to compute dir = -(H \cdot g).
                This is described in page 779 of:
                Jorge Nocedal.
                Updating Quasi-Newton Matrices with Limited Storage.
                Mathematics of Computation, Vol. 35, No. 151,
                pp. 773--782, 1980.
         */
        bound = (m <= k) ? m : k;
        ++k;
        end = (end + 1) % m;

        /* Compute the steepest direction. */
        if (param.orthantwise_c == 0.) {
            /* Compute the negative of gradients. */
            vecncpy(d, g, n);
        } else {
            vecncpy(d, pg, n);
        }

        j = end;
        for (i = 0;i < bound;++i) {
            j = (j + m - 1) % m;    /* if (--j == -1) j = m-1; */
            it = &lm[j];
            /* \alpha_{j} = \rho_{j} s^{t}_{j} \cdot q_{k+1}. */
            vecdot(&it->alpha, it->s, d, n);
            it->alpha /= it->ys;
            /* q_{i} = q_{i+1} - \alpha_{i} y_{i}. */
            vecadd(d, it->y, -it->alpha, n);
        }

        vecscale(d, ys / yy, n);

        for (i = 0;i < bound;++i) {
            it = &lm[j];
            /* \beta_{j} = \rho_{j} y^t_{j} \cdot \gamma_{i}. */
            vecdot(&beta, it->y, d, n);
            beta /= it->ys;
            /* \gamma_{i+1} = \gamma_{i} + (\alpha_{j} - \beta_{j}) s_{j}. */
            vecadd(d, it->s, it->alpha - beta, n);
            j = (j + 1) % m;        /* if (++j == m) j = 0; */
        }

        /*
            Constrain the search direction for orthant-wise updates.
         */
        if (param.orthantwise_c != 0.) {
            for (i = param.orthantwise_start;i < param.orthantwise_end;++i) {
                if (d[i] * pg[i] >= 0) {
                    d[i] = 0;
                }
            }
        }

        /*
            Now the search direction d is ready. We try step = 1 first.
         */
        step = 1.0;
    }

lbfgs_exit:
    /* Return the final value of the objective function. */
    if (ptr_fx != NULL) {
        *ptr_fx = fx;
    }

    vecfree(pf);

    /* Free memory blocks used by this function. */
    if (lm != NULL) {
        for (i = 0;i < m;++i) {
            vecfree(lm[i].s);
            vecfree(lm[i].y);
        }
        vecfree(lm);
    }
    vecfree(pg);
    vecfree(w);
    vecfree(d);
    vecfree(gp);
    vecfree(g);
    vecfree(xp);

    return ret;
}



static int line_search_backtracking(
    int n,
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *f,
    lbfgsfloatval_t *g,
    lbfgsfloatval_t *s,
    lbfgsfloatval_t *stp,
    const lbfgsfloatval_t* xp,
    const lbfgsfloatval_t* gp,
    lbfgsfloatval_t *wp,
    callback_data_t *cd,
    const lbfgs_parameter_t *param
    )
{
    int count = 0;
    lbfgsfloatval_t width, dg;
    lbfgsfloatval_t finit, dginit = 0., dgtest;
    const lbfgsfloatval_t dec = 0.5, inc = 2.1;

    /* Check the input parameters for errors. */
    if (*stp <= 0.) {
        return LBFGSERR_INVALIDPARAMETERS;
    }

    /* Compute the initial gradient in the search direction. */
    vecdot(&dginit, g, s, n);

    /* Make sure that s points to a descent direction. */
    if (0 < dginit) {
        return LBFGSERR_INCREASEGRADIENT;
    }

    /* The initial value of the objective function. */
    finit = *f;
    dgtest = param->ftol * dginit;

    for (;;) {
        veccpy(x, xp, n);
        vecadd(x, s, *stp, n);

        /* Evaluate the function and gradient values. */
        *f = cd->proc_evaluate(cd->instance, x, g, cd->n, *stp);

        ++count;

        if (*f > finit + *stp * dgtest) {
            width = dec;
        } else {
            /* The sufficient decrease condition (Armijo condition). */
            if (param->linesearch == LBFGS_LINESEARCH_BACKTRACKING_ARMIJO) {
                /* Exit with the Armijo condition. */
                return count;
	        }

	        /* Check the Wolfe condition. */
	        vecdot(&dg, g, s, n);
	        if (dg < param->wolfe * dginit) {
    		    width = inc;
	        } else {
		        if(param->linesearch == LBFGS_LINESEARCH_BACKTRACKING_WOLFE) {
		            /* Exit with the regular Wolfe condition. */
		            return count;
		        }

		        /* Check the strong Wolfe condition. */
		        if(dg > -param->wolfe * dginit) {
		            width = dec;
		        } else {
		            /* Exit with the strong Wolfe condition. */
		            return count;
		        }
            }
        }

        if (*stp < param->min_step) {
            /* The step is the minimum value. */
            return LBFGSERR_MINIMUMSTEP;
        }
        if (*stp > param->max_step) {
            /* The step is the maximum value. */
            return LBFGSERR_MAXIMUMSTEP;
        }
        if (param->max_linesearch <= count) {
            /* Maximum number of iteration. */
            return LBFGSERR_MAXIMUMLINESEARCH;
        }

        (*stp) *= width;
    }
}



static int line_search_backtracking_owlqn(
    int n,
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *f,
    lbfgsfloatval_t *g,
    lbfgsfloatval_t *s,
    lbfgsfloatval_t *stp,
    const lbfgsfloatval_t* xp,
    const lbfgsfloatval_t* gp,
    lbfgsfloatval_t *wp,
    callback_data_t *cd,
    const lbfgs_parameter_t *param
    )
{
    int i, count = 0;
    lbfgsfloatval_t width = 0.5, norm = 0.;
    lbfgsfloatval_t finit = *f, dgtest;

    /* Check the input parameters for errors. */
    if (*stp <= 0.) {
        return LBFGSERR_INVALIDPARAMETERS;
    }

    /* Choose the orthant for the new point. */
    for (i = 0;i < n;++i) {
        wp[i] = (xp[i] == 0.) ? -gp[i] : xp[i];
    }

    for (;;) {
        /* Update the current point. */
        veccpy(x, xp, n);
        vecadd(x, s, *stp, n);

        /* The current point is projected onto the orthant. */
        owlqn_project(x, wp, param->orthantwise_start, param->orthantwise_end);

        /* Evaluate the function and gradient values. */
        *f = cd->proc_evaluate(cd->instance, x, g, cd->n, *stp);

        /* Compute the L1 norm of the variables and add it to the object value. */
        norm = owlqn_x1norm(x, param->orthantwise_start, param->orthantwise_end);
        *f += norm * param->orthantwise_c;

        ++count;

        dgtest = 0.;
        for (i = 0;i < n;++i) {
            dgtest += (x[i] - xp[i]) * gp[i];
        }

        if (*f <= finit + param->ftol * dgtest) {
            /* The sufficient decrease condition. */
            return count;
        }

        if (*stp < param->min_step) {
            /* The step is the minimum value. */
            return LBFGSERR_MINIMUMSTEP;
        }
        if (*stp > param->max_step) {
            /* The step is the maximum value. */
            return LBFGSERR_MAXIMUMSTEP;
        }
        if (param->max_linesearch <= count) {
            /* Maximum number of iteration. */
            return LBFGSERR_MAXIMUMLINESEARCH;
        }

        (*stp) *= width;
    }
}



static int line_search_morethuente(
    int n,
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *f,
    lbfgsfloatval_t *g,
    lbfgsfloatval_t *s,
    lbfgsfloatval_t *stp,
    const lbfgsfloatval_t* xp,
    const lbfgsfloatval_t* gp,
    lbfgsfloatval_t *wa,
    callback_data_t *cd,
    const lbfgs_parameter_t *param
    )
{
    int count = 0;
    int brackt, stage1, uinfo = 0;
    lbfgsfloatval_t dg;
    lbfgsfloatval_t stx, fx, dgx;
    lbfgsfloatval_t sty, fy, dgy;
    lbfgsfloatval_t fxm, dgxm, fym, dgym, fm, dgm;
    lbfgsfloatval_t finit, ftest1, dginit, dgtest;
    lbfgsfloatval_t width, prev_width;
    lbfgsfloatval_t stmin, stmax;

    /* Check the input parameters for errors. */
    if (*stp <= 0.) {
        return LBFGSERR_INVALIDPARAMETERS;
    }

    /* Compute the initial gradient in the search direction. */
    vecdot(&dginit, g, s, n);

    /* Make sure that s points to a descent direction. */
    if (0 < dginit) {
        return LBFGSERR_INCREASEGRADIENT;
    }

    /* Initialize local variables. */
    brackt = 0;
    stage1 = 1;
    finit = *f;
    dgtest = param->ftol * dginit;
    width = param->max_step - param->min_step;
    prev_width = 2.0 * width;

    /*
        The variables stx, fx, dgx contain the values of the step,
        function, and directional derivative at the best step.
        The variables sty, fy, dgy contain the value of the step,
        function, and derivative at the other endpoint of
        the interval of uncertainty.
        The variables stp, f, dg contain the values of the step,
        function, and derivative at the current step.
    */
    stx = sty = 0.;
    fx = fy = finit;
    dgx = dgy = dginit;

    for (;;) {
        /*
            Set the minimum and maximum steps to correspond to the
            present interval of uncertainty.
         */
        if (brackt) {
            stmin = min2(stx, sty);
            stmax = max2(stx, sty);
        } else {
            stmin = stx;
            stmax = *stp + 4.0 * (*stp - stx);
        }

        /* Clip the step in the range of [stpmin, stpmax]. */
        if (*stp < param->min_step) *stp = param->min_step;
        if (param->max_step < *stp) *stp = param->max_step;

        /*
            If an unusual termination is to occur then let
            stp be the lowest point obtained so far.
         */
        if ((brackt && ((*stp <= stmin || stmax <= *stp) || param->max_linesearch <= count + 1 || uinfo != 0)) || (brackt && (stmax - stmin <= param->xtol * stmax))) {
            *stp = stx;
        }

        /*
            Compute the current value of x:
                x <- x + (*stp) * s.
         */
        veccpy(x, xp, n);
        vecadd(x, s, *stp, n);

        /* Evaluate the function and gradient values. */
        *f = cd->proc_evaluate(cd->instance, x, g, cd->n, *stp);
        vecdot(&dg, g, s, n);

        ftest1 = finit + *stp * dgtest;
        ++count;

        /* Test for errors and convergence. */
        if (brackt && ((*stp <= stmin || stmax <= *stp) || uinfo != 0)) {
            /* Rounding errors prevent further progress. */
            return LBFGSERR_ROUNDING_ERROR;
        }
        if (*stp == param->max_step && *f <= ftest1 && dg <= dgtest) {
            /* The step is the maximum value. */
            return LBFGSERR_MAXIMUMSTEP;
        }
        if (*stp == param->min_step && (ftest1 < *f || dgtest <= dg)) {
            /* The step is the minimum value. */
            return LBFGSERR_MINIMUMSTEP;
        }
        if (brackt && (stmax - stmin) <= param->xtol * stmax) {
            /* Relative width of the interval of uncertainty is at most xtol. */
            return LBFGSERR_WIDTHTOOSMALL;
        }
        if (param->max_linesearch <= count) {
            /* Maximum number of iteration. */
            return LBFGSERR_MAXIMUMLINESEARCH;
        }
        if (*f <= ftest1 && fabs(dg) <= param->gtol * (-dginit)) {
            /* The sufficient decrease condition and the directional derivative condition hold. */
            return count;
        }

        /*
            In the first stage we seek a step for which the modified
            function has a nonpositive value and nonnegative derivative.
         */
        if (stage1 && *f <= ftest1 && min2(param->ftol, param->gtol) * dginit <= dg) {
            stage1 = 0;
        }

        /*
            A modified function is used to predict the step only if
            we have not obtained a step for which the modified
            function has a nonpositive function value and nonnegative
            derivative, and if a lower function value has been
            obtained but the decrease is not sufficient.
         */
        if (stage1 && ftest1 < *f && *f <= fx) {
            /* Define the modified function and derivative values. */
            fm = *f - *stp * dgtest;
            fxm = fx - stx * dgtest;
            fym = fy - sty * dgtest;
            dgm = dg - dgtest;
            dgxm = dgx - dgtest;
            dgym = dgy - dgtest;

            /*
                Call update_trial_interval() to update the interval of
                uncertainty and to compute the new step.
             */
            uinfo = update_trial_interval(
                &stx, &fxm, &dgxm,
                &sty, &fym, &dgym,
                stp, &fm, &dgm,
                stmin, stmax, &brackt
                );

            /* Reset the function and gradient values for f. */
            fx = fxm + stx * dgtest;
            fy = fym + sty * dgtest;
            dgx = dgxm + dgtest;
            dgy = dgym + dgtest;
        } else {
            /*
                Call update_trial_interval() to update the interval of
                uncertainty and to compute the new step.
             */
            uinfo = update_trial_interval(
                &stx, &fx, &dgx,
                &sty, &fy, &dgy,
                stp, f, &dg,
                stmin, stmax, &brackt
                );
        }

        /*
            Force a sufficient decrease in the interval of uncertainty.
         */
        if (brackt) {
            if (0.66 * prev_width <= fabs(sty - stx)) {
                *stp = stx + 0.5 * (sty - stx);
            }
            prev_width = width;
            width = fabs(sty - stx);
        }
    }

    return LBFGSERR_LOGICERROR;
}



/**
 * Define the local variables for computing minimizers.
 */
#define USES_MINIMIZER \
    lbfgsfloatval_t a, d, gamma, theta, p, q, r, s;

/**
 * Find a minimizer of an interpolated cubic function.
 *  @param  cm      The minimizer of the interpolated cubic.
 *  @param  u       The value of one point, u.
 *  @param  fu      The value of f(u).
 *  @param  du      The value of f'(u).
 *  @param  v       The value of another point, v.
 *  @param  fv      The value of f(v).
 *  @param  du      The value of f'(v).
 */
#define CUBIC_MINIMIZER(cm, u, fu, du, v, fv, dv) \
    d = (v) - (u); \
    theta = ((fu) - (fv)) * 3 / d + (du) + (dv); \
    p = fabs(theta); \
    q = fabs(du); \
    r = fabs(dv); \
    s = max3(p, q, r); \
    /* gamma = s*sqrt((theta/s)**2 - (du/s) * (dv/s)) */ \
    a = theta / s; \
    gamma = s * sqrt(a * a - ((du) / s) * ((dv) / s)); \
    if ((v) < (u)) gamma = -gamma; \
    p = gamma - (du) + theta; \
    q = gamma - (du) + gamma + (dv); \
    r = p / q; \
    (cm) = (u) + r * d;

/**
 * Find a minimizer of an interpolated cubic function.
 *  @param  cm      The minimizer of the interpolated cubic.
 *  @param  u       The value of one point, u.
 *  @param  fu      The value of f(u).
 *  @param  du      The value of f'(u).
 *  @param  v       The value of another point, v.
 *  @param  fv      The value of f(v).
 *  @param  du      The value of f'(v).
 *  @param  xmin    The maximum value.
 *  @param  xmin    The minimum value.
 */
#define CUBIC_MINIMIZER2(cm, u, fu, du, v, fv, dv, xmin, xmax) \
    d = (v) - (u); \
    theta = ((fu) - (fv)) * 3 / d + (du) + (dv); \
    p = fabs(theta); \
    q = fabs(du); \
    r = fabs(dv); \
    s = max3(p, q, r); \
    /* gamma = s*sqrt((theta/s)**2 - (du/s) * (dv/s)) */ \
    a = theta / s; \
    gamma = s * sqrt(max2(0, a * a - ((du) / s) * ((dv) / s))); \
    if ((u) < (v)) gamma = -gamma; \
    p = gamma - (dv) + theta; \
    q = gamma - (dv) + gamma + (du); \
    r = p / q; \
    if (r < 0. && gamma != 0.) { \
        (cm) = (v) - r * d; \
    } else if (a < 0) { \
        (cm) = (xmax); \
    } else { \
        (cm) = (xmin); \
    }

/**
 * Find a minimizer of an interpolated quadratic function.
 *  @param  qm      The minimizer of the interpolated quadratic.
 *  @param  u       The value of one point, u.
 *  @param  fu      The value of f(u).
 *  @param  du      The value of f'(u).
 *  @param  v       The value of another point, v.
 *  @param  fv      The value of f(v).
 */
#define QUARD_MINIMIZER(qm, u, fu, du, v, fv) \
    a = (v) - (u); \
    (qm) = (u) + (du) / (((fu) - (fv)) / a + (du)) / 2 * a;

/**
 * Find a minimizer of an interpolated quadratic function.
 *  @param  qm      The minimizer of the interpolated quadratic.
 *  @param  u       The value of one point, u.
 *  @param  du      The value of f'(u).
 *  @param  v       The value of another point, v.
 *  @param  dv      The value of f'(v).
 */
#define QUARD_MINIMIZER2(qm, u, du, v, dv) \
    a = (u) - (v); \
    (qm) = (v) + (dv) / ((dv) - (du)) * a;

/**
 * Update a safeguarded trial value and interval for line search.
 *
 *  The parameter x represents the step with the least function value.
 *  The parameter t represents the current step. This function assumes
 *  that the derivative at the point of x in the direction of the step.
 *  If the bracket is set to true, the minimizer has been bracketed in
 *  an interval of uncertainty with endpoints between x and y.
 *
 *  @param  x       The pointer to the value of one endpoint.
 *  @param  fx      The pointer to the value of f(x).
 *  @param  dx      The pointer to the value of f'(x).
 *  @param  y       The pointer to the value of another endpoint.
 *  @param  fy      The pointer to the value of f(y).
 *  @param  dy      The pointer to the value of f'(y).
 *  @param  t       The pointer to the value of the trial value, t.
 *  @param  ft      The pointer to the value of f(t).
 *  @param  dt      The pointer to the value of f'(t).
 *  @param  tmin    The minimum value for the trial value, t.
 *  @param  tmax    The maximum value for the trial value, t.
 *  @param  brackt  The pointer to the predicate if the trial value is
 *                  bracketed.
 *  @retval int     Status value. Zero indicates a normal termination.
 *  
 *  @see
 *      Jorge J. More and David J. Thuente. Line search algorithm with
 *      guaranteed sufficient decrease. ACM Transactions on Mathematical
 *      Software (TOMS), Vol 20, No 3, pp. 286-307, 1994.
 */
static int update_trial_interval(
    lbfgsfloatval_t *x,
    lbfgsfloatval_t *fx,
    lbfgsfloatval_t *dx,
    lbfgsfloatval_t *y,
    lbfgsfloatval_t *fy,
    lbfgsfloatval_t *dy,
    lbfgsfloatval_t *t,
    lbfgsfloatval_t *ft,
    lbfgsfloatval_t *dt,
    const lbfgsfloatval_t tmin,
    const lbfgsfloatval_t tmax,
    int *brackt
    )
{
    int bound;
    int dsign = fsigndiff(dt, dx);
    lbfgsfloatval_t mc; /* minimizer of an interpolated cubic. */
    lbfgsfloatval_t mq; /* minimizer of an interpolated quadratic. */
    lbfgsfloatval_t newt;   /* new trial value. */
    USES_MINIMIZER;     /* for CUBIC_MINIMIZER and QUARD_MINIMIZER. */

    /* Check the input parameters for errors. */
    if (*brackt) {
        if (*t <= min2(*x, *y) || max2(*x, *y) <= *t) {
            /* The trival value t is out of the interval. */
            return LBFGSERR_OUTOFINTERVAL;
        }
        if (0. <= *dx * (*t - *x)) {
            /* The function must decrease from x. */
            return LBFGSERR_INCREASEGRADIENT;
        }
        if (tmax < tmin) {
            /* Incorrect tmin and tmax specified. */
            return LBFGSERR_INCORRECT_TMINMAX;
        }
    }

    /*
        Trial value selection.
     */
    if (*fx < *ft) {
        /*
            Case 1: a higher function value.
            The minimum is brackt. If the cubic minimizer is closer
            to x than the quadratic one, the cubic one is taken, else
            the average of the minimizers is taken.
         */
        *brackt = 1;
        bound = 1;
        CUBIC_MINIMIZER(mc, *x, *fx, *dx, *t, *ft, *dt);
        QUARD_MINIMIZER(mq, *x, *fx, *dx, *t, *ft);
        if (fabs(mc - *x) < fabs(mq - *x)) {
            newt = mc;
        } else {
            newt = mc + 0.5 * (mq - mc);
        }
    } else if (dsign) {
        /*
            Case 2: a lower function value and derivatives of
            opposite sign. The minimum is brackt. If the cubic
            minimizer is closer to x than the quadratic (secant) one,
            the cubic one is taken, else the quadratic one is taken.
         */
        *brackt = 1;
        bound = 0;
        CUBIC_MINIMIZER(mc, *x, *fx, *dx, *t, *ft, *dt);
        QUARD_MINIMIZER2(mq, *x, *dx, *t, *dt);
        if (fabs(mc - *t) > fabs(mq - *t)) {
            newt = mc;
        } else {
            newt = mq;
        }
    } else if (fabs(*dt) < fabs(*dx)) {
        /*
            Case 3: a lower function value, derivatives of the
            same sign, and the magnitude of the derivative decreases.
            The cubic minimizer is only used if the cubic tends to
            infinity in the direction of the minimizer or if the minimum
            of the cubic is beyond t. Otherwise the cubic minimizer is
            defined to be either tmin or tmax. The quadratic (secant)
            minimizer is also computed and if the minimum is brackt
            then the the minimizer closest to x is taken, else the one
            farthest away is taken.
         */
        bound = 1;
        CUBIC_MINIMIZER2(mc, *x, *fx, *dx, *t, *ft, *dt, tmin, tmax);
        QUARD_MINIMIZER2(mq, *x, *dx, *t, *dt);
        if (*brackt) {
            if (fabs(*t - mc) < fabs(*t - mq)) {
                newt = mc;
            } else {
                newt = mq;
            }
        } else {
            if (fabs(*t - mc) > fabs(*t - mq)) {
                newt = mc;
            } else {
                newt = mq;
            }
        }
    } else {
        /*
            Case 4: a lower function value, derivatives of the
            same sign, and the magnitude of the derivative does
            not decrease. If the minimum is not brackt, the step
            is either tmin or tmax, else the cubic minimizer is taken.
         */
        bound = 0;
        if (*brackt) {
            CUBIC_MINIMIZER(newt, *t, *ft, *dt, *y, *fy, *dy);
        } else if (*x < *t) {
            newt = tmax;
        } else {
            newt = tmin;
        }
    }

    /*
        Update the interval of uncertainty. This update does not
        depend on the new step or the case analysis above.

        - Case a: if f(x) < f(t),
            x <- x, y <- t.
        - Case b: if f(t) <= f(x) && f'(t)*f'(x) > 0,
            x <- t, y <- y.
        - Case c: if f(t) <= f(x) && f'(t)*f'(x) < 0, 
            x <- t, y <- x.
     */
    if (*fx < *ft) {
        /* Case a */
        *y = *t;
        *fy = *ft;
        *dy = *dt;
    } else {
        /* Case c */
        if (dsign) {
            *y = *x;
            *fy = *fx;
            *dy = *dx;
        }
        /* Cases b and c */
        *x = *t;
        *fx = *ft;
        *dx = *dt;
    }

    /* Clip the new trial value in [tmin, tmax]. */
    if (tmax < newt) newt = tmax;
    if (newt < tmin) newt = tmin;

    /*
        Redefine the new trial value if it is close to the upper bound
        of the interval.
     */
    if (*brackt && bound) {
        mq = *x + 0.66 * (*y - *x);
        if (*x < *y) {
            if (mq < newt) newt = mq;
        } else {
            if (newt < mq) newt = mq;
        }
    }

    /* Return the new trial value. */
    *t = newt;
    return 0;
}





static lbfgsfloatval_t owlqn_x1norm(
    const lbfgsfloatval_t* x,
    const int start,
    const int n
    )
{
    int i;
    lbfgsfloatval_t norm = 0.;

    for (i = start;i < n;++i) {
        norm += fabs(x[i]);
    }

    return norm;
}

static void owlqn_pseudo_gradient(
    lbfgsfloatval_t* pg,
    const lbfgsfloatval_t* x,
    const lbfgsfloatval_t* g,
    const int n,
    const lbfgsfloatval_t c,
    const int start,
    const int end
    )
{
    int i;

    /* Compute the negative of gradients. */
    for (i = 0;i < start;++i) {
        pg[i] = g[i];
    }

    /* Compute the psuedo-gradients. */
    for (i = start;i < end;++i) {
        if (x[i] < 0.) {
            /* Differentiable. */
            pg[i] = g[i] - c;
        } else if (0. < x[i]) {
            /* Differentiable. */
            pg[i] = g[i] + c;
        } else {
            if (g[i] < -c) {
                /* Take the right partial derivative. */
                pg[i] = g[i] + c;
            } else if (c < g[i]) {
                /* Take the left partial derivative. */
                pg[i] = g[i] - c;
            } else {
                pg[i] = 0.;
            }
        }
    }

    for (i = end;i < n;++i) {
        pg[i] = g[i];
    }
}

static void owlqn_project(
    lbfgsfloatval_t* d,
    const lbfgsfloatval_t* sign,
    const int start,
    const int end
    )
{
    int i;

    for (i = start;i < end;++i) {
        if (d[i] * sign[i] <= 0) {
            d[i] = 0;
        }
    }
}
/*

END OF IN-LINE LBFGS CODE Block (Not in Header Files)

*/
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

	static char time_buffer[TIME_SIZE];
	const struct tm *tm;
	size_t len;
	time_t now;

	now = time ( NULL );
	tm = localtime ( &now );

	len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

	printf ( "%s\n", time_buffer );
	fprintf (kp, "%s\n", time_buffer );
	fprintf (lp, "%s\n", time_buffer );

	return;
# undef TIME_SIZE
}
