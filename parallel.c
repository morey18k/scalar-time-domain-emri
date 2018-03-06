#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_odeiv2.h>

//global variables
const long int ip=150000;
double ra[ip],rb[ip],ram2m[ip],rbm2m[ip];
double xa[ip], xb[ip];
double erra[ip],errb[ip];

//function prototypes
double p(int l, int m,double x);
double complex y(int l, int m, double theta, double phii);
double complex ystar(int l, int m, double theta, double phii);
int odefunc (double x, const double *y, double *f,void *params);
double rst2rsch (double rstar,double rsch,double rschm2m,double xm);
void mesh (double xm,double width,double xmin, double xmax,double delx,double sourx,double *r0,int *im, int *imm, int *ismhf,int iwidth,double rsch0,double rstar0);
void initdata (int im,int imm,double delx,double complex *z,double complex *phi,double complex *pi,double *xa);
double vpotential (double rsch,double rschm2m,double xm,int lval);
void laxwen(int im,int imm, double delx,double rcour, double dt, double complex *z, double complex *phi, double complex *pi,double *va, double * vb,double complex zag,double complex zagp1,double complex phiag,double complex phiagp1,double complex piag,double complex piagp1,int isphf,int ismhf, double complex *za, double complex *zb, double complex *pia, double complex *pib, double complex *phib);

int main (void){
    int lmode;

    //would like to keep these variables defined locally in case any change is needed

    double delx=0.03,rcour,sourx=10.0,xm=1.0,dt,rsch0=10.0;
    double width=1000.0;
    int iwidth=floor(width/delx);
    //declares integer values
    int itest,itest1,itest2,itest3,itest4,itest5;
    int im,imm,ismhf,isphf,nperiod;
    //other valuable constants
    double rstar0=rsch0+2*xm*log((rsch0-2*xm)/(2*xm));
    double xmin=rstar0-iwidth*delx-0.5*delx;
    double xmax=rstar0+iwidth*delx+0.5*delx;
    double totflux=0.0,totflux1=0.0,totflux2=0.0,totflux3=0.0,totflux4=0.0,totflux5=0.0;
    double xachk,rsch,rschm2m,r0m2m;
    //energy is flux at infinity, hflux is flux at horizon
    double f0,r0,fp0,vr0;
    double totforce;
    //sets up mesh and returns r0, im, imm, and ismhf
    mesh(xm,width,xmin,xmax,delx,sourx,&r0,&im,&imm,&ismhf,iwidth,rsch0,rstar0);
    
    double freq=sqrt(xm/pow(r0,3));
    double complex totalforcearray[70000];
    double t0 = (2.0*M_PI)/freq;
    
    rcour=0.99;
    dt=rcour*delx;
    nperiod=ceil(t0/dt);
    printf("number of time steps in period=%d\n",nperiod);
    dt=t0/nperiod;
    rcour=dt/delx;
    
    //calculates zone face on other side of source
    isphf=ismhf+1;
    //sets up test value for measuring energy flux and value of wave
    itest=floor((15.0+width-rstar0)/delx)+1;
    itest1=floor((25.0+width-rstar0)/delx)+1;
    itest2=floor((50.0+width-rstar0)/delx)+1;
    itest3=floor((200.0+width-rstar0)/delx)+1;
    itest4=floor((400.0+width-rstar0)/delx)+1;
    itest5=floor((800.0+width-rstar0)/delx)+1;
    
    //calculates various information about source radius used in jump conditions
    f0   =(r0-2*xm)/r0;
    fp0  =(2*xm)/pow(r0,2);
    r0m2m=r0-2*xm;
    vr0  = vpotential(r0,r0m2m,xm,lmode);
    
    int n, ntot=30000; 
    for (n=0;n<ntot;n++){
        totalforcearray[n]=0.0;
    }
    
    //scalar charge of source
    double q=1.0;
    //Keplerian frequency of source
//    #pragma omp parallel for schedule(dynamic) 
    for (lmode=1;lmode<=1;lmode++){
        int mmode;
        for (mmode=1;mmode<=lmode;mmode++){
            char fluxfilename[20];
            sprintf(fluxfilename, "flux%02d,%02dmode.csv",lmode,mmode);

            FILE *forceptr=fopen("force.csv","w");
            FILE *fluxptr=fopen(fluxfilename,"w");
            
            FILE *rxptr= fopen("bhpert.csv","w");
            FILE *ryptr= fopen("energy.txt","w");
            time_t start_t, end_t;
            double diff_t;
            
            time(&start_t);
            //printf("l=%d,m=%d,thread=%d\n", lmode,mmode,omp_get_thread_num());
            //frequency by m mode
            double pitempleft,pitempright;
            double complex * z=malloc(sizeof(double complex)*ip);
            double complex * pi=malloc(sizeof(double complex)*ip);
            double complex * phi=malloc(sizeof(double complex)*ip);
            double complex * zinit=malloc(sizeof(double complex)*ip);
            double complex * za =malloc(sizeof(double complex)*ip);
            double complex * zb =malloc(sizeof(double complex)*ip);
            double complex * pia=malloc(sizeof(double complex)*ip);
            double complex * pib=malloc(sizeof(double complex)*ip);
            double complex * phib=malloc(sizeof(double complex)*ip);
            double t;
            double omega=freq*mmode;
    
            double zonecps;
            double complex gb,dgb,fb,dfb,ddfb;
            double complex zag,zagp1,phiag,phiagp1,piag,piagp1;
            double complex jpsi,jphi,jpi,jdpsi,jdphi,jdpi;
            double sumflux,sumflux1,sumflux2,sumflux3,sumflux4,sumflux5;
            double eaflux,eaflux1,eaflux2,eaflux3,eaflux4,eaflux5;
            double complex forcet; 
            double sumforce;
            double selft;
            int i,n=1,icycle=0,icyclep=0,ifilenum=0;
            double complex energyv,energyv1,energyv2,energyv3,energyv4,energyv5;
            double * energy = malloc(sizeof(double)*60000);
            double * forcetarray = malloc(sizeof(double)*60000);
            //,energy1[60000],energy2[60000],energy3[60000];
            //double energy4[60000],energy5[60000];
            //calculates q_l_m coefficients
            double complex qlm =4.0*M_PI*q*ystar(lmode,mmode,M_PI/2.0,0)*(1.0/r0)*sqrt(1-(3.0*xm)/r0);
            //printf("qlm for l=%02d, m=%02d is %.5f+%.5fI\n",lmode, mmode, creal(qlm),cimag(qlm));        
            //prints frequency of source to the screen
            //calculates a potential array from the array of r values
            FILE *vptr=fopen("potential.csv","w");
            double * va = malloc(sizeof(double)*im);
            double * vb= malloc(sizeof(double)*imm);
            for (i=0;i<=im;i++){
                va[i]=vpotential(ra[i],ram2m[i],xm,lmode);
                if (i==im)
                    break;
                vb[i]=vpotential(rb[i],rbm2m[i],xm,lmode);
                fprintf(vptr,"%d,  %f,  %f,  %f,\n",i,xa[i],xb[i],va[i]);
            }
            fclose(vptr);
            //sets up the initial values of all the arrays
            initdata (im,imm,delx,z,phi,pi,xa);
            //time loop
            for (n=1;n<=ntot;n++){
                if(n%100==0){
                    printf("%d\n",n);
                }
                //calculates time
                t=n*dt;
                
                
                //complex functions of time that multiply delta function
                //use logistic function to reduce transient behavior
                gb   = -qlm*f0*cexp(-I*omega*t)/(1.0+exp(15.0-0.5*t));
                //gb   = -qlm*f0*cexp(-I*omega*t);
                dgb  = -qlm*f0*(0.5*cexp(15.0-0.5*t-I*t*omega)/pow(1.0+exp(15-0.5*t),2.0))-(I*omega*cexp(I*omega*t)/(1.0+exp(15.0-0.5*t)));
                //dgb=I*omega*qlm*f0*cexp(-I*omega*t);
                fb   = 0.0+0.0*I;
                dfb  = 0.0+0.0*I;
                ddfb = 0.0+0.0*I;
                
                //jump conditions
                
                jpsi   = (1.0/pow(f0,2.0))*fb;
                jphi   = fp0*jpsi+((1.0/f0)*gb);
                jpi   = (1.0/pow(f0,2))*dfb;
                jdpsi = jphi;
                jdphi  = (1/pow(f0,2.0))*ddfb+vr0*jpsi;
                jdpi   = fp0*jpi+(1.0/f0)*dgb;
                
                //ghost zone variables
                
                zag=z[ismhf]+jpsi-0.5*jdpsi*delx;
                zagp1    =z[isphf]-jpsi-0.5*jdpsi*delx;
                phiag    =phi[ismhf]+jphi-0.5*jdphi*delx;
                phiagp1  =phi[isphf]-jphi-0.5*jdphi*delx;
                piag     =pi[ismhf]+jpi-0.5*jdpi*delx;
                piagp1   =pi[isphf]-jpi-0.5*jdpi*delx;
                
                //uses ghost zone variables to move wave forward one step
                laxwen (im,imm,delx,rcour,dt,z,phi,pi,va,vb,zag,zagp1,phiag,phiagp1,piag,piagp1,isphf,ismhf,za,zb,pia,pib,phib);
                
                //calculates the energy flux at the test point
                energyv = -(pi[itest]*conj(phi[itest]))/(4.0*M_PI);
                //energyv1 = -(pi[itest1]*conj(phi[itest1]))/(4.0*M_PI);
                //energyv2 = -(pi[itest2]*conj(phi[itest2]))/(4.0*M_PI);
                //energyv3 = -(pi[itest3]*conj(phi[itest3]))/(4.0*M_PI);
                //energyv4 = -(pi[itest4]*conj(phi[itest4]))/(4.0*M_PI);
                //energyv5 = -(pi[itest5]*conj(phi[itest5]))/(4.0*M_PI);
                
                energy[n-1]=creal(energyv);
                //energy1[n-1]=creal(energyv1);
                //energy2[n-1]=creal(energyv2);
                //energy3[n-1]=creal(energyv3);
                //energy4[n-1]=creal(energyv4);
                //energy5[n-1]=creal(energyv5);

                pitempleft=pi[ismhf]+0.5*(pi[ismhf]-pi[ismhf-1]);
                pitempright=pi[isphf]-0.5*(pi[isphf+1]-pi[isphf]);
                
                forcet= (1/r0)*0.5*(pi[ismhf]+pi[isphf])*y(lmode,mmode,0.5*M_PI,omega*t);
                forcetarray[n-1]=creal(forcet);
                totalforcearray[n-1]+=forcet;
                if ((n>nperiod)&&(((t0/dt)-floor(t0/dt))<1e-10)){
                    
                    int k;
                    sumflux=0.0;
                    //sumflux1=0.0;
                    //sumflux2=0.0;
                    //sumflux3=0.0;
                    //sumflux4=0.0;
                    //sumflux5=0.0;
                    sumforce=0.0;
                    for(k=n-nperiod;k<=n-1;k++){
                        sumflux+=energy[k]*dt;
                        //sumflux1+=energy1[k]*dt;
                        //sumflux2+=energy2[k]*dt;
                        //sumflux3+=energy3[k]*dt;
                        //sumflux4+=energy4[k]*dt;
                        //sumflux5+=energy4[k]*dt;
                        sumforce+=forcetarray[k]*dt;
                    }
                    eaflux=sumflux/t0;
                    //eaflux1=sumflux1/t0;
                    //eaflux2=sumflux2/t0;
                    //eaflux3=sumflux3/t0;
                    //eaflux4=sumflux4/t0;
                    //eaflux5=sumflux5/t0;
                    selft=sumforce/t0;
                }
                fprintf(fluxptr,"%5d,%.4f,%.12e\n",n,t,eaflux);
                //prints values of wave and energy flux at test point to file
                //fprintf(ryptr,"%7d,    %+e,    %+e,    %+e,\n",n,time,creal(z[itest1]),cimag(z[itest1]));
                //routine for displaying progress on screen
                icycle++;
                //if (icycle>=100) {
                //    printf("cycle= %d\n",n);
                //    icycle=0;
                //}

                //routine that prints value of wave at constant time through various x positions
                //does this every 2000 iterations or every 59.4 M
                //file names titled dn##.txt, where ## is the file number
                icyclep++;
                if (icyclep>=2000) {
                    ifilenum++;
                    char filename [20];
                    sprintf(filename, "dn%02d,%02d,%02dmode.csv", ifilenum,lmode,mmode);
                    FILE *dnptr= fopen(filename,"w");
                    for (i=2;i<=imm-1;i++){
                        double pp=creal(z[i]);
                        double up=cimag(z[i]);
                        if (fabs(pp)<1e-60||fabs(pp)>4.0)
                            pp=0.0;
                        if (fabs(up)<1e-60||fabs(pp)>4.0)
                            up=0.0;
                        fprintf(dnptr,"%7d,   %10e,   %10e,   %10e,   %10e,\n",i,xa[i],ra[i],pp,up);
                    }
                    fclose(dnptr);
                    icyclep=0;
                }
            }
            fclose(rxptr);
            fclose(ryptr);
            
            totforce+=2.0*selft;
            totflux+=2.0*eaflux;
            //totflux1+=2.0*eaflux1;
            //totflux2+=2.0*eaflux2;
            //totflux3+=2.0*eaflux3;
            //totflux4+=2.0*eaflux4;
            //totflux5+=2.0*eaflux5;
            //fprintf(fluxptr,"%d     %d     %+.12e     %+.12e      %+.12e      %+.12e      %+.12e     %+.12e \n",lmode,mmode,eaflux,eaflux1,eaflux2,eaflux3,eaflux4,eaflux5);
            printf("l=%02d,m=%02d, partial total flux=%.12e, thread=%d\n",lmode, mmode,totflux,omp_get_thread_num());
            //printf("self force, time component=%.12e\n",selft);
            time(&end_t);
            diff_t = difftime(end_t, start_t);
            //printf ("Zone cycles per second = %e\n",((double) im*ntot)/diff_t);
            fclose(forceptr);
            free(va);
            free(vb);
            free(energy);
            free(forcetarray);
            free(z);
            free(phi);
            free(pi);
            free(zinit);
            free(za);
            free(zb);
            free(pia);
            free(pib);
            free(phib);
            fclose(fluxptr);	
        }
    }
    FILE *finalptr=fopen("finalflux.csv","w");
    for (n=0;n<ntot;n++){
        fprintf(finalptr,"%5d,%+.4f,%+.12e,%+.12e\n",n,n*dt,creal(totalforcearray[n]),cimag(totalforcearray[n]));
    }
    fclose(finalptr);
    printf("12.5M point, flux=%.12e\n",totflux);
    printf("25.0M point, flux=%.12e\n",totflux1);
    printf("50.0M point, flux=%.12e\n",totflux2);
    printf("200.M point, flux=%.12e\n",totflux3);
    printf("400.M point, flux=%.12e\n",totflux4);
    printf("800.M point, flux=%.12e\n",totflux5);
}

//end of routine

//calculates assosciated legendre polynomials (all positive to get rid
//of c-s phase factor
double p(int l, int m,double x){
    if(m<0){
        puts("m is negative, cannot evaluate.");
        exit(-1);
    }
    if (l<m){
        puts("l is less than m, cannot evaluate.");
        printf("l=%d, m=%d\n",l,m);
        exit(-1);
    }
    if (x<-1.0){
        puts("Cannot input a value less than -1");   
        exit(-1);
    }
    if (x>1.0){
        puts("Cannot input a value greater than 1");
        exit(-1);
    }
    return pow(-1,m)*gsl_sf_legendre_sphPlm(l, m, x);   
}
//calculates spherical harmonics from assosciated legendre polynomials
double complex y(int l, int m, double theta, double phii){
    if (m>=0)
        return p(l,m,cos(theta))*cexp(I*m*phii);
    else
        //note symmetry of spherical harmonics
        return p(l,-m,cos(theta))*cexp(I*m*phii);
}
//conjugate of spherical harmonics
double complex ystar(int l, int m, double theta, double phii){
    //printf("the l mode is %d\n", l);
    return pow(-1,m)*y(l,-m,theta,phii);

}
//differential equation governing r as a function of rstar
//NOTE: this function is what the GSL routines use to solve the
//differential equation. The routines have functions as their parameters,
//this function being one of them. The parameters to this function must be
//very specific in order to be compatible with the GSL scientific library,
//and the compiler got angry when I tried to change the parameters.
int odefunc (double x, const double *y, double *f,void *params)
{
    *f = *y/(*y+2);
    return GSL_SUCCESS;
}
//converts xmax from rstar to rsch, using a guess and check method
//NOTE: I do not need to pass by reference here. I only need to
//pass the value of xmax to the routine and output its corresponding
//rm2m value for use as the initial value in the integration interval.
//I do not need to update the value when I call the function
double rst2rsch (double rstar,double rsch,double rschm2m,double xm){
    double rsn, drsdr, rtempm2m;
    if (rstar>(4.0*xm)){
        rsch    = rstar;
        rschm2m = rsch-2*xm;
    }
    else {
        rschm2m = 2*xm*exp(rstar/(2.0*xm)-1.0);
        rsch    = rschm2m+2*xm;
    }
//for some reason my c compiler won't let me put i in the for loop.
    int i;
    for (i=1;i<=20;i++){
        //test value for rstar
        rsn     = rschm2m+2.0*xm*(1.0+log(rschm2m/(2.0*xm)));
        drsdr   = (rschm2m+2.0*xm)/rschm2m;
        //adjusts r value
        rtempm2m= rschm2m + (rstar-rsn)/drsdr;
        rschm2m = rtempm2m;
        rsch    = rschm2m + 2.0*xm;
    }
    return rschm2m;
    
}

//this mesh creates corresponding r values to every x value in the mesh.
void mesh (double xm,double width,double xmin, double xmax,double delx,double sourx,double *r0,int *im, int *imm, int *ismhf,int iwidth,double rsch0,double rstar0)
{
    //printf("%d\n",iwidth);
    int dim =1;
    double rschmax=xmax,rschm2m,rschmaxm2m;
    //arrays to be used for integration
    double r[ip],rstar[ip],rm2m[ip];
    double delrstar=0.5*delx;
    
    //sets up the integration tool utilized by GSL.
    gsl_odeiv2_system sys = {odefunc, NULL, dim, NULL};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, -1e-10, 1e-10, 1e-10);
    
    int i,imax,i2,i2n;
    double x0 = xmax,  xf = xmin; /* start and end of integration interval */
    double x = x0;
    if (xm<=1e-10) //if flat spacetime
        rschmaxm2m=xmax;
    else
        rschmaxm2m = rst2rsch(xmax,rschmax,rschm2m,xm);
    rschmax= rschmaxm2m+2*xm;
    double y[1]= {rschmaxm2m};  /* initial value */
    //initializes arrays, with twice as many entries as the final a or b arrays
    //this is because so the integration happens once, with steps delrstar instead of delx
    rstar[0]  = xmax;
    r[0]      = rschmax;
    rm2m[0]   = rschmaxm2m;
    *im= 2*iwidth+1;
    imax      = *im*2+1;
    //integration, evolves y from xi to x_i+1
    for (i = 1; i <= imax-1; i++)
    {
        double xi = x0 - i *delrstar;
        //weird gsl syntax
        int status = gsl_odeiv2_driver_apply (d, &x, xi, y);
        
        if (status != GSL_SUCCESS)
        {
            printf ("error, return value=%d\n", status);
            break;
        }
        /*printf ("%.8e %.8e\n", x, y[0]);*/
        //stores values in arrays
        rstar[i]=x;
        rm2m[i] =y[0];
        r[i]    =rm2m[i]+2*xm;
    }
    //finds the nearest zone face to put the source on approximately near the sourx
    //this is in the index according to arrays r and rstar
    i2       =  2*iwidth+1;
    //calculates the index of the zone center to the left of the source in terms of a array
    i2n      = iwidth;
    //calculates the source radius
    *r0      = r[i2];
    //zone face to the left of the source
    *ismhf   = i2n;
    //frees the gsl routine
    gsl_odeiv2_driver_free (d);

    
    //max i values in the a and b array index
    *imm = *im-1;
    //a arrays (ram2m, ra, xa) correspond to zone centers of the mesh
    //b arrays (rbm2m,rb,xb) correspond to zone faces of the mesh
    for (i = 0; i <= *im; i++)
    {
        xa[i]    = rstar[imax-2*(i)-1];
        ra[i]    = r[imax-2*(i)-1];
        ram2m[i] = rm2m[imax-2*(i)-1];
        if (i==*im)
            break;
        xb[i]    = rstar[imax-2*(i)-2];
        rb[i]    = r[imax-2*(i)-2];
        rbm2m[i] = rm2m[imax-2*(i)-2];
    }
    //sum is used for averaging mesh error
    //checks formulas according to backwards solutions
    //sees if they satisfy tolerances
    double sum =0,xachk=0,xbchk=0;
    int iflag  =0;
    for (i = 1; i <= *imm; i++)
    {
        xachk       = ram2m[i]+2.0*xm*(1.0+log(ram2m[i]/(2.0*xm)));
        xbchk       = rbm2m[i]+2.0*xm*(1.0+log(rbm2m[i]/(2.0*xm)));
        erra[i]     = fabs(xachk-xa[i]);
        errb[i]     = fabs(xbchk-xb[i]);
        if(erra[i]>5.0e-10||errb[i]>5.0e-10)
            iflag     = 1;
        sum += erra[i] + errb[i];
    }
    printf("Mesh Constructed,Average Mesh Error=%e\n",(sum/imax));
    if (iflag==0)
        printf("All Radmesh:  All mesh points constructed satisfy tolerances\n");
    else
        printf("Radmesh:  Mesh has points not satisfying tolerances\n");
}


//initializes arrays to zero
//NOTE: in my compiler the 0.0*I is necessary, because
//if the imaginary part is not initialized to zero, it will
//be filled with a random number in the stack
//even if it may not do it all the time, it did it once and
//I don't want to take a chance
void initdata (int im,int imm,double delx,double complex *z,double complex *phi,double complex *pi,double *xa) {
    int i;
    for (i=0;i<= im;i++){
        z[i]=0.0+0.0*I;
    }
    for (i=0;i<= im;i++){
        phi[i]=0.0+0.0*I;
    }
    for (i=0;i<= im;i++){
        pi[i]=0.0+0.0*I;
    }
}


//calculates the scalar potential as a fucntion of the arguments rsch, rschm2m, and l
//useful in creating the vpot array later
double vpotential (double rsch,double rschm2m,double xm,int lval){
    double spin= 2*xm;
    double ang = lval*(lval+1);
    return (rschm2m/(rschm2m+2.0*xm))*(ang/pow(rsch,2.0)+spin/pow(rsch,3.0));
}


//lax wendroff method, moves wave forward one step and uses ghost zone values
void laxwen(int im,int imm, double delx,double rcour, double dt, double complex *z, double complex *phi, double complex *pi,double *va, double * vb,double complex zag,double complex zagp1,double complex phiag,double complex phiagp1,double complex piag,double complex piagp1,int isphf,int ismhf, double complex *za, double complex *zb, double complex *pia, double complex *pib, double complex *phib){
    int i;
    
    double complex zcp,picp,phicp,zcm,picm,phicm;
    
    //calculates fluxes near source
    //this is necessary as the fluxes near source use
    //ghost zone values as opposed to actual ones
    zcp   = 0.5*(zag + z[isphf])+0.25*dt*(piag + pi[isphf]);
    picp  = 0.5*(piag+pi[isphf])+0.5*rcour*(phi[isphf] - phiag)-0.25*dt*vb[ismhf]*(zag+z[isphf]);
    phicp = 0.5*(phiag + phi[isphf])+0.5*rcour*(pi[isphf] - piag);
    zcm   = 0.5*(z[ismhf] + zagp1)+0.25*dt*(pi[ismhf] + piagp1);
    picm  = 0.5*(pi[ismhf]+piagp1)+0.5*rcour*(phiagp1 - phi[ismhf])-0.25*dt*vb[ismhf]*(z[ismhf] + zagp1);
    phicm = 0.5*(phi[ismhf] + phiagp1)+0.5*rcour*(piagp1 - pi[ismhf]);
    
    //NOTE: I have corrected the incorrect array handling
    //I start loops at one and end them at imm (one less than xmax)
    //impose boundary conditions at i=0 and i=im (xmin and xmax)
    
    
    //fluxes
    for (i=1;i<=imm;i++){
        zb[i]   = 0.5*(z[i] + z[i+1])+0.25*dt*(pi[i] + pi[i+1]);
        pib[i]  = 0.5*(pi[i]+pi[i+1])+0.5*rcour*(phi[i+1] - phi[i])-0.25*dt*vb[i]*(z[i] + z[i+1]);
        phib[i] = 0.5*(phi[i] + phi[i+1])+0.5*rcour*(pi[i+1] - pi[i]);
    }
    
    //average half time results to center, making exceptions for the
    //faces near the center
    for (i=1;i<=imm;i++){
        if (i==isphf){
            pia[i]=0.5*(pib[i]+picp);
            za[i]=0.5*(zb[i]+zcp);
        }
        else if (i==ismhf) {
            pia[i]=0.5*(picm+pib[i-1]);
            za[i]=0.5*(zcm+zb[i-1]);
        }
        else{
            pia[i]=0.5*(pib[i]+pib[i-1]);
            za[i]=0.5*(zb[i]+zb[i-1]);
        }
    }
    //move pi and z forward using averaged results
    for (i=1;i<=imm;i++){
        z[i]+=dt*pia[i];
        if (i==isphf){
            pi[i]+=rcour*(phib[i]-phicp)-dt*va[i]*za[i];
        }
        else if (i==ismhf){
            pi[i]+=rcour*(phicm-phib[i-1])-dt*va[i]*za[i];
        }
        else {
            pi[i]+=rcour*(phib[i]-phib[i-1])-dt*va[i]*za[i];
        }
    }
    //boundary conditions
    z[0]      = 0.0+0.0*I;
    pi[0]     = 0.0+0.0*I;
    z[im]     = 0.0+0.0*I;
    pi[im]    = 0.0+0.0*I;
    
    double rate =0.1;
    double complex phibyz;
    //move phi forward
    for (i=1;i<=imm;i++){
        phibyz=(z[i+1]-z[i-1])/(2*delx);
        if (i==isphf) {
            phi[i]=phi[i]+rcour*(pib[i]-picp);
        }
        else if (i==ismhf) {
            phi[i]=phi[i]+rcour*(picm-pib[i-1]);
        }
        else {
            phi[i]=phi[i]+rcour*(pib[i]-pib[i-1]);
        }
    }
    //more boundary conditions
    phi[0]=0.0+0.0*I;
    phi[im]=0.0+0.0*I;

}

    

