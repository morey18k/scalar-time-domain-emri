#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include <complex.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_odeiv2.h>


typedef struct csts{
    double rho_minus;
    double rho_plus;
    double horizon;
    double scri_plus;
    double xm;    
}Constant;

typedef struct function{
    double * boosta;
    double * boostderiva;
    double * omega2dLa;
    double * omega2dLderiva;
    double * boostb;
    double * boostderivb;
    double * omega2dLb;
    double * omega2dLderivb;
}Function;

//function prototypes
double p(int l, int m,double x);
double complex y(int l, int m, double theta, double phii);
double complex ystar(int l, int m, double theta, double phii);
int odefunc (double x, const double *y, double *f,void *params);
double rst2rsch (double rstar, Constant * constants);
void initdata (int im,int imm,double delrho,double complex *z,double complex *phi,double complex *pi,double *xa);
double vpotential (double rsch,double rschm2m,double xm,int lval);
void mesh (Constant * constants,double width, double delrho,double *r0,int *im, int *imm, int *ismhf,int iwidth,double rsch0,double rstar0, double ** ra, double ** rb, double ** ram2m, double ** rbm2m, double ** xa, double ** xb, double **rhoa, double **rhob, double **erra, double ** errb);
double omega(double rho, Constant * constants);
double omega_prime(double rho, Constant * constants);
double omega_double_prime(double rho, Constant * constants);
double H_bar(double rho, Constant * constants);
double H_bar_prime_rho(double rho, Constant * constants);
double H_bar_prime_x(double rho, Constant * constants);
double H(double rho, Constant * constants);
double H_prime_rho(double rho, Constant * constants);
double H_prime_x(double rho, Constant * constants);
double Lfunc(double rho, Constant * constants);
double Lfunc_prime(double rho, Constant * constants);
double h(double rho, Constant * constants);
void profile_source(double complex * pi, double complex piag, double complex piagp1, int ismhf, int isphf);

int main (void){
    int lmode;
    //would like to keep these variables defined locally in case any change is needed
    //(constants->xm) is black hole mass
    double delrho=0.015;
    double rcour,dt,rsch0=10.0;
    double width=400.0;
    double *ra, * rb,*ram2m,*rbm2m;
    double * xa, *xb;
    double * rhoa, *rhob;
    double * erra, *errb;
    Constant * constants = malloc(sizeof(Constant));
    constants->xm=1.0;
    int iwidth=floor(width/delrho);
    //declares integer values
    int itest,itest1,itest2,itest3,itest4,itest5;
    int im,imm,ismhf,isphf,nperiod;
    //other valuable constants
    double rstar0=rsch0+2*(constants->xm)*log((rsch0-2*(constants->xm))/(2*(constants->xm)));
    constants->horizon=rstar0-iwidth*delrho-0.5*delrho;
    constants->scri_plus=rstar0+iwidth*delrho+0.5*delrho;
    printf("horizon=%f, scri_plus=%f\n", constants->horizon, constants->scri_plus);
    int iwidth2=floor(0.25*width/delrho)+0.5*delrho;
    //declares integer values
    constants->rho_minus=rstar0-iwidth2*delrho;
    constants->rho_plus=rstar0+iwidth2*delrho;
    printf("rho_minus=%f, rho_plus=%f\n",constants->rho_minus, constants->rho_plus);
    double totflux=0.0,totflux1=0.0,totflux2=0.0,totflux3=0.0,totflux4=0.0,totflux5=0.0;
    double xachk,rsch,rschm2m,r0m2m;
    //energy is flux at infinity, hflux is flux at horizon
    double f0,r0,fp0;
    double totforce;
    //sets up mesh and returns r0, im, imm, and ismhf
    mesh(constants,width,delrho,&r0,&im,&imm,&ismhf,iwidth,rsch0,rstar0, &ra, &rb, &ram2m, &rbm2m, &xa, &xb, &rhoa, &rhob, &erra, &errb);
    double freq=sqrt((constants->xm)/pow(r0,3));
    double t0 = (2.0*M_PI)/freq;
    //FILE * meshptra = fopen("mesha.csv", "w");
    //FILE * meshptrb = fopen("meshb.csv", "w");
    //int i;
    //printf("horizon=%f   scri_plus=%f\n", horizon, scri_plus);
    //for (i=0;i<=im;i++){
    //    fprintf(meshptra, "%05d %+.6e %+.6e  %+.6e\n",i, rhoa[i],xa[i],ram2m[i]);
    //    fprintf(meshptrb, "%05d %+.6e %+.6e  %+.6e\n",i, rhob[i],xb[i],rbm2m[i]);
    //}
    //fclose(meshptra);
    //fclose(meshptrb);
    rcour=0.99;
    dt=rcour*delrho;
    nperiod=ceil(t0/dt);
    printf("number of time steps in period=%d\n",nperiod);
    dt=t0/nperiod;
    rcour=dt/delrho;
    
    //calculates zone face on other side of source
    isphf=ismhf+1;
    //sets up test value for measuring energy flux and value of wave
    itest=im;
    itest1=0;
    
    //calculates various information about source radius used in jump conditions
    f0   =(r0-2*(constants->xm))/r0;
    fp0  =(2*(constants->xm))/pow(r0,2);
    r0m2m=r0-2*(constants->xm);
    double vr0  = vpotential(r0,r0-2*(constants->xm),(constants->xm),lmode);
    //vr0=((r0-2*(constants->xm))/(r0))*((lmode*(lmode+1))/pow(r0,2.0)+(2*constants->xm)/pow(r0,3.0));
    printf("pot at source=%+.10e\n",vr0); 
    int n, ntot=70000; 
    
    //scalar charge of source
    double q=1.0;
    //Keplerian frequency of source
    Function * functions=malloc(sizeof(Function));
    
    functions->boosta = malloc(sizeof(double)*(im+1));
    functions->boostderiva = malloc(sizeof(double)*(im+1));
    functions->omega2dLa = malloc(sizeof(double)*(im+1));
    functions->omega2dLderiva = malloc(sizeof(double)*(im+1));
    functions->boostb = malloc(sizeof(double)*(im));
    functions->boostderivb = malloc(sizeof(double)*(im));
    functions->omega2dLb = malloc(sizeof(double)*(im));
    functions->omega2dLderivb = malloc(sizeof(double)*(im));
    int i;
    
    FILE * funcptr = fopen("functions.csv", "w");
    for (i=0;i<=im;i++){
        functions->boosta[i]= H(rhoa[i],constants);
        functions->boostderiva[i]=H_prime_rho(rhoa[i], constants);
        functions->omega2dLa[i] = pow(omega(rhoa[i],constants),2.0)/Lfunc(rhoa[i],constants);
        functions->omega2dLderiva[i] = H_bar_prime_rho(rhoa[i],constants);
        //fprintf(funcptr,"%d, %+.2e, %+.2e, %+.2e, %+.2e, %+.2e, %+.2e, %+.2e\n",
        //i, rhoa[i], functions->boosta[i], functions->boostderiva[i], functions->compressa[i],
        //functions->compressderiva[i], functions->lfunca[i], functions->lfuncderiva[i]);
        if (i==im)
            break;
        functions->boostb[i]= H(rhob[i],constants);
        functions->boostderivb[i]=H_prime_rho(rhob[i],constants);
        functions->omega2dLb[i] = pow(omega(rhob[i],constants),2.0)/Lfunc(rhob[i],constants);
        functions->omega2dLderiva[i] = H_bar_prime_rho(rhob[i],constants);
    }
    printf("%d, %d\n", ismhf, isphf);
    fprintf(funcptr, "a, b, c, d,");
    for (i = 0; i<=im; i+=10){
        fprintf(funcptr,"%e, %e, %e,\n", rhoa[i], omega(rhoa[i],constants), H(rhoa[i],constants)); 
    }
    fclose(funcptr);
//    #pragma omp parallel for schedule(dynamic)   
    for (lmode=1;lmode<=1;lmode++){
        int mmode;
        for (mmode=1;mmode<=lmode;mmode++){
            char fluxfilename[20];
            sprintf(fluxfilename, "flux%02d,%02dmode.csv",lmode,mmode);

            FILE *forceptr=fopen("force.csv","w");
            //FILE *fluxptr=fopen(fluxfilename,"w");
            
            FILE *ryptr= fopen("energy.txt","w");
            time_t start_t, end_t;
            double diff_t;
            
            time(&start_t);
            //printf("l=%d,m=%d,thread=%d\n", lmode,mmode,omp_get_thread_num());
            //frequency by m mode
            double pitempleft,pitempright;
            double complex * z=malloc(sizeof(double complex)*(im+5));
            double complex * pi=malloc(sizeof(double complex)*(im+5));
            double complex * phi=malloc(sizeof(double complex)*(im+5));
            double complex * zinit=malloc(sizeof(double complex)*(im+5));
            double complex * za =malloc(sizeof(double complex)*(im+5));
            double complex * zb =malloc(sizeof(double complex)*(im+5));
            double complex * pia=malloc(sizeof(double complex)*(im+5));
            double complex * pib=malloc(sizeof(double complex)*(im+5));
            double complex * phib=malloc(sizeof(double complex)*(im+5));
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
            double * energy = malloc(sizeof(double)*ntot);
            double * energy1 = malloc(sizeof(double)*ntot);
            double * forcetarray = malloc(sizeof(double)*ntot);
            double * energy2;
            //,energy1[60000],energy2[60000],energy3[60000];
            //double energy4[60000],energy5[60000];
            //calculates q_l_m coefficients
            double complex qlm =-4.0*M_PI*q*ystar(lmode,mmode,M_PI/2.0,0)*(1.0/r0)*sqrt(1-(3.0*(constants->xm))/r0);
            //prints frequency of source to the screen
            //calculates a potential array from the array of r values
            //sets up the initial values of all the arrays
            initdata (im,imm,delrho,z,phi,pi,xa);
            //time loop
            double * va = malloc(sizeof(double)*(im+1));
            double * vb = malloc(sizeof(double)*(im));
            for (i=0;i<=im;i++){
                if (i==0||i==im){
                    va[i]=0;
                }
                else{
                    va[i]=vpotential(ra[i],ram2m[i],(constants->xm),lmode);
                }
                va[i]=vpotential(ra[i],ram2m[i],(constants->xm),lmode);
                if (i==im)
                    break;
                vb[i]=vpotential(rb[i],rbm2m[i], (constants->xm), lmode);
            }
            double complex field, tderiv;
            FILE *rxptr= fopen("bhpert.csv","w");
            FILE *fptr = fopen("sf.csv","w");
            
            for (n=0;n<ntot;n++){
                //calculates time
                t=n*dt;
                
                if (n%10==0){
                    printf("\r%d", n);
                    fflush(stdout);
                }
                //complex functions of time that multiply delta function
                //use logistic function to reduce transient behavior
                //gb   = -qlm*f0*cexp(-I*omega*t)/(1.0+exp(15.0-0.5*t));
                
                gb   = qlm*f0*cexp(-I*omega*t);
                //printf("r0=%e\n",r0);
                //dgb  = -qlm*f0*(0.5*cexp(15.0-0.5*t-I*t*omega)/pow(1.0+exp(15-0.5*t),2.0))-(I*omega*cexp(I*omega*t)/(1.0+exp(15.0-0.5*t)));
                dgb=-I*omega*qlm*f0*cexp(-I*omega*t);
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
                
                zag      =z[ismhf]+jpsi-0.5*jdpsi*delrho;
                zagp1    =z[isphf]-jpsi-0.5*jdpsi*delrho;
                phiag    =phi[ismhf]+jphi-0.5*jdphi*delrho;
                phiagp1  =phi[isphf]-jphi-0.5*jdphi*delrho;
                piag     =pi[ismhf]+jpi-0.5*jdpi*delrho;
                piagp1   =pi[isphf]-jpi-0.5*jdpi*delrho;
                
                laxwen(im,imm, delrho, rcour, dt, z, phi, pi,va, vb,zag,zagp1,phiag,phiagp1, piag,piagp1,isphf,ismhf, za, zb, pia,pib, phib, functions);
           //uses ghost zone variables to move wave forward one step
                fprintf(rxptr,"%d, %e, %e, %e\n",n, t, creal(z[im]), cimag(z[im]));    
                //calculates the energy flux at the test point
                energyv = -(pi[itest]*conj(H_bar(rhoa[itest],constants)*phi[itest]-H(rhoa[itest],constants)*pi[itest]))/(4.0*M_PI);
                energyv1 = pi[itest1]*conj(pi[itest1])/(4.0*M_PI);
                
                energy[n-1]=creal(energyv);
                energy1[n-1]=creal(energyv1);

                pitempleft=pi[ismhf]+0.5*(pi[ismhf]-pi[ismhf-1]);
                pitempright=pi[isphf]-0.5*(pi[isphf+1]-pi[isphf]);
                
                forcet= (1/r0)*0.5*(pi[ismhf]+piagp1)*y(lmode,mmode,0.5*M_PI,omega*t);
                field = (1/r0)*0.5*(z[ismhf]+zagp1)*y(lmode,mmode,0.5*M_PI,omega*t);
                tderiv= (1/r0)*0.5*(pi[ismhf]+piagp1)*y(lmode,mmode,0.5*M_PI,omega*t);
                fprintf(fptr, "%d,    %e,    %e,\n", n, creal(forcet), cimag(forcet));
                if ((n>nperiod)&&(((t0/dt)-floor(t0/dt))<1e-10)){
                    
                    int k;
                    sumflux=0.0;
                    sumflux1=0.0;
                    //sumflux2=0.0;
                    //sumflux3=0.0;
                    //sumflux4=0.0;
                    //sumflux5=0.0;
                    sumforce=0.0;
                    for(k=n-nperiod;k<=n-1;k++){
                        sumflux+=energy[k]*dt;
                        sumflux1+=energy1[k]*dt;
                        sumforce+=forcetarray[k]*dt;
                    }
                    eaflux=sumflux/t0;
                    eaflux1=sumflux1/t0;
                    selft=sumforce/t0;
                }
                //fprintf(fluxptr,"%5d,%.4f,%.12e\n",n,t,eaflux);
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
                /*
                if (icyclep>=2000) {
                    ifilenum++;
                    char filename [20];
                    sprintf(filename, "dn%02d,%02d,%02dmode.csv", ifilenum,lmode,mmode);
                    FILE *dnptr= fopen(filename,"w");
                    for (i=0;i<=im;i++){
                        double pp=creal(pi[i]);
                        double up=cimag(pi[i]);
                        if (fabs(pp)<1e-60)
                            pp=0.0;
                        if (fabs(up)<1e-60)
                            up=0.0;
                        fprintf(dnptr,"%7d,   %10e,   %10e,   %10e,   %10e,\n",i,xa[i],ra[i],pp,up);
                    }
                    fclose(dnptr);
                    icyclep=0;
                }
                */
                //printf("hor=%e\n", eaflux1);
            }

            fclose(rxptr);
            fclose(ryptr);
            
            totforce+=2.0*selft;
            totflux+=2.0*eaflux;
            totflux1+=2.0*eaflux1;
            printf("l=%02d,m=%02d, mode flux= %.12e, %.12e total flux=%.12e,  %.12e, thread=%d\n",lmode, mmode,eaflux,eaflux1,totflux,totflux1,omp_get_thread_num());
            //printf("self force, time component=%.12e\n",selft);
            time(&end_t);
            diff_t = difftime(end_t, start_t);
            //printf ("Zone cycles per second = %e\n",((double) im*ntot)/diff_t);
            fclose(fptr);
            //POINTER SEARCH
            /*
            printf("energy=%p\n",energy);
            printf("energy1=%p\n",energy1);
            printf("forcetarray=%p\n",forcetarray);
            printf("z=%p\n",z);
            printf("phi=%p\n",phi);
            printf("pi=%p\n",pi);
            printf("zinit=%p\n",zinit);
            printf("za=%p\n",za);
            printf("zb=%p\n",zb);
            printf("pia=%p\n",pia);
            printf("pib=%p\n",pib);
            printf("phib=%p\n",phib);
            printf("va=%p\n",va);
            printf("vb=%p\n",vb);
            printf("boosta=%p\n",functions->boosta);
            printf("boostderiva=%p\n",functions->boostderiva);
            printf("omega2dLa=%p\n",functions->omega2dLa);
            printf("omega2dLderiva=%p\n",functions->omega2dLderiva);
            printf("boostb=%p\n",functions->boostb);
            printf("boostderivb=%p\n",functions->boostderivb);
            printf("omega2dLb=%p\n",functions->omega2dLb);
            printf("omega2dLderivb=%p\n",functions->omega2dLderivb);
            printf("functions=%p\n",functions);
            printf("ra=%p\n",ra);
            printf("rb=%p\n",rb);
            printf("ram2m=%p\n",ram2m);
            printf("rbm2m=%p\n",rbm2m);
            printf("xa=%p\n",xa);
            printf("xb=%p\n",xb);
            printf("rhoa=%p\n",rhoa);
            printf("rhob=%p\n",rhob);
            printf("erra=%p\n",erra);
            printf("errb=%p\n",errb);
            printf("constants=%p\n",constants);
            */
            fclose(forceptr);
            free(energy);
            free(energy1);
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
            free(va);
            free(vb);
        }
    }
    free(functions->boosta);
    free(functions->boostderiva);
    free(functions->omega2dLa);
    free(functions->omega2dLderiva);
    free(functions->boostb);
    free(functions->boostderivb);
    free(functions->omega2dLb);
    free(functions->omega2dLderivb);
    free(functions);
    free(ra);
    free(rb);
    free(ram2m);
    free(rbm2m);
    free(xa);
    free(xb);
    free(rhoa);
    free(rhob);
    free(erra);
    free(errb);
    free(constants);
}

//end of routine


//compress function
double omega(double rho, Constant * constants){
    if (rho<(constants->rho_minus)){
        return 1.0-pow(((constants->rho_minus)-rho)/((constants->rho_minus)-(constants->horizon)),4);
    }
    else if (rho>(constants->rho_plus)){
        return 1.0-pow((rho-(constants->rho_plus))/((constants->scri_plus)-(constants->rho_plus)),4);
    }
    else{
        return 1.0;   
    }
}

//derivative of compress function
double omega_prime(double rho, Constant * constants){
    if (rho<(constants->rho_minus)){
        return 4.0*(1.0/((constants->rho_minus)-(constants->horizon)))*pow(((constants->rho_minus)-rho)/((constants->rho_minus)-(constants->horizon)),3);
    }
    else if (rho>(constants->rho_plus)){
        return 4.0*(-1.0/((constants->scri_plus)-(constants->rho_plus)))*pow((rho-(constants->rho_plus))/((constants->scri_plus)-(constants->rho_plus)),3);
    }
    else{
        return 0.0;   
    }
}

double omega_double_prime(double rho, Constant * constants){
    if (rho<(constants->rho_minus)){
        return 12.0*(-1.0/pow(((constants->rho_minus)-(constants->horizon)),2))*pow(((constants->rho_minus)-rho)/((constants->rho_minus)-(constants->horizon)),2);
    }
    else if (rho>(constants->rho_plus)){
        return 12.0*(-1.0/pow(((constants->scri_plus)-(constants->rho_plus)),2))*pow((rho-(constants->rho_plus))/((constants->scri_plus)-(constants->rho_plus)),2);
    }
    else{
        return 0.0;   
    }
}

double H_bar(double rho, Constant * constants){
    return pow(omega(rho, constants),2)/Lfunc(rho, constants);
}

double Lfunc(double rho, Constant * constants){
    return omega(rho, constants)-rho*omega_prime(rho,constants); 
}

double Lfunc_prime(double rho, Constant * constants){
    return -rho*omega_double_prime(rho, constants); 
}

double H_bar_prime_rho(double rho, Constant * constants){
    return ((2.0*omega(rho, constants)*omega_prime(rho,constants))/Lfunc(rho, constants))+((rho*pow(omega(rho, constants),2)*omega_double_prime(rho, constants))/(pow(Lfunc(rho, constants),2)));
}

double H_bar_prime_x(double rho, Constant * constants){
    return H_bar_prime_rho(rho, constants)*(1.0/H_bar(rho, constants));
}


double H(double rho, Constant * constants){
    if (rho<(constants->rho_minus)){
        return H_bar(rho, constants)-1.0;
    }
    else if (rho>(constants->rho_plus)){
        return 1.0-H_bar(rho, constants);
    }
    else{
        return 0.0;
    }
}

double H_prime_rho(double rho, Constant * constants){
    if (rho<(constants->rho_minus)){
        return H_bar_prime_rho(rho, constants);
    }
    else if (rho>(constants->rho_plus)){
        return -1.0*H_bar_prime_rho(rho, constants);
    }
    else{
        return 0.0;
    }
}

double H_prime_x(double rho, Constant * constants){
    if (rho<(constants->rho_minus)){
        return H_bar_prime_x(rho, constants);
    }
    else if (rho>(constants->rho_plus)){
        return -1.0*H_bar_prime_x(rho, constants);
    }
    else{
        return 0.0;
    }
}

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
    return gsl_sf_legendre_sphPlm(l, m, x);   
}
//calculates spherical harmonics from assosciated legendre polynomials
double complex y(int l, int m, double theta, double phii){
    if (m>=0)
        return p(l,m,cos(theta))*cexp(I*m*phii);
    else
        //note symmetry of spherical harmonics
        return pow(-1.0,m)*ystar(l, -m, theta, phii);
}
//conjugate of spherical harmonics
double complex ystar(int l, int m, double theta, double phii){
    //printf("the l mode is %d\n", l);
    return conj(y(l,m,theta,phii));

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
double rst2rsch (double rstar, Constant * constants){
    double rsn, drsdr, rtempm2m, rsch, rschm2m;
    if (rstar>(4.0*(constants->xm))){
        rsch    = rstar;
        rschm2m = rsch-2*(constants->xm);
    }
    else {
        rschm2m = 2*(constants->xm)*exp(rstar/(2.0*(constants->xm))-1.0);
        rsch    = rschm2m+2*(constants->xm);
    }
//for some reason my c compiler won't let me put i in the for loop.
    int i;
    for (i=1;i<=20;i++){
        //test value for rstar
        rsn     = rschm2m+2.0*(constants->xm)*(1.0+log(rschm2m/(2.0*(constants->xm))));
        if (isnan(rsn)){
            puts("Nan error");
            exit(-1);
        }
        drsdr   = (rschm2m+2.0*(constants->xm))/rschm2m;
        //adjusts r value
        rtempm2m= rschm2m + (rstar-rsn)/drsdr;
        rschm2m = rtempm2m;
        rsch    = rschm2m + 2.0*(constants->xm);
    }
    return rschm2m;
    
}

//this mesh creates corresponding r values to every x value in the mesh.
void mesh (Constant * constants,double width, double delrho,double *r0,int *im, int *imm, int *ismhf,int iwidth,double rsch0,double rstar0, double ** ra, double ** rb, double ** ram2m, double ** rbm2m, double ** xa, double ** xb, double **rhoa, double **rhob, double **erra, double ** errb){
    int dim =1;
    //arrays to be used for integration
    double delRstar=0.5*delrho;
    
    //sets up the integration tool utilized by GSL.
    
    int imax,i2,i2n;
    
    *im= 2*iwidth+1;
    imax      = *im*2+1;
    
    double * r= malloc(sizeof(double)*(imax+5));
    double * rstar= malloc(sizeof(double)*(imax+5));
    double * rm2m= malloc(sizeof(double)*(imax+5));
    double * Rstar= malloc(sizeof(double)*(imax+5));

    *ra= malloc(sizeof(double)*(imax+5));
    *rb= malloc(sizeof(double)*(imax+5));
    *ram2m= malloc(sizeof(double)*(imax+5));
    *rbm2m= malloc(sizeof(double)*(imax+5));
    *xa= malloc(sizeof(double)*(imax+5));
    *xb= malloc(sizeof(double)*(imax+5));
    *rhoa= malloc(sizeof(double)*(imax+5));
    *rhob= malloc(sizeof(double)*(imax+5));
    *erra= malloc(sizeof(double)*(imax+5));
    *errb= malloc(sizeof(double)*(imax+5));
    long int i;
    
    double rho0 = (constants->scri_plus),  rhof = (constants->horizon); /* start and end of integration interval */
    //initializes arrays, with twice as many entries as the final a or b arrays
    //this is because so the integration happens once, with steps delrstar instead of delx
    Rstar[0]=rho0;
    
    
    //integration, evolves y from xi to x_i+1
    //qqw
    double rhoi,xi;
    long int j;
    for (i = 1; i < imax-1; i++)
    {
        rhoi = rho0 - i *delRstar;
        //printf("%d, %e\n",i, rhoi);
        xi=rhoi/omega(rhoi,constants);
        if(omega(rhoi,constants)==0){
            printf("inf error at rho=%e, i = %d\n",rhoi, i);
        }
        //stores values in arrays
        Rstar[i]=rhoi;
        rstar[i]=xi;
        if(xi<-1400.0){
            rm2m[i]=0;
        }
        else{
            rm2m[i] =rst2rsch(xi,constants);
        }
        r[i]    =rm2m[i]+2*(constants->xm);
    }
    for (i=0; i<imax; ++i){
        //printf("%e\n", Rstar[i]);
    }
    //finds the nearest zone face to put the source on approximately near the sourx
    //this is in the index according to arrays r and rstar
    i2       =  2*iwidth+1;
    //calculates the index of the zone center to the left of the source in terms of a array
    i2n      = iwidth;
    //calculates the source radius
    *r0      = r[i2];
    Rstar[imax-1]=(constants->horizon);
    //zone face to the left of the source
    *ismhf   = i2n;
    //frees the gsl routine
    
    //max i values in the a and b array index
    *imm = *im-1;
    
    //a arrays (ram2m, ra, xa) correspond to zone centers of the mesh
    //b arrays (rbm2m,rb,xb) correspond to zone faces of the mesh
    for (i = 0; i <= *im; i++){
        (*rhoa)[i]    = Rstar[imax-2*(i)-1];
        (*xa)[i]    = rstar[imax-2*(i)-1];
        (*ra)[i]    = r[imax-2*(i)-1];
        (*ram2m)[i] = rm2m[imax-2*(i)-1];
    }
    for(i=0;i<=*imm;i++){
        (*rhob)[i]    = Rstar[imax-2*(i)-2];
        (*xb)[i]    = rstar[imax-2*(i)-2];
        (*rb)[i]    = r[imax-2*(i)-2];
        (*rbm2m)[i] = rm2m[imax-2*(i)-2];
    }
    
    //sum is used for averaging mesh error
    //checks formulas according to backwards solutions
    //sees if they satisfy tolerances
    double sum =0,xachk=0,xbchk=0;
    int iflag  =0;
    int count=0;
    for (i = 2; i <= *imm; i++)
    {
        if((*ram2m)[i]!=0.0){
            xachk       = (*ram2m)[i]+2.0*(constants->xm)*(1.0+log((*ram2m)[i]/(2.0*(constants->xm))));
            (*erra)[i]     = fabs(xachk-(*xa)[i]);
            if((*erra)[i]>1.0e-10){
                iflag     = 1;
                printf("%e %e, %d\n", (*xa)[i], xachk, i);
            }
            if (isnan((*erra)[i])){
                puts("Nan error in mesh");
                exit(-1);
            }
            sum+=(*erra)[i];
            count++;
        }
        if((*rbm2m)[i]!=0.0){
            xbchk       = (*rbm2m)[i]+2.0*(constants->xm)*(1.0+log((*rbm2m)[i]/(2.0*(constants->xm))));
            (*errb)[i]     = fabs(xbchk-(*xb)[i]);
            if((*errb)[i]>1.0e-10){
                iflag=1;
                printf("%e %e %d\n", (*xb)[i],xbchk, i);
            }
            if (isnan((*errb)[i])){
                puts("Nan error in mesh");
            }
            sum += (*errb)[i];
            count++;
        }
    }
    free(r);
    free(rstar);
    free(rm2m);
    free(Rstar);
    printf("Mesh Constructed,Average Mesh Error=%e\n",(double)sum/count);
    if (iflag==0)
        printf("All Radmesh:  All mesh points constructed satisfy tolerances\n");
    else
        printf("Radmesh:  Mesh has points not satisfying tolerances\n");
    

    //exit(-1);
}


//initializes arrays to zero
//NOTE: in my compiler the 0.0*I is necessary, because
//if the imaginary part is not initialized to zero, it will
//be filled with a random number in the stack
//even if it may not do it all the time, it did it once and
//I don't want to take a chance
void initdata (int im,int imm,double delrho,double complex *z,double complex *phi,double complex *pi,double *xa) {
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
    return (rschm2m/(rschm2m+2.0*xm))*((lval*(lval+1))/pow(rsch,2.0)+(2*xm)/pow(rsch,3.0));
}

void profile_source(double complex * pi, double complex piag, double complex piagp1, int ismhf, int isphf){
    printf("%d    %e     %e\n",ismhf-2,creal(pi[ismhf-2]), cimag(pi[ismhf-2]));
    printf("%d    %e     %e\n",ismhf-1,creal(pi[ismhf-1]), cimag(pi[ismhf-1]));
    printf("%d    %e     %e\n",ismhf,creal(pi[ismhf]), cimag(pi[ismhf]));
    printf("piagp1    %e     %e\n",creal(piagp1), cimag(piagp1));
    printf("piag      %e     %e\n",creal(piag), cimag(piag));
    printf("%d    %e     %e\n",isphf,creal(pi[isphf]), cimag(pi[isphf]));
    printf("%d    %e     %e\n",isphf+1,creal(pi[isphf+1]), cimag(pi[isphf+1]));
    printf("%d    %e     %e\n",isphf+2,creal(pi[isphf+2]), cimag(pi[isphf+2]));
    return;
}


//lax wendroff method, moves wave forward one step and uses ghost zone values
void laxwen(int im,int imm, double delrho,double rcour, double dt, double complex *z, double complex *phi, double complex *pi,double *va, double * vb,double complex zag,double complex zagp1,double complex phiag,double complex phiagp1,double complex piag,double complex piagp1,int isphf,int ismhf, double complex *za, double complex *zb, double complex *pia, double complex *pib, double complex *phib, Function * functions){
    int i;
    
    double complex zcp,picp,phicp,zcm,picm,phicm;
    
    
    zcp   = 0.5*(zag + z[isphf])+0.25*dt*(piag + pi[isphf]);
    zcm   = 0.5*(z[ismhf] + zagp1)+0.25*dt*(pi[ismhf] + piagp1);   
    //calculates fluxes near source
    //this is necessary as the fluxes near source use
    //ghost zone values as opposed to actual ones
    picp  = 0.5*(piag+pi[isphf])+(1.0/(1.0-pow(functions->boostb[ismhf],2.0)))*
    (0.5*rcour*pow(functions->omega2dLb[ismhf],2.0)*(phi[isphf]-phiag)
    +0.25*dt*(functions->omega2dLderivb[ismhf])*(functions->omega2dLb[ismhf])*(phi[isphf]+phiag)
    -rcour*(functions->omega2dLb[ismhf])*(functions->boostb[ismhf])*(pi[isphf]-piag)
    -0.25*dt*(functions->omega2dLb[ismhf])*(functions->boostderivb[ismhf])*(pi[isphf]+piag)-0.25*dt*vb[ismhf]*(z[isphf]+zag)
    );
 //   printf("potential= %+.10e, zag=%+.10e %+.10e\n",vb[ismhf],creal(zag),cimag(zag)); 
    phicp = 0.5*(phiag + phi[isphf])+0.5*rcour*(pi[isphf] - piag);
    picm  = 0.5*(pi[ismhf]+piagp1)+(1.0/(1.0-pow(functions->boostb[ismhf],2.0)))*
    (0.5*rcour*pow(functions->omega2dLb[ismhf],2.0)*(phiagp1-phi[ismhf])
    +0.25*dt*(functions->omega2dLderivb[ismhf])*(functions->omega2dLb[ismhf])*(phiagp1+phi[ismhf])
    -rcour*(functions->omega2dLb[ismhf])*(functions->boostb[ismhf])*(piagp1-pi[ismhf])
    -0.25*dt*(functions->omega2dLb[ismhf])*(functions->boostderivb[ismhf])*(piagp1+pi[ismhf])-0.25*dt*vb[ismhf]*(zagp1+z[ismhf])
    );
    phicm = 0.5*(phi[ismhf] + phiagp1)+0.5*rcour*(piagp1 - pi[ismhf]);
   // printf("phiag=%+.10e%+.10eI, piag=%+.10e%+.10eI\n",creal(phiag),cimag(phiag), creal(piag),cimag(piag));    
    //printf("picp=%+.10e%+.10eI, phicp=,%+.10e%+.10eI, picm=%+.10e%+.10eI, phicm=%+.10e%+.10eI\n", creal(picp),cimag(picp),creal(phicp), cimag(phicp), creal(picm),cimag(phicm), creal(phicm),cimag(phicm)); 
   // printf("boost at source=%.10e\n", functions->boostderivb[ismhf]);
    //NOTE: I have corrected the incorrect array handling
    //I start loops at one and end them at imm (one less than xmax)
    //impose boundary conditions at i=0 and i=im (xmin and xmax)
    
    

    //fluxes
    for (i=0;i<=imm;i++){
        zb[i]   = 0.5*(z[i] + z[i+1])+0.25*dt*(pi[i] + pi[i+1]);
        pib[i]  = 0.5*(pi[i]+pi[i+1])+(1.0/(1.0-pow(functions->boostb[i],2.0)))*
        (0.5*rcour*pow(functions->omega2dLb[i],2.0)*(phi[i+1]-phi[i])
        +0.25*dt*(functions->omega2dLderivb[i])*(functions->omega2dLb[i])*(phi[i+1]+phi[i])
        -rcour*(functions->omega2dLb[i])*(functions->boostb[i])*(pi[i+1]-pi[i])
        -0.25*dt*(functions->omega2dLa[i])*(functions->boostderivb[i])*(pi[i+1]+pi[i])-0.25*dt*vb[i]*(z[i+1]+z[i])
        );
        phib[i] = 0.5*(phi[i] + phi[i+1])+0.5*rcour*(pi[i+1] - pi[i]);
    }
    
    //average half time results to center, making exceptions for the
    //faces near the center
    //move pi and z forward using averaged results
    for (i=1;i<=imm;i++){
        if(i==isphf){
            za[i]=0.5*(zb[i]+zcp);
        }
        else if (i==ismhf){
            za[i]=0.5*(zcm+zb[i-1]);
        }
        else{
            za[i]=0.5*(zb[i]+zb[i-1]);
        }
    }
    
        
    double complex other=0.0+0.0*I;
    for (i=1;i<=imm;i++){
        if (i==isphf){
            z[i]+=0.5*dt*(pib[i]+picp);
            pi[i]  += (1.0/(1.0-pow(functions->boosta[i],2.0)))
            *(rcour*pow(functions->omega2dLa[i],2.0)*(phib[i]-phicp)
            +0.5*dt*(functions->omega2dLderiva[i])*(functions->omega2dLa[i])*(phib[i]+phicp)
            -2.0*rcour*(functions->omega2dLa[i])*(functions->boosta[i])*(pib[i]-picp)
            -0.5*dt*(functions->omega2dLa[i])*(functions->boostderiva[i])*(pib[i]+picp)-dt*va[i]*(za[i])
            );
            phi[i]=phi[i]+rcour*(pib[i]-picp);
        }
        else if (i==ismhf){
            z[i]+=0.5*dt*(picm+pib[i]);
            pi[i]  += (1.0/(1.0-pow(functions->boosta[i],2.0)))*
            (rcour*pow(functions->omega2dLa[i],2.0)*(phicm-phib[i-1])
            +0.5*dt*(functions->omega2dLderiva[i])*(functions->omega2dLa[i])*(phicm+phib[i-1])
            -2.0*rcour*(functions->omega2dLa[i])*(functions->boosta[i])*(picm-pib[i-1])
            -0.5*dt*(functions->omega2dLa[i])*(functions->boostderiva[i])*(picm+pib[i-1])-dt*va[i]*(za[i])
            );
            other+=rcour*(phicm-phib[i-1])-dt*va[i]*za[i];
            phi[i]=phi[i]+rcour*(picm-pib[i-1]);
        }
        else {
            z[i]+=0.5*dt*(pib[i]+pib[i-1]);
            pi[i]  += (1.0/(1.0-pow(functions->boosta[i],2.0)))*
            (rcour*pow(functions->omega2dLa[i],2.0)*(phib[i]-phib[i-1])
            +0.5*dt*(functions->omega2dLderiva[i])*(functions->omega2dLa[i])*(phib[i]+phib[i-1])
            -2.0*rcour*(functions->omega2dLa[i])*(functions->boosta[i])*(pib[i]-pib[i-1])
            -0.5*dt*(functions->omega2dLa[i])*(functions->boostderiva[i])*(pib[i]+pib[i-1])-dt*va[i]*(za[i])
            );
            phi[i]=phi[i]+rcour*(pib[i]-pib[i-1]);
        }
    }
    //printf("%+.14e% +.14eI\n",creal(pi[ismhf]), cimag(pi[ismhf]));

//    printf("%+.15f %+.15fI\n", creal(z[ismhf]),cimag(z[ismhf]));
    //boundary conditions
    
    z[0]      = z[1]-(z[2]-z[1]);
    pi[0]     = pi[1]-(pi[2]-pi[1]);
    z[im]     = z[im-1]+(z[im-1]-z[im-2]);
    pi[im]    = pi[im-1]+(pi[im-1]-pi[im-2]);
    
    //more boundary conditions
    phi[0]      = phi[1]-(phi[2]-phi[1]);
    phi[im]    = phi[im-1]+(phi[im-1]-phi[im-2]);
    
    

}

