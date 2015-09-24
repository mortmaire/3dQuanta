#include <iostream>
#include <fftw3.h>
#include <complex>
#include <png++/png.hpp>
#include <iomanip>

using namespace std;
using namespace png;

#define complex complex<double>

class parameters{
public:
    int n;
    double L;
    double dx;
    double dt;
    parameters(int n=2048,double L=20,double dt=0.1):n(n),L(L),dt(dt){
        dx=L/n;
    }
    void refresh(){
        dx=L/n;
    }
    
}p0;

class ffal{
    complex *tab;
    int n;
    int n2;
    fftw_plan p1,p2;
    fftw_complex *in,*out;
    double L;
    double dx;
    parameters p;
//     double *kk;
    double *rr;
    double *ee;
public:
    ffal(parameters p=p0):p(p){
        p.refresh();
        L=p.L;
        n=p.n;
        n2=n*n;
        dx=p.dx;
        tab=new complex[n2];
        rr=new double[n2];
        for(int i=0;i<n2;i++)tab[i]=1;
        ee=new double[n2];
        in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n2);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n2);
        
        p1  = fftw_plan_dft_2d(n,n,in,out, 1,FFTW_ESTIMATE);
        p2  = fftw_plan_dft_2d(n,n,out,in,-1,FFTW_ESTIMATE);
//         kk=new double[n];
        
    }
    ~ffal(){
    delete[] tab;
    delete[] ee;
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(in);
    fftw_free(out);
//     delete[] kk;
    delete[] rr;
    }
    
        double normuj(){
        double sum=0;
        for(int i=0;i<n2;i++)sum+=norm(tab[i]);
        sum=sqrt(sum)*dx;
        for(int i=0;i<n2;i++)tab[i]/=sum;
        return sum;
    }

public:
    void evolve(int N=10,bool ff=1){
        for(int i=0;i<n2;i++){
            int a=i%n<n/2?i%n:i%n-n;
            int b=i/n<n/2?i/n:i/n-n;
            double kx=pow(a*2*M_PI/L,2)/2;
            double ky=pow(b*2*M_PI/L,2)/2;
            ee[i]=exp(-(kx+ky)*p.dt);
        }
        for(int i=0;i<n2;i++){
            double r0=(((pow(i/n-n/2,2)+pow(i%n-n/2,2))*dx*dx));
            rr[i]=exp(-(r0/2.)*p.dt);
//             if(i%n==0 or i/n==0 or i%n==n-1 or i/n==n-1)rr[i]=0;
        }
        
        
        for(int ii=0;ii<N;ii++){
            int t0=time(0);
                    for(int i=0;i<n2;i++){
                    in[i][0]=tab[i].real();
                    in[i][1]=tab[i].imag();
                }
            fftw_execute(p1);
            for(int i=0;i<n2;i++){
                out[i][0]*=ee[i];
                out[i][1]*=ee[i];
            }
        fftw_execute(p2);
        for(int i=0;i<n2;i++)tab[i]=complex(in[i][0],in[i][1])/complex(n2,0);
            for(int i=0;i<n2;i++){
                tab[i]*=rr[i];
        }

            int t1=time(0);
            cerr<<ii<<"\t"<<normuj()<<endl;
            if(ff)plot(ii);
        }
        
    }
    
      void plot(int ii=0){
     image< rgb_pixel > image(n,n);
     double max=0;
     for(int i=0;i<n2;i++)
             if(abs(tab[i])>max)max=abs(tab[i]);

     for(int i=0;i<n2;i++){
                double val=abs(tab[i])/max;
                int r=sqrt(val)*255;
                int g=val*val*val*val*255;
                int b=val*(0.5-val)*16*255;
                if(b<0)b=0;
                image[i/n][i%n]=rgb_pixel(r,g,b);
            }
     char s[64];
     sprintf(s,"%03d.png",ii);
     image.write(s);
 }
    void copy(complex *&tt){
        normuj();
        complex *temp=tt;
        tt=tab;
        tab=temp;
 }
};


int main(){
    int n=256;
    int n2=n*n;
    p0.n=n;
    p0.dt=0.001;
    p0.L=20;
    
    complex *tt;
    {
    ffal psi;
    psi.evolve(500,0);
    psi.copy(tt);
    }
    ffal psi;
//     psi.evolve

}