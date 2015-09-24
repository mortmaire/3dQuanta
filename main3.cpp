#include <iostream>
#include <fftw3.h>
#include <complex>
#include <png++/png.hpp>
#include <fstream>

using namespace std;
using namespace png;

#define complex complex<double>

class parameters{
public:
    int n;
    double L;
    double dx;
    double dt;
    parameters(int n=512,double L=0.1,double dt=0.1):n(n),L(L),dt(dt){
        dx=L/n;
    }
    void refresh(){
        dx=L/n;
    }
    
}p0;

class ffal{
    complex *tab;
    int n;
    fftw_plan p1,p2;
    fftw_complex *in,*out;
    double L;
    double dx;
    parameters p;
    double *kk;
    double *rr;
public:
    ffal(parameters p=p0):p(p){
        L=p.L;
        n=p.n;
        dx=p.dx;
        tab=new complex[n];
        rr=new double[n];

        for(int i=0;i<n;i++)tab[i]=1;
        in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
        
        p1  = fftw_plan_dft_1d(n,in,out, 1,FFTW_ESTIMATE);
        p2  = fftw_plan_dft_1d(n,out,in,-1,FFTW_ESTIMATE);
        kk=new double[n];
        
    }
    ~ffal(){
    delete[] tab;
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(in);
    fftw_free(out);
    delete[] kk;
    delete[] rr;
    }
private:
    double normuj(){
        double sum=0;
        for(int i=0;i<n;i++)sum+=norm(tab[i]);
        sum=sqrt(sum)*dx;
        for(int i=0;i<n;i++)tab[i]/=sum;
        return sum;
    }

    void mult(){

        }
public:
    void evolve(int N=10){
        for(int k=0;k<n;k++){
        kk[k]=pow(k<n/2?2*k*M_PI/L:(k-n)*2*M_PI/L,2)/2;
        }
        for(int i=0;i<n;i++)
            rr[i]=i==0?20:1/((i%n)*dx);
        double nn=2;
        for(int ii=0;ii<N;ii++){
            int t0=time(0);
                    for(int i=0;i<n;i++){
                    in[i][0]=tab[i].real();
                    in[i][1]=tab[i].imag();
                }
            fftw_execute(p1);
            for(int i=0;i<n;i++){
                double e=exp(-kk[i]*p.dt);
                out[i][0]*=e;
                out[i][1]*=e;
            }
        fftw_execute(p2);
        for(int i=0;i<n;i++)tab[i]=complex(in[i][0]/n,in[i][1]/n);
            for(int i=0;i<n;i++){
                tab[i]*=exp((rr[i]-2)*p.dt);
        }

            int t1=time(0);
            nn=normuj();
            cerr<<nn<<endl;
    
        }
     for(int i=0;i<n;i++)
      cout<<i*dx<<" "<<rr[i]<<endl;
     cout<<endl<<endl;
    }
    void show(){
        for(int i=0;i<n;i++)cout<<i*dx<<" "<<abs(tab[i])<<endl;
        cout<<endl<<endl;
    }
    void plot(int a1){
        image< rgb_pixel > image(n,n);
        double max=0;
        for(int i=0;i<n;i++)if(abs(tab[i])>max)max=abs(tab[i]);
int j=n/2;
     for(int i=0;i<n;i++){
                double val=abs(tab[i])/max;
                int r=sqrt(val)*255;
                int g=val*val*val*val*255;
                int b=val*(0.5-val)*16*255;
                if(b<0)b=0;
                image[i/n%n][i%n]=rgb_pixel(r,g,b);
            }
     char s[64];
     sprintf(s,"%02d.png",a1);
     image.write(s);

 
     
 }
    
    };
    
int main(){
    p0.n=2048;
    p0.L=100;
    p0.dt=0.01;
    p0.refresh();
    ffal psi;
//     cout<<p0.dx<<endl;
    psi.evolve(100);
    psi.show();
}