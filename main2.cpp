#include <iostream>
#include <fftw3.h>
#include <complex>
#include <png++/png.hpp>


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
    complex *ttab;
    int n;
    int n2;
    int n3;
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
        n2=n*n;
        n3=n2*n;
        dx=p.dx;
        tab=new complex[n3];
        ttab=new complex[n3];
        rr=new double[n3];

        for(int i=0;i<n3;i++)tab[i]=1;
        in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n3);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n3);
        
        p1  = fftw_plan_dft_3d(n,n,n,in,out, 1,FFTW_ESTIMATE);
        p2  = fftw_plan_dft_3d(n,n,n,out,in,-1,FFTW_ESTIMATE);
        kk=new double[n];
        
    }
    ~ffal(){
    delete[] tab;
    delete[] ttab;
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
        for(int i=0;i<n3;i++)sum+=norm(tab[i]);
//         cerr<<endl<<sum<<endl;
        sum=sqrt(sum*dx)*dx;
        for(int i=0;i<n3;i++)tab[i]/=sum;
        return sum;
    }

    void mult(){

        }
public:
    void evolve(int N=10){
        for(int k=0;k<n;k++){
        kk[k]=pow(k<n/2?2*k*M_PI/L:(k-n)*2*M_PI/L,2)/2;
        }
        for(int i=0;i<n3;i++)
            rr[i]=(i==n3/2?10:1/(sqrt(pow(i%n-n/2,2)+pow(i/n%n-n/2,2)+pow(i/n/n-n/2,2))*dx));
        double nn=2;
        for(int ii=0;ii<N;ii++){
            int t0=time(0);
                    for(int i=0;i<n3;i++){
                    in[i][0]=tab[i].real();
                    in[i][1]=tab[i].imag();
                }
            fftw_execute(p1);
            for(int i=0;i<n3;i++){
                out[i][0]*=exp(-(kk[i/n2]+kk[i/n%n]+kk[i%n])*p.dt);
                out[i][1]*=exp(-(kk[i/n2]+kk[i/n%n]+kk[i%n])*p.dt);
            }
        fftw_execute(p2);
        for(int i=0;i<n3;i++)tab[i]=complex(in[i][0]/n3,in[i][1]/n3);
            for(int i=0;i<n3;i++){
                tab[i]*=exp((rr[i])*p.dt);
                if(i%n2==0)cerr<<"\r"<<ii<<"\t"<<i;
        }

            int t1=time(0);
            nn=normuj();
            cerr<<"\t"<<nn<<" "<<n3/(double)(t1-t0)/1e6<<" MCycles"<<endl;
            plot(ii);
    
    }
    }
    void show(){
     
    }
    void plot(int a1){
        image< rgb_pixel > image(n,n);
        double max=0;
        for(int i=n3/2;i<n3/2+n2;i++)if(abs(tab[i])>max)max=abs(tab[i]);
int j=n/2;
     for(int i=0;i<n2;i++){
                double val=abs(tab[j*n2+i])/max;
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
    p0.n=256;
    p0.L=4;
    p0.dt=0.01;
    p0.refresh();
    ffal psi;
    cout<<p0.dx<<endl;
    psi.evolve(100);
//     psi.show();
}