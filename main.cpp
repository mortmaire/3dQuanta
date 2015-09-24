#include <iostream>
#include <fftw3.h>
#include <complex>

using namespace std;

#define complex complex<double>

//chcę użyć równania schrodingera dla at wodoru, -\Delta^2\psi-1\r\psi=E\psi, prawd. kombinacja fft i urojonych czasów

class parameters{
public:
    int n;
    double L;
    double dx;
    double dt;
public:
    parameters(int n=1024,double L=0.1,double dt=0.01):n(n),L(L),dt(dt){
     dx=L/n;   
    }
    void refresh(){
     dx=L/n;   
    }
}p0;


class ffal{
    complex **tab;
    complex *ttab;
    int n;
    int n2;
    fftw_plan p1,p2;
    fftw_complex *in,*out;
    double L;
    double dx;
    parameters p;
public:
    ffal(parameters p=p0):p(p){
        n=p.n;
        n2=n*n;
        tab=new complex*[n];
        ttab=new complex[n2];
        for(int i=0;i<n;i++)tab[i]=new complex[n];
//         for(int i=0;i<n;i++)
//             for(int j=0;j<n;j++)tab[i][j]=new complex[n];
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)tab[i][j]=1;
        in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
        p1  = fftw_plan_dft_1d(n,in,out, 1,FFTW_ESTIMATE);
        p2  = fftw_plan_dft_1d(n,out,in,-1,FFTW_ESTIMATE);
        
        dx=p.dx;
    }
private:
    void normuj(){
        double sum=0;
     for(int i=0;i<n2;i++)sum+=norm(tab[i/n][i%n]);
//          if(i<5100)
//      cerr<<i<<" "<<sum<<tab[i/n][i%n]<<endl;}
     sum=sqrt(sum);
     for(int i=0;i<n2;i++)tab[i/n][i%n]/=sum;
    }
    inline double r1(int i){
     return i%n==n/2?0:1/sqrt((i/n-n/2)*(i/n-n/2)+(i%n-n/2)*(i%n-n/2))*dx;   
    }
    void mult(){
     for(int i=0;i<n;i++){
         for(int j=0;j<n;j++){
                in[j][0]=tab[i][j].real();
                in[j][1]=tab[i][j].imag();
         }
         fftw_execute(p1);
         for(int j=0;j<n;j++){
            int k;
            if(j<n/2)k=2*j*M_PI/L;
            else k=(j-n)*2*M_PI/L;
                out[j][0]*=k*k;
                out[j][1]*=k*k;
         }
         fftw_execute(p2);
         for(int j=0;j<n;j++)ttab[i*n+j]=complex(in[j][0]/n2,in[j][1]/n2);
     }
    }
public:
    void evolve(){
        for(int ii=0;ii<10;ii++){
        mult();
        for(int i=0;i<n2;i++){tab[i/n][i%n]*=(1+r1(i)*p.dt);tab[i/n][i%n]-=ttab[i]*p.dt;
            if(i%n==0)cerr<<"\r"<<i;
        }
        normuj();
        cerr<<endl;
        
        }
    }
    void show(){
     for(int i=0;i<n2;i++){
         if(i%n==0)cout<<endl;
         cout<<(i/n-n/2)<<" "<<(i%n-n/2)<<" "<<abs(tab[i/n][i%n])<<endl;   
   
     }
        
    }
};



int main(){
    ffal psi;   
    psi.evolve();
    psi.show();
    
}