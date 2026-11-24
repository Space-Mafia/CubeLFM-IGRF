/* Electrostatic projection for the optional single-foot wall return.
 * The input potential is sampled at traced footprints of a regular volume
 * cache. The output is a single potential for the gas and CT wall stencils.
 * This code does not alter the ionospheric solution or identify open tubes.
 */
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <mkl_dfti.h>

enum { NR=8, NT=91, NP=180, COUNT=NR*NT*NP };

static int index3(int i, int j, int k) {
    return (i*NT+j)*NP+k;
}

/* -1 denotes a missing radial/meridional neighbor. Longitude is periodic. */
static int neighbor(int site, int axis, int stride) {
    int offset = axis==0 ? NT*NP : axis==1 ? NP : 1;
    int n = axis==0 ? NR : axis==1 ? NT : NP;
    int q = (site/offset)%n + stride;
    if(axis==2) return site-site%NP+(q+NP)%NP;
    return q<0 || q>=n ? -1 : site+stride*offset;
}

/* Difference along increasing coordinate on each sheet, extrapolated to the
 * crossing using the nearest within-sheet segment on each available side. */
static int crossing_target(const double *v, const signed char *branch,
                           int a, int b, int axis, double *target) {
    int before=neighbor(a,axis,-1), after=neighbor(b,axis,1), count=0;
    double sum=0;
    if(before>=0 && branch[before]==branch[a]) { sum+=v[a]-v[before]; ++count; }
    if(after>=0 && branch[after]==branch[b]) { sum+=v[after]-v[b]; ++count; }
    if(!count) return -1;
    *target=sum/count;
    return 0;
}

/* Solve the real normal equations with a periodic FFT in longitude,
 * a Neumann cosine basis in radius and a tridiagonal solve in colatitude.
 * Coincident pole copies have one voltage. We pin only the constant mode. */
static int solve(double complex *b, const double *wr, const double *wt, const double *wp) {
    DFTI_DESCRIPTOR_HANDLE fft=NULL;
    int status=1;
    if(DftiCreateDescriptor(&fft,DFTI_DOUBLE,DFTI_COMPLEX,1,(MKL_LONG)NP) ||
       DftiSetValue(fft,DFTI_NUMBER_OF_TRANSFORMS,(MKL_LONG)(NR*NT)) ||
       DftiSetValue(fft,DFTI_INPUT_DISTANCE,(MKL_LONG)NP) ||
       DftiSetValue(fft,DFTI_OUTPUT_DISTANCE,(MKL_LONG)NP) ||
       DftiSetValue(fft,DFTI_BACKWARD_SCALE,1./NP) ||
       DftiCommitDescriptor(fft) || DftiComputeForward(fft,b)) goto done;
    double pi=acos(-1.), basis[NR][NR];
    for(int m=0;m<NR;m++) for(int i=0;i<NR;i++)
        basis[m][i]=sqrt((m?2.:1.)/NR)*cos(pi*m*(i+.5)/NR);
    for(int k=0;k<NP;k++) {
        double complex modes[NR][NT];
        for(int m=0;m<NR;m++) {
            double lr=4*pow(sin(.5*pi*m/NR),2), lp=4*pow(sin(pi*k/NP),2);
            int first=(k || !m)?1:0, last=k?NT-2:NT-1;
            double pivot[NT];
            for(int j=0;j<NT;j++) {
                double complex r=0;
                for(int i=0;i<NR;i++) r+=basis[m][i]*b[index3(i,j,k)];
                modes[m][j]=j<first||j>last ? 0 : r;
            }
            for(int j=first;j<=last;j++) {
                double south=j<NT-1?wt[j]:0, north=j?wt[j-1]:0;
                double d=lr*wr[j]+lp*wp[j]+south+north;
                if(j>first) {
                    d-=north*north/pivot[j-1];
                    modes[m][j]+=north*modes[m][j-1]/pivot[j-1];
                }
                if(!(d>0)) goto done;
                pivot[j]=d;
            }
            for(int j=last;j>=first;j--) {
                if(j<last) modes[m][j]+=wt[j]*modes[m][j+1];
                modes[m][j]/=pivot[j];
            }
        }
        for(int i=0;i<NR;i++) for(int j=0;j<NT;j++) {
            double complex r=0;
            for(int m=0;m<NR;m++) r+=basis[m][i]*modes[m][j];
            b[index3(i,j,k)]=r;
        }
    }
    if(!DftiComputeBackward(fft,b)) status=0;
done:
    if(fft) DftiFreeDescriptor(&fft);
    return status;
}

int cube_project_branch_voltage(double *v, const signed char *branch, double low, double high) {
    if(!v || !branch || !(low>0 && high>low)) return 1;
    double complex *rhs=calloc(COUNT,sizeof(*rhs));
    if(!rhs) return 1;
    double pi=acos(-1.), ds=pi/(NT-1), dr=(high-low)/(NR-1);
    double wr[NT],wt[NT],wp[NT],area[NT];
    int cuts=0, status=1;
    /* Mid-shell radial metric and exact latitude-band area. */
    for(int j=0;j<NT;j++) {
        double th=j*ds;
        area[j]=(cos(fmax(0,th-ds*.5))-cos(fmin(pi,th+ds*.5)))/ds;
        wr[j]=area[j]*pow(.5*(low+high)*ds/dr,2);
        wt[j]=sin(th+.5*ds);
        wp[j]=(j==0 || j==NT-1)?0:ds*ds/(pow(2*pi/NP,2)*sin(th));
    }
    for(int a=0;a<COUNT;a++) {
        if(!isfinite(v[a]) || (branch[a]!=1 && branch[a]!=-1)) goto done;
        int j=(a/NP)%NT;
        for(int d=0;d<3;d++) {
            int b=neighbor(a,d,1);
            double w=d==0?wr[j]:d==1?wt[j]:wp[j];
            if(b<0 || w==0 || branch[a]==branch[b]) continue;
            double target;
            if(crossing_target(v,branch,a,b,d,&target)) goto done;
            double change=w*(target-(v[b]-v[a]));
            rhs[a]-=change; rhs[b]+=change; ++cuts;
        }
    }
    if(cuts) {
        if(solve(rhs,wr,wt,wp)) goto done;
        double sum=0, norm=0;
        for(int a=0;a<COUNT;a++) {
            double weight=area[(a/NP)%NT];
            if(!isfinite(creal(rhs[a]))) goto done;
            sum+=weight*creal(rhs[a]); norm+=weight;
        }
        for(int a=0;a<COUNT;a++) v[a]+=creal(rhs[a])-sum/norm;
    }
    status=0;
done:
    free(rhs);
    return status;
}
