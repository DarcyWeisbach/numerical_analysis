#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//cl /source-charset:utf-8 exercises3-4.c /EHsc
#define N 4 /* N次正方行列 */
void input_matrix(double a[N][N],char c,FILE* fin, FILE* fout);
void input_vector(double b[N],char c,FILE* fin,FILE* fout);
void simple_gauss(double a[N][N],double b[N]);

int main(void){
    FILE *fin, *fout;
    double a[N][N], b[N];
    int i;
	errno_t error;

    if((error=fopen_s(&fin,"input_cannot.dat","r"))!=0) exit(1);
    if((error=fopen_s(&fout,"output_caonot.dat","w"))!=0) exit(1);
    input_matrix(a,'A',fin,fout); input_vector(b,'b',fin,fout);
    simple_gauss(a,b);
    fprintf(fout,"Ax=bの解は次の通りです\n");
	for(i=0;i<N;i++){ fprintf(fout,"%f\n",b[i]); }
    fclose(fin); fclose(fout);
    return 0;
}

void simple_gauss(double a[N][N],double b[N]){
    int i,j,k;
    double alpha, tmp;

    for(k=0;k<N-1;k++){
        for(i=k+1;i<N;i++){ /* 前進消去 */
            alpha=-a[i][k]/a[k][k];
            for(j=k+1;j<N;j++){
                a[i][j]=a[i][j]+alpha*a[k][j];
            }
            b[i]=b[i]+alpha*b[k];
        }
    }
    b[N-1]=b[N-1]/a[N-1][N-1]; /* 後退代入 */
    for(k=N-2;k>=0;k--){
        tmp=b[k];
        for(j=k+1;j<N;j++){
            tmp=tmp-a[k][j]*b[j];
        }
        b[k]=tmp/a[k][k];
    }
}


void input_matrix(double a[N][N],char c,FILE* fin,FILE* fout){
    int i,j;
    fprintf(fout,"行列%cは次の通りです\n",c);
    for(i=0;i<N;++i){
        for(j=0;j<N;++j){
            fscanf_s(fin,"%lf",&(a[i][j]));
            fprintf(fout,"%5.2f\t",a[i][j]);
        }
        fprintf(fout,"\n");
    }
}

void input_vector(double b[N],char c,FILE* fin,FILE* fout){
    int i;
    fprintf(fout,"ベクトル%cは次の通りです\n",c);
    for(i=0;i<N;++i){
        fscanf_s(fin,"%lf",&(b[i]));
        fprintf(fout,"%5.2f\t",b[i]);
        fprintf(fout,"\n");
    }
}
