#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//cl /source-charset:shift-jis exercises4-3.cpp /EHsc

#define N 5 /* THis exercises4-3 ,matric 5x5*/
void input_matrix(double a[N][N],char c,FILE* fin, FILE* fout);
void input_vector(double b[N],char c,FILE* fin,FILE* fout);
void lu_decomp(double a[N][N],int p[N]);
void lu_solve(double a[N][N],double b[N],int p[N]);

int main(void){
    FILE *fin, *fout;
    double a[N][N], b[N];
    int i, p[N]; /* p[1...N-2]を利用, p[N-1]は未使用 */
	errno_t error;
    if((error=fopen_s(&fin,"input_exercises4-3.dat","r"))!=0) exit(1);
    if((error=fopen_s(&fout,"output_exercises4-3.dat","w"))!=0) exit(1);
    input_matrix(a,'A',fin,fout); input_vector(b,'b',fin,fout);
	lu_decomp(a,p);
    lu_solve(a,b,p);
    fprintf(fout,"Ax=bの解は次の通りです\n");
	for(i=0;i<N;i++){ fprintf(fout,"%f\n",b[i]); }
    fclose(fin); fclose(fout);
    return 0;
}

void lu_decomp(double a[N][N],int p[N]){
    int i,j,k,ip;
    double alpha, tmp;
	double amax, eps=pow(2.0,-50.0); /* eps=2^{-50}*/
    for(k=0;k<N-1;k++){
		amax=fabs(a[k][k]); ip=k; /* ピボットの選択 */
		for(i=k+1;i<N;i++){
			if(fabs(a[i][k])>amax){
				amax=fabs(a[i][k]); ip=i;
			}
		}
		if(amax<eps) { /* 正則性の判定 */
			printf("入力した行列は正則ではない!!\n"); exit(1);
		}
		p[k]=ip;
		if(ip!=k){
			for(j=k;j<N;j++){
				tmp=a[k][j]; a[k][j]=a[ip][j]; a[ip][j]=tmp;
			}
		}
        for(i=k+1;i<N;i++){ /* 前進消去 */
            alpha=a[i][k]/a[k][k];
			a[i][k]=alpha;
            for(j=k+1;j<N;j++){
                a[i][j]=a[i][j]-alpha*a[k][j];
            }
        }
    }
}

void lu_solve(double a[N][N],double b[N],int p[N]){
	int i,j,k;
	double tmp;

	for(k=0;k<N-1;k++){
		tmp=b[k]; b[k]=b[p[k]]; b[p[k]]=tmp; /* 右辺の行変換 */
		for(i=k+1;i<N;i++){  /* 前進代入 */
			b[i]=b[i]-a[i][k]*b[k];
		}
	}
	b[N-1]=b[N-1]/a[N-1][N-1];
	for(k=N-2;k>=0;k--){ /* 後退代入 */
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
