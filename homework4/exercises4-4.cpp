#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//cl /source-charset:utf-8 exercises4-4.cpp /EHsc
#define N 4 /* N次正方行列 */
void input_matrix(double a[N][N],char c,FILE* fin, FILE* fout);
void input_vector(double b[N],char c,FILE* vector,FILE* fout);
void lu_decompose(double a[N][N],int p[N]);//LU分解
void lu_solve(double a[N][N],double b[N],int p[N]);
void make_identity_matrix(double identity_matrix[N][N]);/*単位行列つくる*/
void change_to_transposed_matrix(double matrix[N][N]);
/*転置行列　破壊を伴う(破壊的作業を伴う)*/
void make_inverse_matrix(double upper_triangular_matrix[N][N],double identity_matrix[N][N],double inverse_matrix[N][N],int p[N]);
/*逆行列を求める．引数は名前をみればわかるようにしてあります*/
int main(void){
    FILE *fmatrix, *fout,*fvector;
    double a[N][N], b[N],I[N][N],a_inv[N][N];/*逆行列と単位行列定義*/
    int i,j, p[N]; /* p[1...N-2]を利用, p[N-1]は未使用 */
	  errno_t error;
    if((error=fopen_s(&fvector,"input_vector.dat","r"))!=0)exit(1);
    if((error=fopen_s(&fmatrix,"input_regular_matrix.dat","r"))!=0) exit(1);
    if((error=fopen_s(&fout,"output_exercises4-4.dat","w"))!=0) exit(1);
    make_identity_matrix(I);
    input_matrix(a,'A',fmatrix,fout);
    input_vector(b,'b',fvector,fout);
  	lu_decompose(a,p);
    lu_solve(a,b,p);
    change_to_transposed_matrix(I);
    /*単位行列だから転置する必要はないけど，本来列を入れるところに行を入れるため注意*/
    make_inverse_matrix(a,I,a_inv,p);
    change_to_transposed_matrix(a_inv);
    fprintf(fout,"Ax=bの解は次の通りです\n");
  	for(i=0;i<N;i++){
      fprintf(fout,"%f\n",b[i]);
    }
    fprintf(fout, "行列Aの逆行列は次の通りです\n" );
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
          fprintf(fout,"%5.2f ",a_inv[i][j]);
        }
        fprintf(fout, "\n");
      }
    fclose(fmatrix);
    fclose(fvector);
    fclose(fout);
    return 0;
}
void make_inverse_matrix(double upper_triangular_matrix[N][N],double identity_matrix[N][N],double inverse_matrix[N][N],int p[N]){
  int i,j;
  double column_vector[N];//列ベクトル
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      column_vector[j]=identity_matrix[i][j];
    }
    lu_solve(upper_triangular_matrix,column_vector,p);
    for(j=0;j<N;j++){
      inverse_matrix[i][j]=column_vector[j];
    }
  }
}

void make_identity_matrix(double identity_matrix[N][N]){
  int i,j;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if(i==j){
        identity_matrix[i][j]=1;//対角成分が1
      }
      else{
        identity_matrix[i][j]=0;
      }
    }
  }
}

void change_to_transposed_matrix(double matrix[N][N]){
  int i,j;
  double tmp;
  for(i=0;i<N;i++){
    for(j=i+1;j<N;j++){
      tmp=matrix[i][j];
      matrix[i][j]=matrix[j][i];
      matrix[j][i]=tmp;
    }
  }
}

void lu_decompose(double a[N][N],int p[N]){
    int i,j,k,ip;
    double alpha, tmp;
	  double amax, eps=pow(2.0,-50.0); /* eps=2^{-50}とする */
    for(k=0;k<N-1;k++){
      amax=fabs(a[k][k]);
      ip=k; /* ピボットの選択 */
		for(i=k+1;i<N;i++){
      if(fabs(a[i][k])>amax){
				amax=fabs(a[i][k]);
        ip=i;
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

void input_vector(double b[N],char c,FILE* vector,FILE* fout){
    int i;
    fprintf(fout,"ベクトル%cは次の通りです\n",c);
    for(i=0;i<N;++i){
        fscanf_s(vector,"%lf",&(b[i]));
        fprintf(fout,"%5.2f\t",b[i]);
        fprintf(fout,"\n");
    }
}
