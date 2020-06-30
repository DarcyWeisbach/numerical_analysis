#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* 2分法 */
double bisection(double a,double b,double eps,FILE *fout);
/* 関数の定義 */
double f(double x);

int main(void){
  double a,b,x,h,y1,y2,eps=pow(2.0,-30.0);
	int n;
  FILE *fout;
  errno_t error;
  error = fopen_s(&fout, "output2.csv", "w");
  if(error != 0){
    printf("ファイルは見つかりません:output2.csv ￥n");
    exit(1);
  }
	printf("初期区間[a,b]を入力してください．---> a b\n");
	scanf_s("%lf%lf",&a,&b);
	printf("区間の分割数nを入力してください．---> n\n");
	scanf_s("%d",&n);

	/* 対象区間を探索しながら2分法を適用 */
	h=(b-a)/n; y1=f(a);
	for(x=a+h;x<=b;x+=h){
		y2=f(x);
		if(y1*y2<0.0){
			printf("求める答えはx=%fです.\n",bisection(x-h,x,eps,fout));
		}
		y1=y2;
	}

	return 0;
}

/* 2分法 */
double bisection(double a,double b,double eps,FILE *fout){
	double c;
  int m=1;

	do{
		c=0.5*(a+b);
		if(f(a)*f(c)<0){
			b=c;
		}
		else{
			a=c;
		}
    printf("反復回数：%3d，",m);
    printf("反復列：%+4.14f\n",c);
    fprintf(fout,"%d,%lf\n",m,c);
    m++;
	}while(fabs(b-a) >=eps); /* fabs()は絶対値を返す．「C言語入門」p.264 */
	c=0.5*(a+b);
	return c;
}

/* 関数の定義 */
double f(double x){
	return x*(x*x*(x*x-5.0)+4.0);
}
