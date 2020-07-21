#include <stdio.h>
#include <math.h>

/* 2���@ */
double bisection(double a,double b,double eps);
/* �֐��̒�` */
double f(double x);

int main(void){
    double a,b,x,h,y1,y2,eps=pow(2.0,-30.0);
	int n;

	printf("�������[a,b]����͂��Ă��������D---> a b\n");
	scanf_s("%lf%lf",&a,&b);
	printf("��Ԃ̕�����n����͂��Ă��������D---> n\n");
	scanf_s("%d",&n);

	/* �Ώۋ�Ԃ�T�����Ȃ���2���@��K�p */
	h=(b-a)/n; y1=f(a);
	for(x=a+h;x<=b;x+=h){
		y2=f(x);
		if(y1*y2<0.0){
			printf("���߂铚����x=%f�ł�.\n",bisection(x-h,x,eps));
		}
		y1=y2;
	}

	return 0;
}

/* 2���@ */
double bisection(double a,double b,double eps){
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
    printf("�����񐔁F%3d�C",m);
    printf("������F%+4.14f\n",c);
    m++;
	}while(fabs(b-a) >=eps); /* fabs()�͐�Βl��Ԃ��D�uC�������vp.264 */
	c=0.5*(a+b);
	return c;
}

/* �֐��̒�` */
double f(double x){
	return x*(x*x*(x*x-5.0)+4.0);
}
