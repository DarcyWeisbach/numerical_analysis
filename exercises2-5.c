#include <stdio.h>
#include <stdlib.h>

int main(void) {
  double x,y;
  FILE *fout;
  errno_t error;
  error = fopen_s(&fout, "output.csv", "w");
  if(error != 0){
    printf("ファイルは見つかりません:output.csv ￥n");
    exit(1);
  }
  for(x=-3;x<=3;x+=0.2){
    y=x*x+x-2;
    fprintf(fout,"%lf,%lf\n",x,y);
  }
  fclose(fout);
}
