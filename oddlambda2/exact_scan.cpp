#include <boost/multiprecision/cpp_int.hpp>
using boost::multiprecision::cpp_int;
/* Exact signed fraction-free congruence for integer symmetric matrices.
 * If C = q*(A_alpha(G) - (m-1+alpha)*I), n=2m+1,
 * its positive/zero/negative inertia is identical to the shifted spectrum.
 * A_alpha = alpha*D+(1-alpha)*A, alpha=p/q.
 * All arithmetic uses unbounded exact integers (boost::multiprecision::cpp_int).
 * No eigensolver and no tolerance is used.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#define N 9
static const int frac[][2]={{0,1},{1,20},{1,10},{1,8},{1,6},{1,5},{1,4},{1,3},{2,5},{9,20},{12,25},{49,100},{499,1000}};
#define F (sizeof(frac)/sizeof(frac[0]))
typedef cpp_int Big;
static int sgn(const Big &x){return (x>0)-(x<0);}
static void decode(const char *s,int n,int a[N][N]){
 if(s[0]!=n+63){fprintf(stderr,"bad graph6 length/n\n");exit(2);}
 memset(a,0,sizeof(int)*N*N);int k=0;
 for(int v=1;v<n;++v)for(int u=0;u<v;++u){
  int b=((unsigned char)s[1+k/6]-63)>>(5-k%6)&1;
  a[u][v]=a[v][u]=b;++k;
 }
}
static void iswap(Big *a,Big *b){a->swap(*b);}
static void inertia(int n,const int a[N][N],int p,int q,int *positive,int *zero){
 Big c[N][N]={0};int m=(n-1)/2;
 for(int i=0;i<n;++i){int degree=0;for(int j=0;j<n;++j)degree+=a[i][j];
  for(int j=0;j<n;++j)c[i][j]=(i==j)?(Big)p*(degree-1)-(Big)q*(m-1):(Big)(q-p)*a[i][j];
 }
 Big previous=1;
 int pos=0,z=0;
 for(int k=0;k<n;++k){
  int pivot=-1;
  for(int i=k;i<n;++i)if(c[i][i]!=0){pivot=i;break;}
  if(pivot<0){
   int i0=-1,j0=-1;
   for(int i=k;i<n&&i0<0;++i)for(int j=i+1;j<n;++j)if(c[i][j]!=0){i0=i;j0=j;break;}
   if(i0<0){z+=n-k;break;}
   /* congruence: variable i0 is replaced by variable i0+variable j0 */
   c[i0][i0]=2*c[i0][j0]+c[j0][j0];
   for(int j=k;j<n;++j)if(j!=i0){c[i0][j]+=c[j0][j];c[j][i0]=c[i0][j];}
   pivot=i0;
  }
  if(pivot!=k){for(int j=0;j<n;++j)iswap(&c[k][j],&c[pivot][j]);
   for(int i=0;i<n;++i)iswap(&c[i][k],&c[i][pivot]);}
  Big pk=c[k][k];if(sgn(pk)==sgn(previous))++pos;
  for(int i=k+1;i<n;++i)for(int j=i;j<n;++j){
   Big num=pk*c[i][j]-c[i][k]*c[k][j];
   if(num%previous){fprintf(stderr,"fraction-free exact division failed at n=%d p=%d q=%d k=%d i=%d j=%d\n",n,p,q,k,i,j);exit(3);}
   c[i][j]=num/previous;c[j][i]=c[i][j];
  }
  previous=pk;
 }
 *positive=pos;*zero=z;
}
int main(int argc,char **argv){
 if(argc!=3){fprintf(stderr,"usage: exact_scan n graphsN.g6\n");return 2;}
 int n=atoi(argv[1]);if(n!=5&&n!=7&&n!=9)return 2;
 FILE *f=fopen(argv[2],"r");if(!f){perror(argv[2]);return 2;}
 size_t total=0,violation[F]={0},ties[F]={0},nonpositive[F]={0};
 char first[F][32]={0},line[128];int a[N][N];
 while(fgets(line,sizeof(line),f)){
  if(line[0]=='\n')continue;
  int len=(int)strcspn(line,"\r\n");line[len]='\0';
  if(len!=1+(n*(n-1)/2+5)/6){fprintf(stderr,"bad line: %s\n",line);return 2;}
  decode(line,n,a);++total;
  for(size_t h=0;h<F;++h){int pos,z;inertia(n,a,frac[h][0],frac[h][1],&pos,&z);
   if(pos>=2)++violation[h];
   if(pos==0)++nonpositive[h];
   if((pos==1&&z>=1)||(pos==0&&z>=2)){
    ++ties[h];if(!first[h][0])strcpy(first[h],line);
   }
  }
 }
 fclose(f);
 printf("n=%d graph_count=%zu expected=%d\n",n,total,n==5?34:n==7?1044:274668);
 if(total!=(size_t)(n==5?34:n==7?1044:274668))return 3;
 for(size_t h=0;h<F;++h){
  printf("n=%d alpha=%d/%d bound=%d+%d/%d violations=%zu equality=%zu pos0=%zu first=%s\n",
   n,frac[h][0],frac[h][1],(n-3)/2,frac[h][0],frac[h][1],
   violation[h],ties[h],nonpositive[h],first[h][0]?first[h]:"-");
  if(violation[h])return 4;
 }
 return 0;
}
