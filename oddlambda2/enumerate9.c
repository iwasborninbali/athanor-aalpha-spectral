/* Exact orderly generation of unlabeled graphs up to order 9.
 * Start with the networkx seven-vertex atlas, supplied as graph6 lines.
 * Every graph of order k+1 is an extension of a graph of order k;
 * canon() takes the minimum binary adjacency encoding over all vertex
 * permutations preserving a *canonical* vertex-invariant color signature.
 * Sorting the signature classes is isomorphism-invariant, while exploring
 * all permutations inside classes makes the encoding complete.
 * The signature is (degree, neighbor-degree histogram).
 * No nauty, floating point, or other paid services are used.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#define H8 (1<<15)
#define H9 (1<<20)
#define EMPTY UINT64_MAX
static uint64_t hash8[H8],hash9[H9],graphs8[12346],graphs9[274668];
static size_t count8=0,count9=0;
static uint16_t adj[9];
static int ncur,order[9];
static uint64_t group[9];
static uint64_t sig[9],best;
static int found;
static uint64_t sigv(int v,int n) {
 int d=__builtin_popcount((unsigned)adj[v]);
 uint64_t s=(uint64_t)d<<36;
 for(int w=0;w<n;++w) if(adj[v]>>w&1) {
  int dw=__builtin_popcount((unsigned)adj[w]);
  s += (uint64_t)1 << (4*dw);
 }
 return s;
}
static void recurse(int pos,int used,uint64_t code) {
 if (pos==ncur) {if(!found||code<best){best=code;found=1;} return;}
 for(int v=0;v<ncur;++v) if(!(used>>v&1) && sig[v]==group[pos]) {
  uint64_t c=code;
  for(int j=0;j<pos;++j) c=(c<<1)|((adj[v]>>order[j])&1);
  if(found) { int unfilled=(ncur*(ncur-1)-pos*(pos+1))/2;
    if (c > (best>>unfilled)) continue; }
  order[pos]=v;recurse(pos+1,used|(1<<v),c);
 }
}
static uint64_t canon(int n) {
 ncur=n;
 for(int i=0;i<n;++i){sig[i]=sigv(i,n);group[i]=sig[i];}
 for(int i=1;i<n;++i) for(int j=i;j>0&&group[j-1]>group[j];--j) {
  uint64_t t=group[j];group[j]=group[j-1];group[j-1]=t;
 }
 best=EMPTY;found=0;recurse(0,0,0);
 if(!found){fprintf(stderr,"canonicalization failure\n");exit(2);}
 return best;
}
static void decode_code(uint64_t code,int n) {
 memset(adj,0,sizeof(adj));int k=n*(n-1)/2;
 for(int v=1;v<n;++v)for(int u=0;u<v;++u) {
  if((code>>(--k))&1){adj[u]|=(1<<v);adj[v]|=(1<<u);}
 }
}
static void decode_g6(const char *s,int n){
 if((unsigned char)s[0]!=(unsigned char)(n+63)) {fprintf(stderr,"bad graph6 order\n");exit(2);}
 memset(adj,0,sizeof(adj));int k=0;
 for(int v=1;v<n;++v) for(int u=0;u<v;++u){
  int b=((unsigned char)s[1+k/6]-63)>>(5-k%6)&1;
  if(b){adj[u]|=(1<<v);adj[v]|=(1<<u);}k++;
 }
}
static void encode_g6(FILE *f,uint64_t code,int n){
 decode_code(code,n);
 fputc(n+63,f);int k=0,part=0;
 for(int v=1;v<n;++v)for(int u=0;u<v;++u){
  part=(part<<1)|((adj[u]>>v)&1);
  if(++k%6==0){fputc(63+part,f);part=0;}
 }
 if(k%6)fputc(63+(part<<(6-k%6)),f);
 fputc('\n',f);
}
static uint64_t hash_key(uint64_t x){
 x^=x>>30;x*=UINT64_C(0xbf58476d1ce4e5b9);
 x^=x>>27;x*=UINT64_C(0x94d049bb133111eb);return x^(x>>31);
}
static int insert(uint64_t *tab,size_t cap,uint64_t key){
 size_t i=hash_key(key)&(cap-1);
 for(;;){if(tab[i]==EMPTY){tab[i]=key;return 1;}
  if(tab[i]==key)return 0;i=(i+1)&(cap-1);}
}
static int cmp64(const void *a,const void *b){uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b;return (x>y)-(x<y);}
int main(int argc,char **argv){
 if(argc!=3){fprintf(stderr,"usage: enumerate9 atlas7.g6 output9.g6\n");return 2;}
 FILE *f=fopen(argv[1],"r");if(!f){perror(argv[1]);return 2;}
 memset(hash8,255,sizeof(hash8));memset(hash9,255,sizeof(hash9));
 char line[128];size_t seven=0;
 while(fgets(line,sizeof(line),f)){
  if(line[0]=='#'||line[0]=='\n')continue;decode_g6(line,7);++seven;
  for(int sub=0;sub<128;++sub){
   for(int v=0;v<7;++v)if(sub>>v&1)adj[v]|=(1<<7);
   adj[7]=sub;
   uint64_t code=canon(8);
   if(insert(hash8,H8,code)){if(count8>=12346){fprintf(stderr,"unexpected n8 count\n");return 2;}graphs8[count8++]=code;}
   for(int v=0;v<7;++v)adj[v]&=~(1<<7);
   adj[7]=0;
  }
 }
 fclose(f);
 fprintf(stderr,"atlas7=%zu n8=%zu\n",seven,count8);
 if(seven!=1044||count8!=12346)return 3;
 qsort(graphs8,count8,sizeof(uint64_t),cmp64);
 for(size_t i=0;i<count8;++i){
  decode_code(graphs8[i],8);
  for(int sub=0;sub<256;++sub){
   for(int v=0;v<8;++v)if(sub>>v&1)adj[v]|=(1<<8);
   adj[8]=sub;
   uint64_t code=canon(9);
   if(insert(hash9,H9,code)){
    if(count9>=274668){fprintf(stderr,"unexpected n9 count\n");return 4;}
    graphs9[count9++]=code;
   }
   for(int v=0;v<8;++v)adj[v]&=~(1<<8);
   adj[8]=0;
  }
 }
 fprintf(stderr,"n9=%zu\n",count9);
 if(count9!=274668)return 5;
 qsort(graphs9,count9,sizeof(uint64_t),cmp64);
 FILE *out=fopen(argv[2],"w");if(!out){perror(argv[2]);return 6;}
 for(size_t i=0;i<count9;++i)encode_g6(out,graphs9[i],9);
 fclose(out);
 return 0;
}
