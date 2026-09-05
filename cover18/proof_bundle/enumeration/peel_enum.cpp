// Exhaustive rooted peeling of simple disk triangulations, with exact necessary angle pruning.
// No floating-point arithmetic and no hash used as equality. See enumeration proof in PROOF_zh.md.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <cassert>
#include <cstdio>
using namespace std;
constexpr int VMAX=20,FMAX=36,TSMAX=40;
int N,B,I,FMAXOUT,CA,CC4,CC6,CD=100000; bool PRUNE=true;
using Matrix=array<array<int,VMAX>,VMAX>;
array<uint32_t,VMAX> adj;
array<uint32_t,FMAX> faces;
uint64_t calls=0,choice=0,pruned=0,complete=0,limit=0;
chrono::steady_clock::time_point t0;
struct Task {array<uint8_t,VMAX> p{}; int n=0,k=0;};
array<Task,TSMAX> tasks;
map<vector<uint32_t>,uint64_t> reps;
map<vector<uint32_t>,int> stabs;

bool addbound(Matrix&d,int i,int j,int w){
 if(w>=d[i][j])return true;
 if(d[j][i]+w<0)return false;
 // d is already a closed, consistent DBM. One new edge can be used at most once
 // in a shortest path; more traversals contain a nonnegative cycle.
 for(int a=0;a<B;a++)for(int b=0;b<B;b++)
   if(d[a][i]+w+d[j][b]<d[a][b])d[a][b]=d[a][i]+w+d[j][b];
 return true;
}
bool tighten(const Matrix&old,Matrix&d){
 d=old;if(!PRUNE)return true;
 uint32_t r2[VMAX];
 for(int u=0;u<N;u++){uint32_t z=adj[u];r2[u]=z;while(z){int v=__builtin_ctz(z);z&=z-1;r2[u]|=adj[v];}}
 for(int i=0;i<B;i++){
  uint32_t one=adj[i]|adj[(i+1)%B],two=r2[i]|r2[(i+1)%B];
  for(int j=i+1;j<B;j++){
   uint32_t ends=(1u<<j)|(1u<<((j+1)%B));int cap;
   if(one&ends)cap=CC4;else if(two&ends)cap=CC6;else continue;
   int k=j-i;if(min(k,B-k)*CA<=cap)continue;
   bool l=k*CA<CD-cap,r=(B-k)*CA<CD-cap;
   if(l&&r)return false;
   if(l&&!addbound(d,i,j,cap))return false;
   if(r&&!addbound(d,j,i,cap-CD))return false;
  }
 }
 return true;
}
vector<uint32_t> rooted(int nf,int start,int dir){
 int lab[VMAX];fill(lab,lab+N,-1);
 for(int j=0;j<B;j++)lab[(start+dir*j+B)%B]=j;
 int ef[VMAX][VMAX][2]{};int ct[VMAX][VMAX]{};
 for(int f=0;f<nf;f++){
  uint32_t z=faces[f];int a=__builtin_ctz(z);z&=z-1;int b=__builtin_ctz(z);z&=z-1;int c=__builtin_ctz(z);
  int vs[3]={a,b,c};for(int i=0;i<3;i++)for(int j=i+1;j<3;j++){
   int u=min(vs[i],vs[j]),v=max(vs[i],vs[j]);if(ct[u][v]>=2)throw runtime_error("more than two faces");ef[u][v][ct[u][v]++]=f;
  }
 }
 int a=start,b=(start+dir+B)%B;int u=min(a,b),v=max(a,b);if(ct[u][v]!=1)throw runtime_error("bad boundary root");
 array<array<int,3>,3*FMAX> stack;int sp=0;stack[sp++]={ef[u][v][0],a,b};uint64_t seen=0;int nxt=B;
 vector<uint32_t> out;
 while(sp){auto e=stack[--sp];int f=e[0],a=e[1],b=e[2];if(seen>>f&1)continue;seen|=1ULL<<f;
  uint32_t z=faces[f]^(1u<<a)^(1u<<b);if(__builtin_popcount(z)!=1)throw runtime_error("not incident");int c=__builtin_ctz(z);
  if(lab[c]<0)lab[c]=nxt++;if(lab[a]<0||lab[b]<0)throw runtime_error("unseen endpoint");
  out.push_back((1u<<lab[a])|(1u<<lab[b])|(1u<<lab[c]));
  int aa[2]={c,b},bb[2]={a,c};for(int j=0;j<2;j++){
   int x=min(aa[j],bb[j]),y=max(aa[j],bb[j]);for(int k=0;k<ct[x][y];k++)if(ef[x][y][k]!=f&&!(seen>>ef[x][y][k]&1))stack[sp++]={ef[x][y][k],aa[j],bb[j]};
  }
 }
 if(nxt!=N||out.size()!=size_t(nf))throw runtime_error("bad traversal");sort(out.begin(),out.end());return out;
}
void emit(int nf,int next){
 if(nf!=FMAXOUT||next!=N)throw runtime_error("incomplete leaf");complete++;
 vector<uint32_t> best;int stab=0;
 for(int dr:{1,-1})for(int st=0;st<B;st++){
  auto s=rooted(nf,st,dr);if(best.empty()||s<best){best=s;stab=1;}else if(s==best)stab++;
 }
 reps[best]++;stabs[best]=stab;
}
void search(int nt,int next,int nf,const Matrix&db){
 if(nt<0||nt>=TSMAX||nf<0||nf>FMAXOUT)throw runtime_error("recursion storage bound");
 calls++;
 if(limit&&calls>limit)throw runtime_error("TEST LIMIT: no complete result");
 if(calls%10000000==0)cerr<<"calls "<<calls<<" complete "<<complete<<" orbits "<<reps.size()<<" sec "<<chrono::duration<double>(chrono::steady_clock::now()-t0).count()<<"\n";
 if(nt==0){emit(nf,next);return;}
 Task cur=tasks[nt-1];int p=cur.n,k=cur.k;
 if(p==2){if(k==0)search(nt-1,next,nf,db);return;}
 if(p<3||nf>=FMAXOUT)throw runtime_error("invalid recursion");
 int u=cur.p[0],v=cur.p[1];
 // Case A: the root face's third vertex already belongs to this region's boundary.
 for(int h=2;h<p;h++){
  int w=cur.p[h];bool fresh_uw=(h!=p-1),fresh_vw=(h!=2);
  if((fresh_uw&&(adj[u]>>w&1))||(fresh_vw&&(adj[v]>>w&1)))continue;
  uint32_t au=adj[u],av=adj[v],aw=adj[w];
  adj[u]|=1u<<w;adj[v]|=1u<<w;adj[w]|=(1u<<u)|(1u<<v);faces[nf]=(1u<<u)|(1u<<v)|(1u<<w);
  choice++;Matrix d;
  if(tighten(db,d)){
   Task left,right;left.n=h;for(int j=1;j<=h;j++)left.p[j-1]=cur.p[j];right.n=p-h+1;right.p[0]=cur.p[0];for(int j=h;j<p;j++)right.p[j-h+1]=cur.p[j];
   for(int kl=0;kl<=k;kl++){
    if((left.n==2&&kl)||(right.n==2&&k-kl))continue;
    left.k=kl;right.k=k-kl;int nnt=nt-1;
    if(nnt+2>=TSMAX)throw runtime_error("region stack bound");
    if(right.n>2)tasks[nnt++]=right;
    if(left.n>2)tasks[nnt++]=left;
    search(nnt,next,nf+1,d);
   }
  }else pruned++;
  adj[u]=au;adj[v]=av;adj[w]=aw;
 }
 // Case B: the third vertex is new. First-exposure order supplies its label.
 if(k>0){
  if(next>=N||p+1>VMAX)throw runtime_error("new vertex bound");int w=next;uint32_t au=adj[u],av=adj[v],aw=adj[w];
  adj[u]|=1u<<w;adj[v]|=1u<<w;adj[w]|=(1u<<u)|(1u<<v);faces[nf]=(1u<<u)|(1u<<v)|(1u<<w);
  choice++;Matrix d;
  if(tighten(db,d)){
   Task t;t.n=p+1;t.k=k-1;t.p[0]=u;t.p[1]=w;for(int j=1;j<p;j++)t.p[j+1]=cur.p[j];tasks[nt-1]=t;search(nt,next+1,nf+1,d);
  }else pruned++;
  adj[u]=au;adj[v]=av;adj[w]=aw;
 }
 tasks[nt-1]=cur;
}
int main(int argc,char**argv){
 if(argc<8){cerr<<"usage N B A C4 C6 prefix prune[0|1] [max_calls]\n";return 2;}
 N=stoi(argv[1]);B=stoi(argv[2]);I=N-B;FMAXOUT=2*N-B-2;CA=stoi(argv[3]);CC4=stoi(argv[4]);CC6=stoi(argv[5]);string pre=argv[6];PRUNE=stoi(argv[7]);if(argc>8)limit=stoull(argv[8]);
 if(N>VMAX||B<3||B>N||FMAXOUT>FMAX||CA<=0||CA>=CD/2||CC4<=0||CC4>=CD/2||CC6<CC4||CC6>=CD/2)throw runtime_error("invalid sizes or constants");
 std::remove((pre+".json").c_str());std::remove((pre+".txt").c_str());
 for(int i=0;i<N;i++)adj[i]=1u<<i;for(int i=0;i<B;i++){adj[i]|=1u<<((i+1)%B);adj[(i+1)%B]|=1u<<i;}
 Matrix d;for(int i=0;i<B;i++)for(int j=0;j<B;j++)d[i][j]=(i==j?0:10000000);
 for(int i=0;i<B-1;i++){d[i][i+1]=CA;d[i+1][i]=0;}d[B-1][0]=CA-CD;d[0][B-1]=CD;
 for(int k=0;k<B;k++)for(int i=0;i<B;i++)for(int j=0;j<B;j++)d[i][j]=min(d[i][j],d[i][k]+d[k][j]);
 tasks[0].n=B;tasks[0].k=I;for(int j=0;j<B;j++)tasks[0].p[j]=j;t0=chrono::steady_clock::now();
 for(int i=0;i<B;i++)if(d[i][i]<0&&PRUNE){cerr<<"base infeasible\n";return 3;}
 search(1,B,0,d);
 // Every surviving orbit must appear once for every fixed boundary labeling.
 // This is a check, not the completeness premise of peeling.
 for(auto const&[v,m]:reps){int st=stabs[v];if(m*st!=uint64_t(2*B))throw runtime_error("orbit multiplicity mismatch "+to_string(m)+" * "+to_string(st));}
 ofstream out(pre+".txt");if(!out)throw runtime_error("cannot write manifest");for(auto const&[v,m]:reps){out<<stabs[v];for(auto f:v)out<<" "<<f;out<<"\n";}out.close();if(!out)throw runtime_error("manifest write failed");
 ofstream js(pre+".json");if(!js)throw runtime_error("cannot write report");js<<"{\"N\":"<<N<<",\"B\":"<<B<<",\"I\":"<<I<<",\"A\":"<<CA<<",\"C4\":"<<CC4<<",\"C6\":"<<CC6<<",\"denominator\":"<<CD<<",\"pruning\":"<<(PRUNE?"true":"false")<<",\"calls\":"<<calls<<",\"choices\":"<<choice<<",\"pruned\":"<<pruned<<",\"rooted_survivors\":"<<complete<<",\"orbits\":"<<reps.size()<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-t0).count()<<",\"complete\":true}\n";
 js.close();if(!js)throw runtime_error("report write failed");
 cout<<"COMPLETE B "<<B<<" I "<<I<<" calls "<<calls<<" roots "<<complete<<" orbits "<<reps.size()<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-t0).count()<<"\n";
}
