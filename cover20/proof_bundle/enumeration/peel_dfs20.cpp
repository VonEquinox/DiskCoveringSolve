// Exhaustive root-face peeling with exact necessary angular inequalities.
// Complete finite verifier. See PROOF_zh.md for the pending-region induction.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <stdexcept>
using namespace std;
constexpr int NN=20, FF=27;
int N=20,B,I,FA;
long long SCALE=10000000, CAP[3]={873043,1822179,3019184};
bool screen=true;
uint64_t calls=0,prunes=0,complete=0,duplicates=0,invalid=0,checks=0;
uint64_t limit=0;
chrono::steady_clock::time_point start;
struct Region{ vector<int> p; int inner; };
struct State{array<uint32_t,NN> adj{};array<uint32_t,FF> faces{};int next=0,nf=0;};
struct Edge{int u,v;long long w;};
map<vector<uint32_t>,vector<uint32_t>> survivors;
int canonical_stab=0;
vector<uint32_t> canonical(State const&s){
 // Boundary labels fixed by root; internal labels assigned at first DFS encounter.
 uint32_t inc[NN][NN]{}; int tr[FF][3]{};
 for(int f=0;f<s.nf;f++){uint32_t m=s.faces[f];for(int j=0;j<3;j++){tr[f][j]=__builtin_ctz(m);m&=m-1;}for(int j=0;j<3;j++)for(int k=j+1;k<3;k++){int u=tr[f][j],v=tr[f][k];inc[u][v]|=1u<<f;inc[v][u]|=1u<<f;}}
 vector<int> profile;vector<pair<int,int>> roots;
 for(int dr:{1,-1})for(int a=0;a<B;a++){
  vector<int> p;for(int k=0;k<B;k++)p.push_back(__builtin_popcount(s.adj[(a+dr*k+B)%B]));
  if(profile.empty()||p<profile){profile=p;roots.clear();}if(p==profile)roots.emplace_back(a,dr);
 }
 vector<uint32_t> best;int stab=0;
 for(auto [a,dr]:roots){
  int mp[NN];fill(mp,mp+NN,-1);for(int j=0;j<B;j++)mp[(a+dr*j+B)%B]=j;
  int b=(a+dr+B)%B;uint32_t root=inc[a][b];if(__builtin_popcount(root)!=1)throw runtime_error("root face incidence");
  struct Q{int f,u,v;};vector<Q> st{{__builtin_ctz(root),a,b}};uint32_t seen=0;int next=B;vector<uint32_t> code;
  while(!st.empty()){
   auto q=st.back();st.pop_back();if(seen>>q.f&1)continue;seen|=1u<<q.f;
   int c=tr[q.f][0]^tr[q.f][1]^tr[q.f][2]^q.u^q.v;
   if(mp[c]<0)mp[c]=next++;if(mp[q.u]<0||mp[q.v]<0)throw runtime_error("unassigned");
   code.push_back((1u<<mp[q.u])|(1u<<mp[q.v])|(1u<<mp[c]));
   for(auto [u,v]:vector<pair<int,int>>{{c,q.u},{q.v,c}}){uint32_t fs=inc[u][v]&~seen;if(fs){if(__builtin_popcount(fs)!=1)throw runtime_error("bad incidence");st.push_back({__builtin_ctz(fs),u,v});}}
  }
  if(next!=N||(int)code.size()!=FA)throw runtime_error("canonical incomplete");sort(code.begin(),code.end());
  if(best.empty()||code<best){best=code;stab=1;}else if(best==code)stab++;
 }
 canonical_stab=stab;return best;
}
bool feasible(State const&s){
 if(!screen)return true;checks++;
 uint32_t one[NN],two[NN];for(int i=0;i<N;i++)one[i]=s.adj[i]|(1u<<i);
 for(int i=0;i<N;i++){uint32_t m=one[i],x=m;while(m){int j=__builtin_ctz(m);m&=m-1;x|=one[j];}two[i]=x;}
 Edge e[512];int ne=0;
 for(int i=0;i<B;i++){
  int j=(i+1)%B;e[ne++]={i,j,CAP[0]-(j==0?SCALE:0)};e[ne++]={j,i,j==0?SCALE:0};
 }
 for(int i=0;i<B;i++){
  uint32_t r1=one[i]|one[(i+1)%B],r2=two[i]|two[(i+1)%B];
  for(int j=i+1;j<B;j++){
   uint32_t ends=(1u<<j)|(1u<<((j+1)%B));long long cap;
   if(r1&ends)cap=CAP[1];else if(r2&ends)cap=CAP[2];else continue;
   int k=j-i;if(min(k,B-k)*CAP[0]<=cap)continue;
   bool f=k*CAP[0]<SCALE-cap,back=(B-k)*CAP[0]<SCALE-cap;
   if(f&&back)return false;
   if(f)e[ne++]={i,j,cap};else if(back)e[ne++]={j,i,cap-SCALE};
  }
 }
 long long dist[NN]{};
 for(int t=0;t<B;t++){
  bool changed=false;for(int k=0;k<ne;k++){auto a=e[k];if(dist[a.v]>dist[a.u]+a.w){dist[a.v]=dist[a.u]+a.w;changed=true;}}
  if(!changed)return true;
 }
 return false;
}
void progress(){
 if(calls%1000000==0){cerr<<"B "<<B<<" calls "<<calls<<" prunes "<<prunes<<" full "<<complete<<" distinct "<<survivors.size()<<" sec "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n";}
 if(limit&&calls>=limit)throw runtime_error("TEST LIMIT, not complete");
}
void visit(State const&s,vector<Region> const&pending){
 calls++;progress();
 if(pending.empty()){
  if(s.next!=N||s.nf!=FA)throw runtime_error("completion size");complete++;
  vector<uint32_t> code=canonical(s);vector<uint32_t> fa(s.faces.begin(),s.faces.begin()+s.nf);sort(fa.begin(),fa.end());
  if(!survivors.emplace(code,fa).second)duplicates++;
  return;
 }
 Region reg=pending.back();auto p=reg.p;int m=p.size(),ni=reg.inner;
 vector<Region> rest=pending;rest.pop_back();
 if(m==2){if(ni==0)visit(s,rest);else invalid++;return;}
 if(m<3||ni<0||m+ni>N)throw runtime_error("invalid region");
 int a=p[0],b=p[1];
 for(int k=2;k<m;k++){
  int c=p[k];if(k>2&&(s.adj[b]>>c&1)){invalid++;continue;}if(k<m-1&&(s.adj[a]>>c&1)){invalid++;continue;}
  State ns=s;bool change=false;
  for(auto [u,v]:vector<pair<int,int>>{{a,c},{b,c}}){if(!(ns.adj[u]>>v&1))change=true;ns.adj[u]|=1u<<v;ns.adj[v]|=1u<<u;}
  if(ns.nf>=FA)throw runtime_error("too many faces");ns.faces[ns.nf++]=(1u<<a)|(1u<<b)|(1u<<c);
  if(change&&!feasible(ns)){prunes++;continue;}
  vector<int> left{c};for(int j=1;j<k;j++)left.push_back(p[j]);
  vector<int> right{a};for(int j=k;j<m;j++)right.push_back(p[j]);
  for(int nleft=0;nleft<=ni;nleft++){
   if((left.size()==2&&nleft)||(right.size()==2&&ni-nleft))continue;
   auto work=rest;work.push_back({right,ni-nleft});work.push_back({left,nleft});visit(ns,work);
  }
 }
 if(ni>0){
  State ns=s;int v=ns.next++;if(v>=N)throw runtime_error("too many vertices");
  ns.adj[v]|=(1u<<a)|(1u<<b);ns.adj[a]|=1u<<v;ns.adj[b]|=1u<<v;
  ns.faces[ns.nf++]=(1u<<a)|(1u<<b)|(1u<<v);
  if(!feasible(ns)){prunes++;return;}
  vector<int> poly{a,v};for(int j=1;j<m;j++)poly.push_back(p[j]);auto work=rest;work.push_back({poly,ni-1});visit(ns,work);
 }
}
int main(int argc,char**argv){try{
 if(argc<4)throw runtime_error("usage peel N B output [limit] [noscreen]");N=stoi(argv[1]);B=stoi(argv[2]);I=N-B;FA=2*N-B-2;string path=argv[3];if(N>NN||B<3||I<0||FA>FF)throw runtime_error("parameters");limit=argc>4?stoull(argv[4]):0;screen=!(argc>5&&string(argv[5])=="noscreen");
 State s;s.next=B;for(int i=0;i<B;i++){int j=(i+1)%B;s.adj[i]|=1u<<j;s.adj[j]|=1u<<i;}
 vector<int> p;for(int i=0;i<B;i++)p.push_back(i);start=chrono::steady_clock::now();
 if(feasible(s))visit(s,{{p,I}});else prunes++;
 ofstream f(path);if(!f)throw runtime_error("output open failed");f<<"{\"version\":1,\"screen\":"<<(screen?"true":"false")<<",\"angle_scale\":"<<SCALE<<",\"angle_caps\":["<<CAP[0]<<","<<CAP[1]<<","<<CAP[2]<<"],\"N\":"<<N<<",\"B\":"<<B<<",\"I\":"<<I<<",\"calls\":"<<calls<<",\"prunes\":"<<prunes<<",\"complete_rooted\":"<<complete<<",\"distinct_survivors\":"<<survivors.size()<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<",\"survivors\":[";bool first=true;
 for(auto const&[code,fa]:survivors){if(!first)f<<",";first=false;f<<"[";for(int j=0;j<(int)code.size();j++){if(j)f<<",";f<<code[j];}f<<"]";}f<<"]}\n";f.close();if(!f)throw runtime_error("output write failed");
 cerr<<"COMPLETE N "<<N<<" B "<<B<<" calls "<<calls<<" prunes "<<prunes<<" full "<<complete<<" distinct "<<survivors.size()<<" sec "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n";
}catch(exception const&e){cerr<<"REJECTED: "<<e.what()<<"\n";return 1;}}
