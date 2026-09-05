// Independent exact audit of a full D_B x S_I census, N=15.
// No floating point, no reliance on the enumerator's canonical labels or flags.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
constexpr int N=15, MF=18;
using Face=array<int,3>;
using Code=array<uint16_t,MF>;
int B,I,FN;
unsigned long long fac(int n){unsigned long long v=1;for(int k=2;k<=n;k++)v*=k;return v;}
int pop(unsigned x){return __builtin_popcount(x);}
void require(bool p,string s){if(!p)throw runtime_error(s);}
struct Complex {
 array<unsigned,N> adj{};
 array<array<unsigned,N>,N> link{};
 array<array<unsigned char,N>,N> count{};
 vector<Face> faces;
 array<int,N> deg{};
 explicit Complex(vector<Face> fs):faces(move(fs)) {
  require((int)faces.size()==FN,"face count");
  auto sf=faces;sort(sf.begin(),sf.end());require(adjacent_find(sf.begin(),sf.end())==sf.end(),"duplicate face");
  for(auto t:faces){int a=t[0],b=t[1],c=t[2];require(0<=a&&a<b&&b<c&&c<N,"vertex range");
   for(int k=0;k<3;k++){int u=t[k],v=t[(k+1)%3],w=t[(k+2)%3];adj[u]|=1u<<v;adj[u]|=1u<<w;link[u][v]|=1u<<w;link[u][w]|=1u<<v;count[min(u,v)][max(u,v)]++;}
  }
  int E=0;
  for(int u=0;u<N;u++){
   require(adj[u]!=0,"unused vertex");deg[u]=pop(adj[u]);
   for(int v=u+1;v<N;v++){
    bool border=(u<B&&v<B&&(v==u+1||(u==0&&v==B-1)));
    if(border)require(count[u][v]==1,"missing boundary edge");
    if(count[u][v]){E++;require(count[u][v]==(border?1:2),"edge incidence");}
   }
  }
  require(N-E+FN==1,"Euler characteristic");
  unsigned seen=1, todo=1;
  while(todo){int u=__builtin_ctz(todo);todo&=todo-1;unsigned fresh=adj[u]&~seen;seen|=fresh;todo|=fresh;}
  require(seen==(1u<<N)-1,"complex disconnected");
  for(int v=0;v<N;v++){
   unsigned neighbors=adj[v];
   for(int u=0;u<N;u++)if(neighbors>>u&1){
    int expected=(v<B&&(u==(v+B-1)%B||u==(v+1)%B))?1:2;
    require(pop(link[v][u])==expected,"link degrees/endpoints");
   }
   seen=neighbors&-neighbors;todo=seen;
   while(todo){int u=__builtin_ctz(todo);todo&=todo-1;unsigned fresh=link[v][u]&~seen;seen|=fresh;todo|=fresh;}
   require(seen==neighbors,"link disconnected");
  }
 }
 // Independent invariant: boundary degree sequence only. Then internal keys
 // are (boundary-neighbor mask, total degree), not the enumerator's keys.
 pair<Code,unsigned> canonical()const{
  vector<array<int,N>> maps;vector<int> minseq;bool seqhave=false;
  for(int rev=0;rev<2;rev++)for(int sh=0;sh<B;sh++){
   array<int,N> m{};vector<int> seq(B);
   for(int u=0;u<B;u++){m[u]=rev?(sh-u+B)%B:(sh+u)%B;seq[m[u]]=deg[u];}
   if(!seqhave||seq<minseq){seqhave=true;minseq=seq;maps.assign(1,m);}else if(seq==minseq)maps.push_back(m);
  }
  Code best{};bool have=false;unsigned mult=0;
  for(auto base:maps){
   vector<pair<pair<unsigned,int>,int>> keys;
   for(int v=B;v<N;v++){unsigned mask=0;for(int u=0;u<B;u++)if(adj[v]>>u&1)mask|=1u<<base[u];keys.push_back({{mask,deg[v]},v});}
   sort(keys.begin(),keys.end());vector<int> perm;for(auto k:keys)perm.push_back(k.second);
   vector<pair<int,int>> groups;
   for(int a=0;a<I;){int b=a+1;while(b<I&&keys[a].first==keys[b].first)b++;groups.push_back({a,b});a=b;}
   bool more=true;
   while(more){
    auto m=base;for(int k=0;k<I;k++)m[perm[k]]=B+k;
    Code z{};fill(z.begin(),z.end(),65535);
    for(int k=0;k<FN;k++){auto t=faces[k];for(int &v:t)v=m[v];sort(t.begin(),t.end());z[k]=uint16_t(t[0]*N*N+t[1]*N+t[2]);}
    sort(z.begin(),z.begin()+FN);
    if(!have||z<best){have=true;best=z;mult=1;}else if(z==best)mult++;
    more=false;for(auto [a,b]:groups)if(next_permutation(perm.begin()+a,perm.begin()+b)){more=true;break;}
   }
  }
  require(have&&mult,"empty orbit");return {best,mult};
 }
 // Direct stabilizer traversal for sparse independent spot checks of pruning.
 // Not needed for completeness: canonical's mathematical pruning is exact.
 unsigned brute_stab()const{
  vector<Face> orig=faces;sort(orig.begin(),orig.end());unsigned count=0;
  for(int rev=0;rev<2;rev++)for(int sh=0;sh<B;sh++){
   vector<int> p(I);iota(p.begin(),p.end(),B);
   do{array<int,N> m{};for(int u=0;u<B;u++)m[u]=rev?(sh-u+B)%B:(sh+u)%B;for(int u=0;u<I;u++)m[B+u]=p[u];vector<Face> fs=faces;
    for(auto &t:fs){for(int &v:t)v=m[v];sort(t.begin(),t.end());}sort(fs.begin(),fs.end());if(fs==orig)count++;
   }while(next_permutation(p.begin(),p.end()));
  }
  return count;
 }
};
int main(int argc,char**argv){
 try{
  require(argc==3,"usage: audit15 states.txt report.json");auto t0=chrono::steady_clock::now();ifstream f(argv[1]);require(bool(f),"open input");unsigned long long n;
  require(bool(f>>B>>I>>FN>>n),"header");require(B+I==N&&10<=B&&B<=14&&FN==28-B,"header values");require(n>0&&n<20000000,"header number");
  vector<Face> triples;for(int a=0;a<N;a++)for(int b=a+1;b<N;b++)for(int c=b+1;c<N;c++)triples.push_back({a,b,c});
  vector<Code> canonical;canonical.reserve(n);unsigned long long mass=0,group=2*B*fac(I);unsigned minsta=~0u,maxsta=0,brutes=0;
  for(unsigned long long k=0;k<n;k++){
   vector<Face> fs;fs.reserve(FN);
   for(int j=0;j<FN;j++){int id;require(bool(f>>id),"truncated state");require(0<=id&&id<(int)triples.size(),"triple index");fs.push_back(triples[id]);}
   Complex C(move(fs));auto [co,sta]=C.canonical();require(group%sta==0,"stabilizer divisor");mass+=group/sta;minsta=min(minsta,sta);maxsta=max(maxsta,sta);canonical.push_back(co);
   if(k<3||k%250000==0){require(C.brute_stab()==sta,"full-group stabilizer mismatch");brutes++;}
   if(k%500000==0)cerr<<"audited "<<k<<" / "<<n<<"\n";
  }
  string extra;require(!(f>>extra),"trailing states");
  sort(canonical.begin(),canonical.end());require(adjacent_find(canonical.begin(),canonical.end())==canonical.end(),"duplicate orbit");
  double sec=chrono::duration<double>(chrono::steady_clock::now()-t0).count();ofstream o(argv[2]);require(bool(o),"open output");
  o<<"{\n  \"verified\": true,\n  \"B\": "<<B<<",\n  \"I\": "<<I<<",\n  \"representatives\": "<<n<<",\n  \"orbit_mass\": "<<mass<<",\n  \"group_size\": "<<group<<",\n  \"minimum_stabilizer\": "<<minsta<<",\n  \"maximum_stabilizer\": "<<maxsta<<",\n  \"additional_full_group_checks\": "<<brutes<<",\n  \"seconds\": "<<sec<<"\n}\n";
  cout<<"INDEPENDENT CENSUS VERIFIED B="<<B<<" I="<<I<<" representatives="<<n<<" mass="<<mass<<" seconds="<<sec<<"\n";return 0;
 }catch(exception const&e){cerr<<"AUDIT FAILED: "<<e.what()<<"\n";return 1;}
}
