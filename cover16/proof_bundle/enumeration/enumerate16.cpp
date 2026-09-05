// Canonical enumeration of simple triangulations of a disk.
// Interior labels are eliminated by deterministic rooted dual traversal.
// Validity is preserved by diagonal flips; orbit-mass equality proves completeness.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std; using boost::multiprecision::cpp_int;
constexpr int MAXN=16,MAXF=19;
int N,B,I,F;
struct State {array<uint16_t,MAXF> t{};};
bool eq(const State&a,const State&b){return a.t==b.t;}
uint64_t hashstate(const State&s){uint64_t h=1469598103934665603ULL;for(int i=0;i<F;i++){h^=s.t[i];h*=1099511628211ULL;} h^=h>>33;h*=0xff51afd7ed558ccdULL;h^=h>>33;return h;}
struct Topo{
 uint16_t adj[MAXN]{};
 uint32_t edgefaces[MAXN*MAXN]{};
 int a[MAXF],b[MAXF],c[MAXF];
 uint64_t sig[MAXN]{};
 void init(const State&s){
  for(int f=0;f<F;f++){
   unsigned m=s.t[f]; a[f]=__builtin_ctz(m);m&=m-1;b[f]=__builtin_ctz(m);m&=m-1;c[f]=__builtin_ctz(m);
   int x[3]={a[f],b[f],c[f]};
   for(int i=0;i<3;i++){adj[x[i]]|=s.t[f]^(1u<<x[i]);for(int j=i+1;j<3;j++)edgefaces[x[i]*N+x[j]]|=1u<<f;}
  }
  for(int u=0;u<B;u++){
   unsigned im=adj[u] & (~((1u<<B)-1));int ds[5],ni=0;
   while(im){int v=__builtin_ctz(im);im&=im-1;ds[ni++]=__builtin_popcount(adj[v]);}
   sort(ds,ds+ni);uint64_t s=(__builtin_popcount(adj[u])*8+ni);
   for(int k=0;k<5;k++)s=(s<<4)|(k<ni?ds[k]:0);
   sig[u]=s;
  }
 }
 uint32_t ef(int u,int v) const{if(u>v)swap(u,v);return edgefaces[u*N+v];}
};
State rooted(const State&s,const Topo&T,int st,int dr){
 int mp[MAXN];fill(mp,mp+N,-1);for(int k=0;k<B;k++)mp[(st+dr*k+B)%B]=k;
 int a=st,b=(st+dr+B)%B;uint32_t e=T.ef(a,b);if(__builtin_popcount(e)!=1)throw runtime_error("root edge");
 int sf[3*MAXF],su[3*MAXF],sv[3*MAXF],sp=0;sf[sp]=__builtin_ctz(e);su[sp]=a;sv[sp++]=b;
 uint32_t visited=0;int nn=B,nf=0;State out;
 while(sp){--sp;int f=sf[sp],u=su[sp],v=sv[sp];if(visited&(1u<<f))continue;visited|=1u<<f;
  int w=T.a[f]^T.b[f]^T.c[f]^u^v;if(mp[w]<0)mp[w]=nn++;
  if(mp[u]<0||mp[v]<0)throw runtime_error("traversal label");out.t[nf++]=(1u<<mp[u])|(1u<<mp[v])|(1u<<mp[w]);
  int us[2]={w,v},vs[2]={u,w};
  for(int j=0;j<2;j++){uint32_t fs=T.ef(us[j],vs[j])&~(1u<<f)&~visited;if(fs){if(__builtin_popcount(fs)!=1)throw runtime_error("nonmanifold");sf[sp]=__builtin_ctz(fs);su[sp]=us[j];sv[sp++]=vs[j];}}
 }
 if(nn!=N||nf!=F)throw runtime_error("disconnected");sort(out.t.begin(),out.t.begin()+F);return out;
}
pair<State,int> canonical(const State&s){
 Topo T;T.init(s);int ks[2*MAXN],ds[2*MAXN],count=0;uint64_t best[MAXN]{};
 for(int dr:{1,-1})for(int st=0;st<B;st++){
  int cmp=count?0:-1;
  if(count)for(int j=0;j<B;j++){auto z=T.sig[(st+dr*j+B)%B];if(z<best[j]){cmp=-1;break;}if(z>best[j]){cmp=1;break;}}
  if(cmp<0){count=0;for(int j=0;j<B;j++)best[j]=T.sig[(st+dr*j+B)%B];}
  if(cmp<=0){ks[count]=st;ds[count++]=dr;}
 }
 State bests;int stab=0;
 for(int j=0;j<count;j++){auto z=rooted(s,T,ks[j],ds[j]);if(!stab||z.t<bests.t){bests=z;stab=1;}else if(eq(z,bests))stab++;}
 return {bests,stab};
}
struct Store{
 vector<State> states;vector<uint32_t> table;size_t mask;
 Store(size_t estimate){states.reserve(estimate);size_t sz=1;while(sz<estimate*1.6)sz*=2;table.assign(sz,0);mask=sz-1;}
 void rehash(){vector<uint32_t> nt(table.size()*2,0);size_t m=nt.size()-1;for(uint32_t i=0;i<states.size();i++){size_t h=hashstate(states[i])&m;while(nt[h])h=(h+1)&m;nt[h]=i+1;}table.swap(nt);mask=m;}
 bool insert(const State&s){if((states.size()+1)*10>table.size()*7)rehash();size_t h=hashstate(s)&mask;while(table[h]){if(eq(states[table[h]-1],s))return false;h=(h+1)&mask;}states.push_back(s);table[h]=states.size();return true;}
};
cpp_int fact(int n){cpp_int z=1;for(int i=2;i<=n;i++)z*=i;return z;}
uint64_t brown(){cpp_int x=2*fact(2*B-3)*fact(4*I+2*B-5),y=fact(B-1)*fact(B-3)*fact(3*I+2*B-3);if(x%y)throw runtime_error("noninteger");return (uint64_t)(x/y);}
bool validate(const State&s){
 if(!is_sorted(s.t.begin(),s.t.begin()+F))return false;Topo T;T.init(s);unsigned used=0;int ne=0;
 for(int f=0;f<F;f++){if(__builtin_popcount(s.t[f])!=3||s.t[f]>=(1u<<N)||(f&&s.t[f]==s.t[f-1]))return false;used|=s.t[f];}
 if(used!=(1u<<N)-1)return false;
 for(int u=0;u<N;u++)for(int v=u+1;v<N;v++){
  int n=__builtin_popcount(T.ef(u,v));bool bd=u<B&&v<B&&(v==u+1||(u==0&&v==B-1));if(bd&&n!=1)return false;if(n){ne++;if(n!=(bd?1:2))return false;}
 }
 if(N-ne+F!=1)return false;
 for(int v=0;v<N;v++){
  unsigned adj[MAXN]{},vs=0;
  for(int f=0;f<F;f++)if(s.t[f]>>v&1){unsigned p=s.t[f]^(1u<<v);int a=__builtin_ctz(p);p&=p-1;int b=__builtin_ctz(p);adj[a]|=1u<<b;adj[b]|=1u<<a;vs|=(1u<<a)|(1u<<b);}
  int ends=0;for(int j=0;j<N;j++)if(vs>>j&1){int d=__builtin_popcount(adj[j]);if(d==1)ends++;else if(d!=2)return false;}
  if(ends!=(v<B?2:0))return false;
  if(v<B){int p=(v+B-1)%B,q=(v+1)%B;if(__builtin_popcount(adj[p])!=1||__builtin_popcount(adj[q])!=1)return false;}
  unsigned reached=1u<<__builtin_ctz(vs),old=0;while(old!=reached){old=reached;for(int j=0;j<N;j++)if(reached>>j&1)reached|=adj[j];}if(reached!=vs)return false;
 }
 return true;
}
int main(int argc,char**argv){
 if(argc<4){cerr<<"usage: enum B I output_prefix [max_states_for_test]\n";return 2;}
 B=stoi(argv[1]);I=stoi(argv[2]);N=B+I;F=B+2*I-2;string outp=argv[3];if(B<3||I<0||I>5||N>MAXN||F>MAXF)throw runtime_error("size");
 uint64_t total=brown(),group=2*B;for(int k=2;k<=I;k++)group*=k;
 size_t estimate=total/group+total/group/50+10000;Store st(estimate);
 State seed;int nf=0;for(int j=1;j<B-1;j++)seed.t[nf++]=(1u<<0)|(1u<<j)|(1u<<(j+1));
 for(int v=B;v<N;v++){unsigned m=seed.t[0];int a=__builtin_ctz(m);m&=m-1;int b=__builtin_ctz(m);m&=m-1;int c=__builtin_ctz(m);seed.t[0]=(1u<<v)|(1u<<a)|(1u<<b);seed.t[nf++]=(1u<<v)|(1u<<b)|(1u<<c);seed.t[nf++]=(1u<<v)|(1u<<c)|(1u<<a);sort(seed.t.begin(),seed.t.begin()+nf);}
 if(nf!=F||!validate(seed))throw runtime_error("seed invalid");auto s0=canonical(seed);st.insert(s0.first);uint64_t mass=group/s0.second,attempts=0;auto t0=chrono::steady_clock::now();
 uint64_t lim=argc>4?stoull(argv[4]):0;
 for(size_t idx=0;idx<st.states.size();idx++){
  State s=st.states[idx];Topo T;T.init(s);
  for(int u=0;u<N;u++)for(int v=u+1;v<N;v++){
   auto ef=T.ef(u,v);if(__builtin_popcount(ef)!=2)continue;int f=__builtin_ctz(ef);ef&=ef-1;int g=__builtin_ctz(ef);unsigned uv=(1u<<u)|(1u<<v),cm=s.t[f]^uv,dm=s.t[g]^uv;int c=__builtin_ctz(cm),d=__builtin_ctz(dm);if(T.adj[c]>>d&1)continue;
   State ns=s;ns.t[f]=cm|dm|(1u<<u);ns.t[g]=cm|dm|(1u<<v);sort(ns.t.begin(),ns.t.begin()+F);auto z=canonical(ns);attempts++;
   if(st.insert(z.first)){if(group%z.second)throw runtime_error("stabilizer divide");mass+=group/z.second;}
  }
  if(idx%10000==0){double sec=chrono::duration<double>(chrono::steady_clock::now()-t0).count();cerr<<"processed "<<idx<<" states "<<st.states.size()<<" /~ "<<estimate<<" sec "<<sec<<"\n";}
  if(lim&&idx>=lim){cerr<<"TEST LIMIT reached, NO COMPLETE certificate\n";return 4;}
 }
 if(mass!=total)throw runtime_error("orbit mass mismatch "+to_string(mass)+" vs "+to_string(total));
 uint64_t nr=st.states.size();ofstream bin(outp+".masks",ios::binary);uint32_t hdr[5]={0x4449534b,(uint32_t)B,(uint32_t)I,(uint32_t)F,(uint32_t)nr};bin.write((char*)hdr,sizeof(hdr));
 for(auto const&s:st.states){if(!validate(s))throw runtime_error("invalid representative");bin.write((char*)s.t.data(),F*2);}bin.close();
 double sec=chrono::duration<double>(chrono::steady_clock::now()-t0).count();ofstream j(outp+".json");j<<"{\"verified\":true,\"B\":"<<B<<",\"I\":"<<I<<",\"N\":"<<N<<",\"F\":"<<F<<",\"orbits\":"<<nr<<",\"orbit_mass\":"<<mass<<",\"expected\":"<<total<<",\"group_size\":"<<group<<",\"attempts\":"<<attempts<<",\"seconds\":"<<sec<<"}\n";cout<<"COMPLETE "<<B<<" "<<I<<" orbits "<<nr<<" mass "<<mass<<" seconds "<<sec<<"\n";
}
