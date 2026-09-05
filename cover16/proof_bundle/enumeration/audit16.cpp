// Independent breadth-first rooted normalization and exact metric partition.
// Only the fixed-width storage/hash container is shared with discovery.
#include "state_store.hpp"
#include <map>
#include <sstream>
#include <unordered_map>
struct CheckedGraph {
 array<array<int,3>,MAXF> tri{};
 array<uint16_t,MAXN> adj{};
 array<uint32_t,MAXN*MAXN> inc{};
 uint32_t edge(int u,int v)const{if(u>v)swap(u,v);return inc[u*MAXN+v];}
 explicit CheckedGraph(State const&s){
  if(!is_sorted(s.t.begin(),s.t.begin()+F))throw runtime_error("unsorted faces");unsigned used=0;
  for(int f=0;f<F;f++){
   unsigned mask=s.t[f];if(__builtin_popcount(mask)!=3||mask>=(1u<<N)||(f&&mask==s.t[f-1]))throw runtime_error("bad triangle");used|=mask;
   for(int k=0;k<3;k++){tri[f][k]=__builtin_ctz(mask);mask&=mask-1;}
   for(int j=0;j<3;j++)for(int k=j+1;k<3;k++){int u=tri[f][j],v=tri[f][k];inc[u*MAXN+v]|=1u<<f;adj[u]|=1u<<v;adj[v]|=1u<<u;}
  }
  if(used!=(1u<<N)-1)throw runtime_error("unused vertex");int E=0;
  for(int u=0;u<N;u++)for(int v=u+1;v<N;v++){
   int n=__builtin_popcount(edge(u,v));bool bd=u<B&&v<B&&(v==u+1||(u==0&&v==B-1));
   if(bd){if(n!=1)throw runtime_error("boundary incidence");}
   else if(n&&n!=2)throw runtime_error("interior incidence");if(n)E++;
  }
  if(N-E+F!=1)throw runtime_error("Euler");
  uint32_t seen=1,old=0;while(old!=seen){old=seen;for(int f=0;f<F;f++)if(seen>>f&1)for(int j=0;j<3;j++)seen|=edge(tri[f][j],tri[f][(j+1)%3]);}
  if(seen!=(1u<<F)-1)throw runtime_error("disconnected dual");
  for(int v=0;v<N;v++){
   uint16_t la[MAXN]{};unsigned vv=0;
   for(int f=0;f<F;f++)if(s.t[f]>>v&1){int t[2],k=0;for(int u:tri[f])if(u!=v)t[k++]=u;int a=t[0],b=t[1];la[a]|=1u<<b;la[b]|=1u<<a;vv|=(1u<<a)|(1u<<b);}
   unsigned ends=0;for(int u=0;u<N;u++)if(vv>>u&1){int d=__builtin_popcount(la[u]);if(d==1)ends|=1u<<u;else if(d!=2)throw runtime_error("bad link degree");}
   unsigned need=v<B?((1u<<((v+B-1)%B))|(1u<<((v+1)%B))):0;if(ends!=need)throw runtime_error("bad link endpoints");
   unsigned reached=1u<<__builtin_ctz(vv),prior=0;while(reached!=prior){prior=reached;for(int u=0;u<N;u++)if(reached>>u&1)reached|=la[u];}if(reached!=vv)throw runtime_error("disconnected link");
  }
 }
 State code(int start,int dir)const {
  int label[MAXN];fill(label,label+N,-1);for(int j=0;j<B;j++)label[(start+dir*j+B)%B]=j;
  int qf[MAXF],qu[MAXF],qv[MAXF],front=0,back=1;int u=start,v=(start+dir+B)%B;uint32_t root=edge(u,v);if(__builtin_popcount(root)!=1)throw runtime_error("root");
  qf[0]=__builtin_ctz(root);qu[0]=u;qv[0]=v;uint32_t scheduled=root;int next=B;State out;int nf=0;
  while(front<back){int f=qf[front],a=qu[front],b=qv[front];front++;int c=tri[f][0]^tri[f][1]^tri[f][2]^a^b;
   if(label[c]<0)label[c]=next++;if(label[a]<0||label[b]<0)throw runtime_error("BFS label");
   out.t[nf++]=(1u<<label[a])|(1u<<label[b])|(1u<<label[c]);
   int es[2][2]={{b,c},{c,a}};
   for(auto const&e:es){uint32_t nb=edge(e[0],e[1])&~scheduled;if(nb){if(__builtin_popcount(nb)!=1||back>=MAXF)throw runtime_error("BFS incidence");scheduled|=nb;qf[back]=__builtin_ctz(nb);qu[back]=e[0];qv[back]=e[1];back++;}}
  }
  if(nf!=F||next!=N)throw runtime_error("BFS incomplete");sort(out.t.begin(),out.t.begin()+F);return out;
 }
 pair<State,int> normalize()const {
  array<int,MAXN> profile{};bool have=false;vector<pair<int,int>> roots;
  for(int dir:{1,-1})for(int start=0;start<B;start++){
   array<int,MAXN> cur{};for(int k=0;k<B;k++)cur[k]=__builtin_popcount(adj[(start+dir*k+B)%B]);
   if(!have||cur<profile){profile=cur;roots.clear();have=true;}
   if(cur==profile)roots.emplace_back(start,dir);
  }
  State best;int multiplicity=0;
  for(auto [start,dir]:roots){State c=code(start,dir);if(!multiplicity||c.t<best.t){best=c;multiplicity=1;}else if(c.t==best.t)multiplicity++;}
  return {best,multiplicity};
 }
 int full_stabilizer(State const&s)const {
  vector<int> inner;for(int j=B;j<N;j++)inner.push_back(j);int count=0;
  do {for(int sg:{1,-1})for(int sh=0;sh<B;sh++){
   int mp[MAXN];for(int j=0;j<B;j++)mp[j]=(sh+sg*j+B)%B;for(int j=B;j<N;j++)mp[j]=inner[j-B];
   State image;for(int f=0;f<F;f++){unsigned m=0;for(int u:tri[f])m|=1u<<mp[u];image.t[f]=m;}
   sort(image.t.begin(),image.t.begin()+F);if(eq(image,s))count++;
  }}while(next_permutation(inner.begin(),inner.end()));
  return count;
 }
 bool signature(vector<int>&masks,int A,int C,int D)const {
  constexpr int C6=37566;
  masks.clear();
  for(int i=0;i<B;i++){
   unsigned reach=adj[i]|adj[(i+1)%B]|(1u<<i)|(1u<<((i+1)%B)), reach2=reach;
   for(int u=0;u<N;u++)if(reach>>u&1)reach2|=adj[u];
   for(int j=i+1;j<B;j++){
    unsigned ends=(1u<<j)|(1u<<((j+1)%B));int cap,tag;
    if(reach&ends){cap=C;tag=0;}else if(reach2&ends){cap=C6;tag=1u<<B;}else continue;
    int k=j-i;if(min(k,B-k)*A<=cap)continue;
    bool f=k*A<D-cap,b=(B-k)*A<D-cap;if(f&&b)return false;
    int mask=((1u<<j)-1)^((1u<<i)-1);
    if(f)masks.push_back(tag|mask);else if(b)masks.push_back(tag|(((1u<<B)-1)^mask));
   }
  }
  sort(masks.begin(),masks.end());masks.erase(unique(masks.begin(),masks.end()),masks.end());return true;
 }
};
int main(int argc,char**argv){try{
 if(argc!=5)throw runtime_error("usage: audit16 masks farkas expected_labelled output_prefix");
 ifstream in(argv[1],ios::binary);uint32_t hdr[5];if(!in.read((char*)hdr,sizeof hdr))throw runtime_error("header");
 if(hdr[0]!=0x4449534b)throw runtime_error("magic");if(hdr[1]<11||hdr[1]>15||hdr[2]!=16-hdr[1]||hdr[3]!=30-hdr[1]||!hdr[4]||hdr[4]>50000000)throw runtime_error("dimensions");
 B=(int)hdr[1];I=(int)hdr[2];F=(int)hdr[3];N=B+I;uint32_t nr=hdr[4];
 ifstream cf(argv[2]);int AA,CC,DD,K;cpp_int SS;if(!(cf>>AA>>CC>>DD>>SS>>K)||AA!=9974||CC!=21143||DD!=100000||SS!=cpp_int(1000000000000LL)||K<0||K>200000)throw runtime_error("Farkas header");
 map<vector<int>,bool> cert;
 for(int k=0;k<K;k++){
  int q;if(!(cf>>q)||q<0||q>B*B)throw runtime_error("signature size");vector<int> masks(q);
  for(int &m:masks)if(!(cf>>m)||m<=0||m>=(2<<B)||!(m&((1<<B)-1)))throw runtime_error("signature mask");
  if(!is_sorted(masks.begin(),masks.end())||adjacent_find(masks.begin(),masks.end())!=masks.end())throw runtime_error("signature order");
  vector<cpp_int> vals(B+q);for(auto &v:vals)if(!(cf>>v)||v<0)throw runtime_error("Farkas coefficient");
  for(int i=0;i<B;i++){cpp_int cover=vals[i];for(int j=0;j<q;j++)if(masks[j]>>i&1)cover+=vals[B+j];if(cover<SS)throw runtime_error("Farkas coverage");}
  cpp_int cost=0;for(int i=0;i<B+q;i++)cost+=vals[i]*(i<B?AA:((masks[i-B]>>B)?37566:CC));if(cost>=DD*SS)throw runtime_error("Farkas strict margin");
  if(!cert.emplace(masks,false).second)throw runtime_error("duplicate certificate");
 }
 string extra;if(cf>>extra)throw runtime_error("extra certificate data");
 uint64_t group=2*B;for(int j=2;j<=I;j++)group*=j;cpp_int target(argv[3]),mass=0;Store store(nr);uint32_t ndirect=0,nfark=0,nres=0;int minstab=100000,maxstab=0,nfull=0;vector<int> sig;
 string prefix=argv[4];ofstream out(prefix+".residuals",ios::binary);uint32_t oh[5]={0x52534544,(uint32_t)B,(uint32_t)I,(uint32_t)F,0};out.write((char*)oh,sizeof oh);auto st=chrono::steady_clock::now();
 for(uint32_t idx=0;idx<nr;idx++){
  State s;if(!in.read((char*)s.t.data(),F*2))throw runtime_error("truncated archive");CheckedGraph g(s);auto [code,stab]=g.normalize();if(stab<=0||group%stab)throw runtime_error("stabilizer");
  if(!store.insert(code))throw runtime_error("duplicate orbit");if(idx==0||idx==nr/2||idx+1==nr){if(g.full_stabilizer(s)!=stab)throw runtime_error("full group stabilizer");nfull++;}mass+=group/stab;minstab=min(minstab,stab);maxstab=max(maxstab,stab);
  if(!g.signature(sig,AA,CC,DD))ndirect++;
  else{auto it=cert.find(sig);if(it!=cert.end()){nfark++;it->second=true;}else{nres++;out.write((char*)&idx,sizeof idx);out.write((char*)s.t.data(),F*2);}}
 }
 char ch;if(in.read(&ch,1))throw runtime_error("extra archive bytes");if(mass!=target)throw runtime_error("orbit mass");
 for(auto const&kv:cert)if(!kv.second)throw runtime_error("unused Farkas certificate");
 out.seekp(16);out.write((char*)&nres,sizeof nres);out.close();double sec=chrono::duration<double>(chrono::steady_clock::now()-st).count();
 ofstream report(prefix+".json");report<<"{\"verified\":true,\"B\":"<<B<<",\"I\":"<<I<<",\"representatives\":"<<nr<<",\"orbit_mass\":\""<<mass<<"\",\"group_size\":"<<group<<",\"min_stabilizer\":"<<minstab<<",\"max_stabilizer\":"<<maxstab<<",\"direct\":"<<ndirect<<",\"farkas\":"<<nfark<<",\"residual\":"<<nres<<",\"full_group_crosschecks\":"<<nfull<<",\"seconds\":"<<sec<<"}\n";
 cout<<"AUDITED B="<<B<<" I="<<I<<" orbits="<<nr<<" mass="<<mass<<" direct="<<ndirect<<" farkas="<<nfark<<" residual="<<nres<<" seconds="<<sec<<"\n";return 0;
 }catch(exception const&e){cerr<<"AUDIT REJECTED: "<<e.what()<<"\n";return 1;}}
