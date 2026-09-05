/* Independent root-face enumeration and integer angular pruning, n=18.
 * Uses a different root edge, opposite subregion order, breadth-first
 * canonical codes, and Floyd--Warshall instead of Bellman--Ford.
 * No archive of alleged excluded topologies is trusted.
 */
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
constexpr int MAXV=18,MAXT=23;
constexpr long long UNIT=10000000, ARC=937111, TWO=1970786, THREE=3362065;
int n,b,interior,faces_required;
bool do_screen=true;
uint64_t calls=0,negative=0,completed=0,forbidden=0,metric_checks=0,full_group_tests=0;
struct Region{vector<int> rim;int budget;};
struct Board{array<uint32_t,MAXV> neighbors{};vector<uint32_t> triangles;int fresh;};
struct Orbit{uint64_t visits=0;int stabilizer=0;};
map<vector<uint32_t>,Orbit> output;
void require(bool v,const string&why){if(!v)throw runtime_error(why);}

bool angularly_possible(Board const&s){
 if(!do_screen)return true;
 metric_checks++;
 const long long INF=(1LL<<50);
 long long d[MAXV][MAXV];
 for(int i=0;i<b;i++)for(int j=0;j<b;j++)d[i][j]=(i==j?0:INF);
 for(int i=0;i<b;i++){
  int j=(i+1)%b;
  d[i][j]=min(d[i][j],ARC-(j==0?UNIT:0));
  d[j][i]=min(d[j][i],j==0?UNIT:0);
 }
 for(int i=0;i<b;i++){
  uint32_t reach[3];reach[0]=(1u<<i)|(1u<<((i+1)%b));
  for(int h=1;h<=2;h++){
   reach[h]=reach[h-1];
   for(int v=0;v<s.fresh;v++)if(reach[h-1]>>v&1)reach[h]|=s.neighbors[v];
  }
  for(int j=i+1;j<b;j++){
   uint32_t terminal=(1u<<j)|(1u<<((j+1)%b));int hop=0;
   while(hop<=2&&!(reach[hop]&terminal))hop++;
   if(hop==3)continue;
   long long cap=hop==0?ARC:(hop==1?TWO:THREE);
   int k=j-i;
   if(min(k,b-k)*ARC<=cap)continue;
   bool forward=k*ARC<UNIT-cap,backward=(b-k)*ARC<UNIT-cap;
   if(forward&&backward)return false;
   if(forward)d[i][j]=min(d[i][j],cap);
   if(backward)d[j][i]=min(d[j][i],cap-UNIT);
  }
 }
 for(int k=0;k<b;k++){
  for(int i=0;i<b;i++)if(d[i][k]!=INF)for(int j=0;j<b;j++)if(d[k][j]!=INF)d[i][j]=min(d[i][j],d[i][k]+d[k][j]);
  for(int i=0;i<b;i++)if(d[i][i]<0)return false;
 }
 return true;
}

struct Complex {
 vector<array<int,3>> tri;
 uint32_t inc[MAXV][MAXV]{};
 array<uint32_t,MAXV> adj{};
 explicit Complex(Board const&s){
  require(s.fresh==n && (int)s.triangles.size()==faces_required,"wrong complete size");
  set<uint32_t> unique;
  uint32_t used=0;
  for(uint32_t mask:s.triangles){
   require(__builtin_popcount(mask)==3 && mask<(1u<<n),"invalid triangle mask");
   require(unique.insert(mask).second,"repeated face");used|=mask;
   array<int,3> t{};uint32_t rest=mask;
   for(int k=0;k<3;k++){t[k]=__builtin_ctz(rest);rest&=rest-1;}
   int f=tri.size();tri.push_back(t);
   for(int j=0;j<3;j++)for(int k=j+1;k<3;k++){
    int u=t[j],v=t[k];inc[u][v]|=1u<<f;inc[v][u]|=1u<<f;adj[u]|=1u<<v;adj[v]|=1u<<u;
   }
  }
  require(used==(1u<<n)-1,"unused vertex");
  require(adj==s.neighbors,"edge/face disagreement");
  int edges=0;
  for(int u=0;u<n;u++)for(int v=u+1;v<n;v++){
   int c=__builtin_popcount(inc[u][v]);bool boundary=u<b&&v<b&&(v==u+1||(u==0&&v==b-1));
   require(c==(boundary?1:(c?2:0)),"edge multiplicity");if(c)edges++;
  }
  require(n-edges+faces_required==1,"Euler characteristic");
  uint32_t seen=1,old=0;
  while(seen!=old){old=seen;for(int f=0;f<faces_required;f++)if(seen>>f&1)for(int j=0;j<3;j++)seen|=inc[tri[f][j]][tri[f][(j+1)%3]];}
  require(seen==(1u<<faces_required)-1,"disconnected dual");
  for(int v=0;v<n;v++){
   uint32_t link[MAXV]{},vertices=0;
   for(auto t:tri)if(t[0]==v||t[1]==v||t[2]==v){
    int z[2],at=0;for(int u:t)if(u!=v)z[at++]=u;
    link[z[0]]|=1u<<z[1];link[z[1]]|=1u<<z[0];vertices|=(1u<<z[0])|(1u<<z[1]);
   }
   uint32_t ends=0;
   for(int u=0;u<n;u++)if(vertices>>u&1){int deg=__builtin_popcount(link[u]);require(deg==1||deg==2,"link degree");if(deg==1)ends|=1u<<u;}
   uint32_t expected=v<b?((1u<<((v+b-1)%b))|(1u<<((v+1)%b))):0;
   require(ends==expected,"link ends");
   uint32_t reach=1u<<__builtin_ctz(vertices),before=0;
   while(reach!=before){before=reach;for(int u=0;u<n;u++)if(reach>>u&1)reach|=link[u];}
   require(reach==vertices,"disconnected link");
  }
 }
 vector<uint32_t> rooted_code(int origin,int direction)const{
  int label[MAXV];fill(label,label+MAXV,-1);
  for(int k=0;k<b;k++)label[(origin+direction*k+b)%b]=k;
  struct Entry{int f,u,v;};
  int u=origin,v=(origin+direction+b)%b;uint32_t r=inc[u][v];require(__builtin_popcount(r)==1,"root incidence");
  vector<Entry> queue{{__builtin_ctz(r),u,v}};uint32_t scheduled=r;int next=b;vector<uint32_t> code;
  for(size_t head=0;head<queue.size();head++){
   auto e=queue[head];int w=tri[e.f][0]^tri[e.f][1]^tri[e.f][2]^e.u^e.v;
   if(label[w]<0)label[w]=next++;
   require(label[e.u]>=0&&label[e.v]>=0,"BFS labels");
   code.push_back((1u<<label[e.u])|(1u<<label[e.v])|(1u<<label[w]));
   const int es[2][2]={{e.v,w},{w,e.u}};
   for(auto const&a:es){uint32_t nb=inc[a[0]][a[1]]&~scheduled;if(nb){require(__builtin_popcount(nb)==1,"BFS incidence");scheduled|=nb;queue.push_back({__builtin_ctz(nb),a[0],a[1]});}}
  }
  require(next==n&&(int)code.size()==faces_required,"BFS incomplete");sort(code.begin(),code.end());return code;
 }
 pair<vector<uint32_t>,int> normalize()const{
  vector<uint32_t> best;int stab=0;
  // No profile filter: inspect every oriented boundary root.
  for(int direction:{1,-1})for(int origin=0;origin<b;origin++){
   auto z=rooted_code(origin,direction);
   if(best.empty()||z<best){best=z;stab=1;}else if(best==z)stab++;
  }
  return {best,stab};
 }
 int permutation_stabilizer()const{
  vector<int> p;for(int v=b;v<n;v++)p.push_back(v);
  vector<uint32_t> base;for(auto t:tri)base.push_back((1u<<t[0])|(1u<<t[1])|(1u<<t[2]));sort(base.begin(),base.end());int count=0;
  do{
   for(int sg:{1,-1})for(int sh=0;sh<b;sh++){
    int mp[MAXV];for(int i=0;i<b;i++)mp[i]=(sh+sg*i+b)%b;for(int i=b;i<n;i++)mp[i]=p[i-b];
    vector<uint32_t> image;for(auto t:tri)image.push_back((1u<<mp[t[0]])|(1u<<mp[t[1]])|(1u<<mp[t[2]]));sort(image.begin(),image.end());
    if(image==base)count++;
   }
  }while(next_permutation(p.begin(),p.end()));
  return count;
 }
};

bool boundary_edge(vector<int> const&p,int u,int v){
 for(size_t j=0;j<p.size();j++)if((p[j]==u&&p[(j+1)%p.size()]==v)||(p[j]==v&&p[(j+1)%p.size()]==u))return true;
 return false;
}
bool add_triangle(Board&s,vector<int>const&poly,int a,int c,int d){
 const int e[3][2]={{a,c},{c,d},{d,a}};
 for(auto const&edge:e){int u=edge[0],v=edge[1];if((s.neighbors[u]>>v&1)&&!boundary_edge(poly,u,v))return false;}
 for(auto const&edge:e){int u=edge[0],v=edge[1];s.neighbors[u]|=1u<<v;s.neighbors[v]|=1u<<u;}
 s.triangles.push_back((1u<<a)|(1u<<c)|(1u<<d));return true;
}

void enumerate(Board const&s,vector<Region> const&tasks){
 calls++;
 int pending_interior=0,remaining_faces=0;
 for(auto const&r:tasks){
  require(r.rim.size()>=2&&r.budget>=0,"bad region");
  require(set<int>(r.rim.begin(),r.rim.end()).size()==r.rim.size(),"repeated rim vertex");
  for(size_t i=0;i<r.rim.size();i++){
   int u=r.rim[i],v=r.rim[(i+1)%r.rim.size()];
   require(u>=0&&u<s.fresh&&(s.neighbors[u]>>v&1),"missing rim edge");
  }
  pending_interior+=r.budget;remaining_faces+=2*r.budget+(int)r.rim.size()-2;
 }
 require(s.fresh+pending_interior==n,"lost interior budget");
 require((int)s.triangles.size()+remaining_faces==faces_required,"lost face budget");
 if(tasks.empty()){
  completed++;Complex c(s);auto [code,stab]=c.normalize();
  auto &o=output[code];if(o.visits==0){o.stabilizer=stab;
   // Explicit full group audit at the first unique complete map in each family.
   if(output.size()==1){require(c.permutation_stabilizer()==stab,"full group disagreement");full_group_tests++;}
  }else require(o.stabilizer==stab,"stabilizer disagreement");
  o.visits++;return;
 }
 auto current=tasks.front();vector<Region> tail(tasks.begin()+1,tasks.end());
 int k=current.budget,L=current.rim.size();
 if(L==2){if(k==0)enumerate(s,tail);else forbidden++;return;}
 vector<int> p{current.rim.back()};p.insert(p.end(),current.rim.begin(),current.rim.end()-1);
 int a=p[0],v=p[1];
 // Existing third boundary vertex; all allocations, opposite subregion order.
 for(int index=L-1;index>=2;index--){
  int w=p[index];Board ns=s;
  if(!add_triangle(ns,p,a,v,w)){forbidden++;continue;}
  if(!angularly_possible(ns)){negative++;continue;}
  vector<int> left{w};left.insert(left.end(),p.begin()+1,p.begin()+index);
  vector<int> right{a};right.insert(right.end(),p.begin()+index,p.end());
  for(int j=k;j>=0;j--){
   if((left.size()==2&&j>0)||(right.size()==2&&k-j>0))continue;
   vector<Region> work{{right,k-j},{left,j}};work.insert(work.end(),tail.begin(),tail.end());enumerate(ns,work);
  }
 }
 // A new vertex in the current region. No previously exposed vertex is interior.
 if(k>0){
  Board ns=s;int w=ns.fresh++;require(w<n,"vertex budget exceeded");
  require(add_triangle(ns,p,a,v,w),"fresh vertex edge already present");
  if(!angularly_possible(ns)){negative++;return;}
  vector<int> rim{a,w};rim.insert(rim.end(),p.begin()+1,p.end());
  vector<Region> work{{rim,k-1}};work.insert(work.end(),tail.begin(),tail.end());enumerate(ns,work);
 }
}
int main(int argc,char**argv){try{
 require(argc==4||argc==5,"usage: audit_peeling N B output [no-screen]");
 n=stoi(argv[1]);b=stoi(argv[2]);interior=n-b;faces_required=2*n-b-2;
 require(3<=b&&b<=n&&n<=MAXV&&faces_required<=MAXT,"dimensions");
 if(argc==5){require(string(argv[4])=="no-screen","unknown flag");do_screen=false;}
 Board s;s.fresh=b;vector<int> rim;for(int i=0;i<b;i++){rim.push_back(i);int j=(i+1)%b;s.neighbors[i]|=1u<<j;s.neighbors[j]|=1u<<i;}
 auto t=chrono::steady_clock::now();
 if(angularly_possible(s))enumerate(s,{{rim,interior}});else negative++;
 for(auto const&[code,o]:output){require(o.stabilizer>0&&(2*b)%o.stabilizer==0,"stabilizer divisor");require(o.visits==(uint64_t)(2*b/o.stabilizer),"rooted multiplicity mismatch");}
 ofstream f(argv[3]);require((bool)f,"output open");
 f<<"{\"verified\":true,\"N\":"<<n<<",\"B\":"<<b<<",\"I\":"<<interior<<",\"angle_scale\":"<<UNIT<<",\"angle_caps\":["<<ARC<<","<<TWO<<","<<THREE<<"],\"screen\":"<<(do_screen?"true":"false")<<",\"calls\":"<<calls<<",\"metric_rejections\":"<<negative<<",\"metric_checks\":"<<metric_checks<<",\"rooted_survivors\":"<<completed<<",\"surviving_orbits\":"<<output.size()<<",\"full_group_crosschecks\":"<<full_group_tests<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<",\"orbits\":[";
 bool first=true;for(auto const&[code,o]:output){if(!first)f<<",";first=false;f<<"{\"masks\":[";for(size_t i=0;i<code.size();i++){if(i)f<<",";f<<code[i];}f<<"],\"stabilizer\":"<<o.stabilizer<<",\"rooted_multiplicity\":"<<o.visits<<"}";}f<<"]}\n";f.close();require((bool)f,"output write");
 cout<<"INDEPENDENT PEELING VERIFIED N="<<n<<" B="<<b<<" calls="<<calls<<" rooted="<<completed<<" orbits="<<output.size()<<" seconds="<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<"\n";
 return 0;
}catch(exception const&e){cerr<<"PEELING REJECTED: "<<e.what()<<"\n";return 1;}}
