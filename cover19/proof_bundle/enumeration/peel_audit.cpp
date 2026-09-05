// Independent queue-based peeling audit. Regions are processed FIFO, the root
// edge is rotated, and interior labels run downwards. Angle closure is rebuilt
// from scratch by integer Floyd-Warshall; canonical labeling is breadth-first.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
int n,b,aa,c4,c6; const long long den=100000;
bool pruning=true; uint64_t visits=0,cuts=0,leaves=0;
struct Region {vector<int> boundary; int interior;};
using Face=array<int,3>;using Faces=vector<Face>;using Code=vector<uint32_t>;
map<Code,uint64_t> generated;map<Code,int> stabilizers;
array<uint32_t,20> adjacency(const Faces&fs){
 array<uint32_t,20> a{};for(int i=0;i<n;i++)a[i]=1u<<i;
 for(int i=0;i<b;i++){a[i]|=1u<<((i+1)%b);a[(i+1)%b]|=1u<<i;}
 for(auto f:fs){uint32_t m=(1u<<f[0])|(1u<<f[1])|(1u<<f[2]);for(int v:f)a[v]|=m;}
 return a;
}
bool possible(const array<uint32_t,20>&a){
 if(!pruning)return true;
 long long d[20][20];for(int i=0;i<b;i++)for(int j=0;j<b;j++)d[i][j]=i==j?0:10000000;
 for(int i=0;i<b-1;i++){d[i][i+1]=aa;d[i+1][i]=0;}d[b-1][0]=aa-den;d[0][b-1]=den;
 uint32_t reach[20];for(int v=0;v<n;v++){reach[v]=0;for(int w=0;w<n;w++)if(a[v]>>w&1)reach[v]|=a[w];}
 for(int i=0;i<b;i++)for(int j=i+1;j<b;j++){
  uint32_t ends=(1u<<j)|(1u<<((j+1)%b));long long cap;
  if((a[i]|a[(i+1)%b])&ends)cap=c4;
  else if((reach[i]|reach[(i+1)%b])&ends)cap=c6;
  else continue;
  bool forward=(j-i)*aa<den-cap,backward=(b-j+i)*aa<den-cap;
  if(forward&&backward)return false;
  if(forward)d[i][j]=min(d[i][j],cap);
  if(backward)d[j][i]=min(d[j][i],cap-den);
 }
 for(int k=0;k<b;k++){
  for(int i=0;i<b;i++)for(int j=0;j<b;j++)d[i][j]=min(d[i][j],d[i][k]+d[k][j]);
  for(int i=0;i<b;i++)if(d[i][i]<0)return false;
 }
 return true;
}
void validate(const Faces&fs){
 if(int(fs.size())!=2*n-b-2)throw runtime_error("face count");
 int incidence[20][20]{};uint32_t used=0;set<Face> distinct;
 for(auto f:fs){sort(f.begin(),f.end());if(f[0]<0||f[2]>=n||f[0]==f[1]||f[1]==f[2]||!distinct.insert(f).second)throw runtime_error("invalid face");
  for(int v:f)used|=1u<<v;for(int j=0;j<3;j++)for(int k=j+1;k<3;k++)incidence[f[j]][f[k]]++;
 }
 if(used!=(1u<<n)-1)throw runtime_error("unused vertex");int ne=0;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
  bool bound=j<b&&(j==i+1||(i==0&&j==b-1));int count=incidence[i][j];
  if(bound&&count!=1)throw runtime_error("boundary");if(count){ne++;if(count!=(bound?1:2))throw runtime_error("incidence");}
 }
 if(n-ne+int(fs.size())!=1)throw runtime_error("Euler characteristic");
 for(int v=0;v<n;v++){
  uint32_t link[20]{};uint32_t lv=0;for(auto f:fs)if(find(f.begin(),f.end(),v)!=f.end()){
   vector<int> q;for(int x:f)if(x!=v)q.push_back(x);link[q[0]]|=1u<<q[1];link[q[1]]|=1u<<q[0];lv|=(1u<<q[0])|(1u<<q[1]);
  }
  vector<int> ends;for(int i=0;i<n;i++)if(lv>>i&1){int deg=__builtin_popcount(link[i]);if(deg==1)ends.push_back(i);else if(deg!=2)throw runtime_error("link degree");}
  vector<int> expected;if(v<b){expected={(v+b-1)%b,(v+1)%b};sort(expected.begin(),expected.end());}if(ends!=expected)throw runtime_error("link endpoints");
  uint32_t seen=1u<<__builtin_ctz(lv),old=0;while(old!=seen){old=seen;for(int i=0;i<n;i++)if(seen>>i&1)seen|=link[i];}if(seen!=lv)throw runtime_error("link disconnected");
 }
}
pair<Code,int> canonical(const Faces&fs){
 int ef[20][20][2];int cnt[20][20]{};
 for(int k=0;k<int(fs.size());k++)for(int j=0;j<3;j++)for(int h=j+1;h<3;h++){
  int u=min(fs[k][j],fs[k][h]),v=max(fs[k][j],fs[k][h]);ef[u][v][cnt[u][v]++]=k;
 }
 Code best;int st=0;
 for(int direction:{-1,1})for(int origin=0;origin<b;origin++){
  int label[20];fill(label,label+n,-1);for(int j=0;j<b;j++)label[(origin+direction*j+b)%b]=j;int next=b;
  int u=origin,v=(origin+direction+b)%b;vector<array<int,3>> q{{ef[min(u,v)][max(u,v)][0],u,v}};uint64_t seen=0;Code code;
  for(size_t it=0;it<q.size();it++){
   int f=q[it][0],u=q[it][1],v=q[it][2];if(seen>>f&1)continue;seen|=1ULL<<f;
   int w=-1;for(int x:fs[f])if(x!=u&&x!=v)w=x;if(w<0)throw runtime_error("missing third vertex");
   if(label[w]<0)label[w]=next++;if(label[u]<0||label[v]<0)throw runtime_error("BFS dependency");
   code.push_back((1u<<label[u])|(1u<<label[v])|(1u<<label[w]));
   int us[2]={v,w},vs[2]={w,u};for(int j=0;j<2;j++){
    int x=min(us[j],vs[j]),y=max(us[j],vs[j]);for(int k=0;k<cnt[x][y];k++)if(ef[x][y][k]!=f)q.push_back({ef[x][y][k],us[j],vs[j]});
   }
  }
  if(next!=n||code.size()!=fs.size())throw runtime_error("canonical connectivity");sort(code.begin(),code.end());
  if(best.empty()||code<best){best=code;st=1;}else if(code==best)st++;
 }
 return {best,st};
}
void enumerate(const vector<Region>&regions,int remaining_label,const Faces&faces){
 visits++;auto a=adjacency(faces);if(!possible(a)){cuts++;return;}
 if(regions.empty()){
  if(remaining_label!=b-1)throw runtime_error("unused allotted vertices");validate(faces);auto [code,stab]=canonical(faces);generated[code]++;stabilizers[code]=stab;leaves++;return;
 }
 Region current=regions.front();vector<Region> later(regions.begin()+1,regions.end());
 auto p=current.boundary;int size=p.size(),inside=current.interior;if(size<3)throw runtime_error("nonempty digon");
 // Rotate to use the last edge, rather than the first edge of the discovery implementation.
 rotate(p.begin(),p.end()-1,p.end());int u=p[0],v=p[1];
 for(int j=2;j<size;j++){
  int w=p[j];if(j!=size-1&&(a[u]>>w&1))continue;if(j!=2&&(a[v]>>w&1))continue;
  Faces nextfaces=faces;nextfaces.push_back({u,v,w});
  vector<int> left(p.begin()+1,p.begin()+j+1),right{u};right.insert(right.end(),p.begin()+j,p.end());
  for(int count=inside;count>=0;count--){
   if((left.size()==2&&count)||(right.size()==2&&inside-count))continue;
   auto pending=later;
   if(right.size()>2)pending.push_back({right,inside-count});
   if(left.size()>2)pending.push_back({left,count});
   enumerate(pending,remaining_label,nextfaces);
  }
 }
 if(inside){
  int w=remaining_label;if(w<b)throw runtime_error("interior budget");Faces nextfaces=faces;nextfaces.push_back({u,v,w});
  vector<int> boundary{u,w};boundary.insert(boundary.end(),p.begin()+1,p.end());auto pending=later;pending.push_back({boundary,inside-1});enumerate(pending,w-1,nextfaces);
 }
}
int main(int argc,char**argv){
 if(argc<9){cerr<<"usage N B A C4 C6 expected.txt report.json prune\n";return 2;}
 std::remove(argv[7]);
 n=stoi(argv[1]);b=stoi(argv[2]);aa=stoi(argv[3]);c4=stoi(argv[4]);c6=stoi(argv[5]);pruning=stoi(argv[8]);
 if(n>20||b<3||b>n||aa<=0||aa>=den/2||c4<=0||c4>=den/2||c6<c4||c6>=den/2)throw runtime_error("dimensions or constants");std::remove(argv[7]);auto start=chrono::steady_clock::now();
 map<Code,int> expected;ifstream in(argv[6]);if(!in)throw runtime_error("missing manifest");string line;
 while(getline(in,line)){
  if(line.empty())throw runtime_error("empty graph record");istringstream ss(line);int stab;if(!(ss>>stab)||stab<=0)throw runtime_error("missing stabilizer");uint32_t mask;Faces fs;string token;
  while(ss>>token){if(token.empty()||token.find_first_not_of("0123456789")!=string::npos)throw runtime_error("invalid integer token");unsigned long long raw=stoull(token);if(raw>=(1ULL<<n))throw runtime_error("mask range");mask=uint32_t(raw);Face f;int k=0;for(int j=0;j<n;j++)if(mask>>j&1){if(k==3)throw runtime_error("not triangle");f[k++]=j;}if(k!=3||mask>=(1u<<n))throw runtime_error("mask");fs.push_back(f);}
  if(!ss.eof())throw runtime_error("trailing noninteger input");
  validate(fs);if(pruning&&!possible(adjacency(fs)))throw runtime_error("infeasible expected graph");auto [code,s]=canonical(fs);if(stab!=s||!expected.emplace(code,s).second)throw runtime_error("duplicate or bad stabilizer");
 }
 if(in.bad())throw runtime_error("manifest read failed");
 vector<int> boundary{0};for(int j=b-1;j>0;j--)boundary.push_back(j);enumerate({{boundary,n-b}},n-1,{});
 if(expected.size()!=generated.size())throw runtime_error("case count mismatch "+to_string(expected.size())+" vs "+to_string(generated.size()));
 for(auto const&[code,s]:expected){auto it=generated.find(code);if(it==generated.end()||it->second*s!=uint64_t(2*b)||stabilizers[code]!=s)throw runtime_error("missing case or multiplicity");}
 ofstream out(argv[7]);if(!out)throw runtime_error("cannot write report");out<<"{\"verified\":true,\"N\":"<<n<<",\"B\":"<<b<<",\"A\":"<<aa<<",\"C4\":"<<c4<<",\"C6\":"<<c6<<",\"visits\":"<<visits<<",\"pruned\":"<<cuts<<",\"rooted_survivors\":"<<leaves<<",\"surviving_orbits\":"<<generated.size()<<",\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";
 out.close();if(!out)throw runtime_error("report write failed");
 cout<<"INDEPENDENT PEEL AUDIT VERIFIED "<<n<<" "<<b<<" "<<generated.size()<<" orbits visits "<<visits<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"\n";
}
