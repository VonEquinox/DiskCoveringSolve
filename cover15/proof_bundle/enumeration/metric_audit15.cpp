// Exact independent metric/Farkas partition. Candidate generation is not trusted.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;using boost::multiprecision::cpp_int;
constexpr int N=15,A=10306,C=21954,D=100000;
void ck(bool p,string s){if(!p)throw runtime_error(s);}
string code(vector<int>const&v){string s;for(int m:v){s+=to_string(m);s.push_back(',');}return s;}
struct Entry{bool exclude;bool used=false;};
int main(int argc,char**argv){try{
 ck(argc==5,"usage: metric_audit15 states.txt classes.txt residual.txt summary.txt");
 ifstream states(argv[1]),cert(argv[2]);ck(bool(states)&&bool(cert),"input missing");
 int B,I,FN;long long n;ck(bool(states>>B>>I>>FN>>n),"state header");ck(B+I==N&&10<=B&&B<=14&&FN==28-B&&n>0,"state header values");
 int cb,nclasses;ck(bool(cert>>cb>>nclasses)&&cb==B&&nclasses>=0,"class header");
 unordered_map<string,Entry> table;
 cpp_int S=cpp_int(1000000000000LL);
 for(int row=0;row<nclasses;row++){
  int tag,k;ck(bool(cert>>tag>>k)&&(tag==0||tag==1)&&0<=k&&k<=B*B,"class tag");vector<int> masks(k);
  for(int &m:masks)ck(bool(cert>>m)&&m>0&&m<(1<<B),"class mask");
  ck(is_sorted(masks.begin(),masks.end())&&adjacent_find(masks.begin(),masks.end())==masks.end(),"mask ordering");
  if(tag==1){
   vector<cpp_int> w(B+k);for(auto &v:w)ck(bool(cert>>v)&&v>=0,"negative/missing multiplier");
   for(int i=0;i<B;i++){cpp_int t=w[i];for(int j=0;j<k;j++)if(masks[j]>>i&1)t+=w[B+j];ck(t>=S,"Farkas coordinate deficit");}
   cpp_int cost=0;for(int i=0;i<B+k;i++)cost+=(i<B?A:C)*w[i];ck(cost<D*S,"non-strict Farkas inequality");
  }
  ck(table.emplace(code(masks),Entry{tag==1,false}).second,"duplicate classification");
 }
 string extra;ck(!(cert>>extra),"trailing classification");
 vector<array<int,3>> triples;for(int a=0;a<N;a++)for(int b=a+1;b<N;b++)for(int c=b+1;c<N;c++)triples.push_back({a,b,c});
 ofstream out(argv[3]);ck(bool(out),"residual output");long long direct=0,farkas=0,residual=0;
 for(long long idx=0;idx<n;idx++){
  array<unsigned,15> adj{};for(int i=0;i<B;i++)adj[i]=1u<<i;vector<int> ids(FN);
  for(int &id:ids){ck(bool(states>>id)&&0<=id&&id<(int)triples.size(),"bad triple id");auto t=triples[id];for(int u:t)if(u<B)for(int v:t)if(v<B)adj[u]|=1u<<v;}
  set<int> ms;bool bad=false;
  for(int i=0;i<B&&!bad;i++){
   unsigned reach=adj[i]|adj[(i+1)%B];
   for(int j=i+1;j<B;j++){
    if(!(reach&((1u<<j)|(1u<<((j+1)%B)))))continue;
    int k=j-i;if(min(k,B-k)<=2)continue;
    bool left=k*A<D-C,right=(B-k)*A<D-C;
    if(left&&right){bad=true;break;}
    int mask=((1<<j)-1)^((1<<i)-1);
    if(left)ms.insert(mask);else if(right)ms.insert(((1<<B)-1)^mask);
   }
  }
  if(bad){direct++;continue;}
  vector<int> masks(ms.begin(),ms.end());auto it=table.find(code(masks));ck(it!=table.end(),"unclassified signature");it->second.used=true;
  if(it->second.exclude){farkas++;continue;}
  residual++;out<<idx<<" "<<masks.size();for(int m:masks)out<<" "<<m;for(int id:ids)out<<" "<<id;out<<"\n";
 }
 ck(!(states>>extra),"trailing states");for(auto const &e:table)if(!e.second.exclude)ck(e.second.used,"unused residual signature");
 ck(direct+farkas+residual==n,"partition sum");ofstream report(argv[4]);ck(bool(report),"summary output");report<<B<<" "<<I<<" "<<n<<" "<<direct<<" "<<farkas<<" "<<residual<<"\n";
 cout<<"EXACT METRIC PARTITION B="<<B<<" N="<<n<<" direct="<<direct<<" farkas="<<farkas<<" residual="<<residual<<"\n";return 0;
 }catch(exception const&e){cerr<<"METRIC AUDIT FAILED: "<<e.what()<<"\n";return 1;}}
