#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using boost::multiprecision::cpp_int;
using Tri=array<uint8_t,3>;

struct State{ vector<Tri> t; };
static int B,I,N,F;

static inline Tri mk(int a,int b,int c){ Tri z{(uint8_t)a,(uint8_t)b,(uint8_t)c}; sort(z.begin(),z.end()); return z; }
static inline void norm(vector<Tri>&t){ for(auto &z:t)sort(z.begin(),z.end()); sort(t.begin(),t.end()); }
static inline string key_of(const vector<Tri>&t){ string s; s.resize(3*t.size()); size_t p=0; for(auto z:t)for(int k=0;k<3;k++)s[p++]=char(z[k]); return s; }

static vector<vector<int>> perms_vec(vector<int> a){
 sort(a.begin(),a.end()); vector<vector<int>> out; do{out.push_back(a);}while(next_permutation(a.begin(),a.end())); return out;
}

struct Canon{State s; string key; int stab;};

static Canon canonical(const State& in){
 vector<int> deg(N,0);
 vector<unordered_set<int>> nb(N);
 for(auto z:in.t){
  for(int a=0;a<3;a++)for(int b=a+1;b<3;b++){nb[z[a]].insert(z[b]);nb[z[b]].insert(z[a]);}
 }
 for(int v=0;v<N;v++)deg[v]=(int)nb[v].size();
 // boundary dihedral maps old->new, retain only minimal transformed boundary degree word.
 vector<vector<int>> bmaps; vector<int> bestbd;
 for(int rev=0;rev<2;rev++)for(int sh=0;sh<B;sh++){
  vector<int> mp(B),word(B);
  for(int x=0;x<B;x++){int y=rev? (sh-x)%B:(sh+x)%B;if(y<0)y+=B;mp[x]=y;word[y]=deg[x];}
  if(bestbd.empty()||word<bestbd){bestbd=word;bmaps.clear();bmaps.push_back(mp);} else if(word==bestbd)bmaps.push_back(mp);
 }
 // Interior degree slots are sorted ascending. Enumerate only within equal-degree groups.
 map<int,vector<int>> groups;
 for(int x=B;x<N;x++)groups[deg[x]].push_back(x);
 vector<pair<vector<int>,vector<vector<int>>>> gp;
 int slot=B;
 for(auto &kv:groups){
  vector<int> slots;for(size_t k=0;k<kv.second.size();k++)slots.push_back(slot++);
  gp.push_back({slots,perms_vec(kv.second)}); // permutation lists old vertices assigned in slot order
 }
 vector<Tri> best; int count=0;
 vector<int> mp(N);
 function<void(int,const vector<int>&)> rec=[&](int gi,const vector<int>& bmap){
  if(gi==(int)gp.size()){
   vector<Tri> tt;tt.reserve(in.t.size());
   for(auto z:in.t)tt.push_back(mk(mp[z[0]],mp[z[1]],mp[z[2]]));norm(tt);
   if(best.empty()||tt<best){best=move(tt);count=1;}else if(tt==best)count++;
   return;
  }
  auto &slots=gp[gi].first;
  for(auto &oldperm:gp[gi].second){
   for(size_t k=0;k<slots.size();k++)mp[oldperm[k]]=slots[k];
   rec(gi+1,bmap);
  }
 };
 for(auto &bm:bmaps){for(int x=0;x<B;x++)mp[x]=bm[x];rec(0,bm);}
 return {State{best},key_of(best),count};
}

static vector<State> neighbors(const State&s){
 struct Inc{int ti,opp;};
 unordered_map<int,vector<Inc>> inc; inc.reserve(3*F*2);
 unordered_set<int> edges;edges.reserve(3*F*2);
 auto ek=[](int a,int b){if(a>b)swap(a,b);return a*N+b;};
 for(int ti=0;ti<(int)s.t.size();ti++){
  auto z=s.t[ti];
  for(int k=0;k<3;k++){int a=z[(k+1)%3],b=z[(k+2)%3],o=z[k];int e=ek(a,b);inc[e].push_back({ti,o});edges.insert(e);}
 }
 vector<State> out;
 for(auto &kv:inc){auto &v=kv.second;if(v.size()!=2)continue;
  int e=kv.first,a=e/N,b=e%N,c=v[0].opp,d=v[1].opp;
  // never flip prescribed outer boundary edge (it has one incident anyway); link condition.
  if(c==d||edges.count(ek(c,d)))continue;
  State q=s;q.t[v[0].ti]=mk(c,d,a);q.t[v[1].ti]=mk(c,d,b);norm(q.t);out.push_back(move(q));
 }
 return out;
}

static bool validate(const State&s,string &why){
 if((int)s.t.size()!=F){why="face count";return false;}
 set<Tri> ts(s.t.begin(),s.t.end());if((int)ts.size()!=F){why="duplicate face";return false;}
 vector<int> seen(N);map<pair<int,int>,int> ec;
 for(auto z:s.t){if(z[0]==z[1]||z[1]==z[2]){why="degenerate";return false;}for(int x:z)seen[x]=1;
  for(int k=0;k<3;k++){int a=z[k],b=z[(k+1)%3];if(a>b)swap(a,b);ec[{a,b}]++;}}
 if(accumulate(seen.begin(),seen.end(),0)!=N){why="missing vertex";return false;}
 for(auto &x:ec)if(x.second<1||x.second>2){why="edge incidence";return false;}
 set<pair<int,int>> bd;for(int i=0;i<B;i++){int j=(i+1)%B;bd.insert(minmax(i,j));}
 set<pair<int,int>> got;for(auto &x:ec)if(x.second==1)got.insert(x.first);
 if(got!=bd){why="boundary";return false;}
 int E=ec.size();if(N-E+F!=1){why="Euler";return false;}
 return true;
}

static cpp_int fact(int n){cpp_int x=1;for(int i=2;i<=n;i++)x*=i;return x;}
static cpp_int brown(int b,int i){return 2*fact(2*b-3)*fact(4*i+2*b-5)/(fact(b-1)*fact(b-3)*fact(3*i+2*b-3));}

static State seed_state(){
 State s;
 for(int k=1;k<=B-2;k++)s.t.push_back(mk(0,k,k+1));
 for(int v=B;v<N;v++){
  // insert into lexicographically first face
  norm(s.t);Tri z=s.t.front();s.t.erase(s.t.begin());
  s.t.push_back(mk(z[0],z[1],v));s.t.push_back(mk(z[1],z[2],v));s.t.push_back(mk(z[2],z[0],v));
 }
 norm(s.t);return s;
}

int main(int argc,char**argv){
 if(argc<4){cerr<<"usage: enum B I output.jsonl [summary.json]\n";return 2;}
 B=stoi(argv[1]);I=stoi(argv[2]);N=B+I;F=B+2*I-2;
 string outpath=argv[3],summ=argc>4?argv[4]:outpath+".summary.json";
 auto t0=chrono::steady_clock::now();
 Canon c0=canonical(seed_state());
 vector<State> states;vector<int> stabs;states.push_back(c0.s);stabs.push_back(c0.stab);
 unordered_map<string,int> seen;seen.reserve(1000000);seen.emplace(c0.key,0);
 size_t head=0;unsigned long long attempts=0;
 while(head<states.size()){
  State cur=states[head++];
  auto ns=neighbors(cur);attempts+=ns.size();
  for(auto &q:ns){Canon c=canonical(q);auto [it,ok]=seen.emplace(c.key,(int)states.size());if(ok){states.push_back(move(c.s));stabs.push_back(c.stab);}}
  if(head%10000==0){double sec=chrono::duration<double>(chrono::steady_clock::now()-t0).count();cerr<<"states "<<head<<" total "<<states.size()<<" attempts "<<attempts<<" sec "<<sec<<"\n";}
 }
 cpp_int osum=0; long long G=2LL*B;for(int k=2;k<=I;k++)G*=k;
 map<int,int> sh;
 for(int s:stabs){if(G%s){cerr<<"bad stab "<<s<<" group "<<G<<"\n";return 3;}osum+=G/s;sh[s]++;}
 cpp_int br=brown(B,I);bool match=(osum==br);
 // validate and emit records
 ofstream fo(outpath);if(!fo){cerr<<"cannot open output\n";return 4;}
 for(size_t idx=0;idx<states.size();idx++){
  string why;if(!validate(states[idx],why)){cerr<<"invalid "<<idx<<" "<<why<<"\n";return 5;}
  fo<<"{\"idx\":"<<idx<<",\"stab\":"<<stabs[idx]<<",\"faces\":[";
  for(size_t k=0;k<states[idx].t.size();k++){if(k)fo<<',';auto z=states[idx].t[k];fo<<'['<<(int)z[0]<<','<<(int)z[1]<<','<<(int)z[2]<<']';}
  fo<<"]}\n";
 }
 fo.close();
 { vector<Tri> trip; map<Tri,int> id;
 for(int a=0;a<N;a++)for(int b=a+1;b<N;b++)for(int c=b+1;c<N;c++){Tri z=mk(a,b,c);id[z]=trip.size();trip.push_back(z);}
 ofstream txt(outpath+".txt");txt<<B<<" "<<I<<" "<<F<<" "<<states.size()<<"\n";
 for(const auto& st:states){for(size_t k=0;k<st.t.size();k++){if(k)txt<<",";txt<<id[st.t[k]];}txt<<"\n";}
 }
 double sec=chrono::duration<double>(chrono::steady_clock::now()-t0).count();
 ofstream fs(summ);fs<<"{\n  \"B\": "<<B<<", \"I\": "<<I<<",\n  \"orbits\": "<<states.size()<<",\n  \"group_size\": "<<G<<",\n  \"orbit_size_sum\": \""<<osum<<"\",\n  \"brown_labeled\": \""<<br<<"\",\n  \"brown_match\": "<<(match?"true":"false")<<",\n  \"attempted_flips\": "<<attempts<<",\n  \"stabilizer_histogram\": {";
 bool first=true;for(auto &kv:sh){if(!first)fs<<',';first=false;fs<<"\""<<kv.first<<"\":"<<kv.second;}fs<<"},\n  \"seconds\": "<<sec<<"\n}\n";fs.close();
 cerr<<"DONE B "<<B<<" I "<<I<<" orbits "<<states.size()<<" orbit_sum "<<osum<<" brown "<<br<<" match "<<match<<" sec "<<sec<<"\n";
 return match?0:6;
}
