#define PEEL_AS_LIBRARY 1
#include "peel_bfs20.cpp"

int main(int argc,char**argv){try{
 require(argc==5,"usage compare N B dfs_rows bfs_rows");n=stoi(argv[1]);b=stoi(argv[2]);interior=n-b;faces_required=2*n-b-2;require(n==20&&12<=b&&b<20,"dimensions");
 ifstream src(argv[3]);require(bool(src),"DFS input missing");size_t count;src>>count;require(bool(src),"DFS count missing");
 map<vector<uint32_t>,int> expected;
 for(size_t j=0;j<count;j++){
  Board s;s.fresh=n;s.triangles.resize(faces_required);
  for(auto &m:s.triangles){src>>m;require(bool(src),"DFS truncated");require(__builtin_popcount(m)==3&&m<(1u<<n),"invalid mask");for(int v=0;v<n;v++)if(m>>v&1)s.neighbors[v]|=m^(1u<<v);}
  Complex c(s);auto [key,stab]=c.normalize();require(expected.emplace(key,stab).second,"duplicate DFS orbit");
 }
 string junk;require(!(src>>junk),"DFS trailing tokens");ifstream obs(argv[4]);require(bool(obs),"BFS input missing");size_t count2;obs>>count2;require(count==count2,"different complete orbit counts");uint64_t rooted=0;
 for(size_t j=0;j<count2;j++){
  int stab;uint64_t mult;obs>>stab>>mult;require(bool(obs),"BFS record truncated");vector<uint32_t> key(faces_required);
  for(auto &m:key){obs>>m;require(bool(obs),"BFS faces truncated");}
  auto it=expected.find(key);require(it!=expected.end(),"BFS orbit absent or repeated");require(it->second==stab&&stab>0&&(2*b)%stab==0&&mult==uint64_t(2*b/stab),"stabilizer or multiplicity mismatch");rooted+=mult;expected.erase(it);
 }
 require(!(obs>>junk),"BFS trailing tokens");require(expected.empty(),"missing BFS orbit");
 cout<<"{\"verified\":true,\"N\":"<<n<<",\"B\":"<<b<<",\"orbits\":"<<count<<",\"rooted\":"<<rooted<<"}\n";
}catch(exception const&e){cerr<<"REJECTED: "<<e.what()<<"\n";return 1;}}
