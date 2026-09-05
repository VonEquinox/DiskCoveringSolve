// Exact orbit-mass auditor for disk triangulation representatives.
// Input format:
//   first line: B I NF NSTATE
//   each following line: comma-separated lexicographic triple indices.
// All arithmetic and comparisons are integral.

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;
using Tri = std::array<int,3>;
using Code = std::array<uint64_t,6>; // C(14,3)=364

static int B,I,N,NF;
static std::vector<Tri> triples;
static int tid[14][14][14];

static Tri normtri(int a,int b,int c){
    Tri t{a,b,c}; std::sort(t.begin(),t.end()); return t;
}
static bool code_less(const Code&a,const Code&b){
    // Compare as a 384-bit integer in a stable lexicographic order.
    for(int k=5;k>=0;--k){
        if(a[k]!=b[k]) return a[k]<b[k];
    }
    return false;
}
static std::string code_key(const Code&c){
    return std::string(reinterpret_cast<const char*>(c.data()),sizeof(Code));
}
static cpp_int fact(int n){ cpp_int z=1; for(int k=2;k<=n;k++)z*=k; return z; }
static cpp_int brown(int b,int i){
    // 2(2B-3)!(4I+2B-5)! / ((B-1)!(B-3)!(3I+2B-3)!)
    return cpp_int(2)*fact(2*b-3)*fact(4*i+2*b-5)/
           (fact(b-1)*fact(b-3)*fact(3*i+2*b-3));
}
struct CodeHash{
    size_t operator()(Code const&c)const noexcept{
        uint64_t h=0x9e3779b97f4a7c15ULL;
        for(auto x:c){x^=x>>30;x*=0xbf58476d1ce4e5b9ULL;x^=x>>27;
            x*=0x94d049bb133111ebULL;x^=x>>31;h^=x+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);}
        return size_t(h);
    }
};
static Code make_code(const std::vector<int>&ids){
    Code c{}; for(int x:ids)c[x>>6]|=uint64_t(1)<<(x&63); return c;
}
static std::vector<int> parse_ids(const std::string&s){
    std::vector<int>a; std::string z;
    for(char ch:s){
        if(ch==','||ch==' '||ch=='\t'){if(!z.empty()){a.push_back(std::stoi(z));z.clear();}}
        else z.push_back(ch);
    }
    if(!z.empty())a.push_back(std::stoi(z));
    return a;
}
static uint64_t epair(int a,int b){if(a>b)std::swap(a,b);return (uint64_t(a)<<32)|uint32_t(b);}
static bool is_boundary_edge(int a,int b){
    if(a>=B||b>=B)return false;
    int d=std::abs(a-b); return d==1||d==B-1;
}
static void validate(const std::vector<int>&ids){
    if((int)ids.size()!=NF)throw std::runtime_error("wrong face count");
    std::set<int> uniq(ids.begin(),ids.end());
    if((int)uniq.size()!=NF)throw std::runtime_error("duplicate face");
    std::unordered_map<uint64_t,int> ec;
    std::vector<std::vector<int>> adj(N);
    std::vector<int> used(N,0);
    for(int id:ids){
        if(id<0||id>=(int)triples.size())throw std::runtime_error("bad triple id");
        auto t=triples[id]; for(int v:t)used[v]=1;
        for(int k=0;k<3;k++){
            int a=t[k],b=t[(k+1)%3];ec[epair(a,b)]++;
        }
    }
    for(int v=0;v<N;v++)if(!used[v])throw std::runtime_error("unused vertex");
    for(auto [e,c]:ec){
        int a=int(e>>32),b=int(uint32_t(e));
        if(is_boundary_edge(a,b)){if(c!=1)throw std::runtime_error("boundary edge multiplicity");}
        else if(c!=2)throw std::runtime_error("interior edge multiplicity");
        adj[a].push_back(b);adj[b].push_back(a);
    }
    int E=int(ec.size());
    if(N-E+NF!=1)throw std::runtime_error("Euler failure");
    // 1-skeleton connectivity.
    std::vector<int>vis(N);std::queue<int>q;q.push(0);vis[0]=1;
    while(!q.empty()){int u=q.front();q.pop();for(int v:adj[u])if(!vis[v])vis[v]=1,q.push(v);}
    if(std::accumulate(vis.begin(),vis.end(),0)!=N)throw std::runtime_error("disconnected");
    // Vertex links.
    for(int v=0;v<N;v++){
        std::unordered_map<int,std::vector<int>> la;
        std::set<int> ln;
        for(int id:ids){
            auto t=triples[id];
            int pos=-1;for(int k=0;k<3;k++)if(t[k]==v)pos=k;
            if(pos>=0){
                int a=t[(pos+1)%3],b=t[(pos+2)%3];
                la[a].push_back(b);la[b].push_back(a);ln.insert(a);ln.insert(b);
            }
        }
        if(ln.empty())throw std::runtime_error("empty link");
        int deg1=0;for(int x:ln){
            int d=int(la[x].size());
            if(v<B){if(d==1)deg1++;else if(d!=2)throw std::runtime_error("bad boundary link degree");}
            else if(d!=2)throw std::runtime_error("bad interior link degree");
        }
        if(v<B&&deg1!=2)throw std::runtime_error("boundary link not path");
        if(v>=B&&deg1!=0)throw std::runtime_error("interior link not cycle");
        int start=*ln.begin(),cnt=0;std::set<int>seen;std::queue<int>qq;qq.push(start);seen.insert(start);
        while(!qq.empty()){int u=qq.front();qq.pop();cnt++;for(int w:la[u])if(!seen.count(w))seen.insert(w),qq.push(w);}
        if(cnt!=(int)ln.size())throw std::runtime_error("disconnected link");
        if(v<B){
            int pm=(v+B-1)%B,pp=(v+1)%B;
            if(la[pm].size()!=1||la[pp].size()!=1)throw std::runtime_error("wrong boundary link endpoints");
        }
    }
}
struct GMap{std::array<int,14> p; std::vector<int> tmap;};
int main(int argc,char**argv){
    if(argc<2){std::cerr<<"usage: audit states_file [report.json]\n";return 2;}
    std::ifstream f(argv[1]);if(!f)throw std::runtime_error("open failed");
    long long declared;
    f>>B>>I>>NF>>declared;std::string line;std::getline(f,line);
    N=B+I;if(N!=14||B<10||B>13||I!=14-B||NF!=26-B||declared<=0)throw std::runtime_error("bad family/header");
    for(int a=0;a<N;a++)for(int b=a+1;b<N;b++)for(int c=b+1;c<N;c++){
        tid[a][b][c]=int(triples.size());triples.push_back({a,b,c});
    }
    // All group elements D_B x S_I, with precomputed triple maps.
    std::vector<int> ip(I);std::iota(ip.begin(),ip.end(),0);
    std::vector<std::vector<int>> perms;
    do{perms.push_back(ip);}while(std::next_permutation(ip.begin(),ip.end()));
    std::vector<GMap> G;
    for(int refl=0;refl<2;refl++)for(int sh=0;sh<B;sh++)for(auto const&pm:perms){
        GMap g;
        for(int x=0;x<N;x++){
            if(x<B){int y=refl?(sh-x)%B:(x+sh)%B;if(y<0)y+=B;g.p[x]=y;}
            else g.p[x]=B+pm[x-B];
        }
        g.tmap.resize(triples.size());
        for(int k=0;k<(int)triples.size();k++){
            auto t=triples[k];auto u=normtri(g.p[t[0]],g.p[t[1]],g.p[t[2]]);
            g.tmap[k]=tid[u[0]][u[1]][u[2]];
        }
        G.push_back(std::move(g));
    }
    const long long gsize=2LL*B*(long long)perms.size();
    if((long long)G.size()!=gsize)throw std::runtime_error("group size");
    cpp_int mass=0;
    std::unordered_set<Code,CodeHash> canonical_seen;
    long long nstate=0,minstab=1e18,maxstab=0;
    while(std::getline(f,line)){
        if(line.empty())continue;
        auto ids=parse_ids(line);validate(ids);
        Code orig=make_code(ids),best{};bool have=false;long long stab=0;
        for(auto const&g:G){
            Code c{};for(int id:ids){int x=g.tmap[id];c[x>>6]|=uint64_t(1)<<(x&63);}
            if(c==orig)stab++;
            if(!have||code_less(c,best)){best=c;have=true;}
        }
        if(stab<=0||gsize%stab)throw std::runtime_error("bad stabilizer");
        if(!canonical_seen.insert(best).second)throw std::runtime_error("duplicate orbit representative");
        mass+=gsize/stab;nstate++;minstab=std::min(minstab,stab);maxstab=std::max(maxstab,stab);
        if(nstate%10000==0)std::cerr<<"checked "<<nstate<<"\n";
    }
    if(declared>=0&&nstate!=declared)throw std::runtime_error("declared count mismatch");
    cpp_int target=brown(B,I);bool ok=(mass==target);
    std::ostringstream js;
    js<<"{\n  \"verified\": "<<(ok?"true":"false")
      <<",\n  \"B\": "<<B<<",\n  \"I\": "<<I
      <<",\n  \"representatives\": "<<nstate
      <<",\n  \"group_size\": "<<gsize
      <<",\n  \"orbit_mass\": \""<<mass<<"\""
      <<",\n  \"brown_total\": \""<<target<<"\""
      <<",\n  \"min_stabilizer\": "<<minstab
      <<",\n  \"max_stabilizer\": "<<maxstab<<"\n}\n";
    std::cout<<js.str();
    if(argc>=3){std::ofstream o(argv[2]);o<<js.str();}
    return ok?0:1;
}
