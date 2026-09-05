#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
using boost::multiprecision::cpp_int;
using namespace std;

static constexpr int N=14, MAXF=16;
struct State {
    array<uint16_t,MAXF> f{};
    bool operator==(State const&o) const noexcept {return f==o.f;}
    bool operator<(State const&o) const noexcept {return f<o.f;}
};
struct HState { size_t operator()(State const&s) const noexcept {
    uint64_t h=1469598103934665603ULL;
    for(auto x:s.f){h^=x;h*=1099511628211ULL;}
    return (size_t)(h^(h>>32));
}};

int B,I,FN;
vector<array<int,N>> group_maps;
vector<array<int,N>> boundary_maps;
vector<array<int,4>> iperms;

inline uint16_t enc3(int a,int b,int c){
    if(a>b)swap(a,b); if(b>c)swap(b,c); if(a>b)swap(a,b);
    return (uint16_t)(a*N*N+b*N+c);
}
inline array<int,3> dec3(uint16_t z){
    int a=z/(N*N); z%=N*N; int b=z/N,c=z%N; return {a,b,c};
}
void sort_state(State& s){sort(s.f.begin(),s.f.begin()+FN); for(int i=FN;i<MAXF;i++)s.f[i]=65535;}

void gen_groups(){
    iperms.clear(); array<int,4> p={0,1,2,3};
    do {iperms.push_back(p);} while(next_permutation(p.begin(),p.begin()+I));
    boundary_maps.clear(); group_maps.clear();
    for(int refl=0;refl<2;refl++) for(int sh=0;sh<B;sh++){
        array<int,N> bm{};
        for(int v=0;v<B;v++) bm[v]= refl ? ((sh-v)%B+B)%B : (v+sh)%B;
        for(int v=B;v<N;v++) bm[v]=v;
        boundary_maps.push_back(bm);
        for(auto const& pp:iperms){
            auto m=bm; for(int k=0;k<I;k++)m[B+k]=B+pp[k]; group_maps.push_back(m);
        }
    }
}

using Sig=array<int,10>;
vector<Sig> boundary_signature(State const&s){
    array<int,N> deg{}; array<array<uint8_t,N>,N> adj{};
    array<array<int,3>,N> ftype{};
    for(int k=0;k<FN;k++){
        auto t=dec3(s.f[k]);
        int ni=(t[0]>=B)+(t[1]>=B)+(t[2]>=B);
        for(int u:t) if(u<B) ftype[u][ni]++;
        for(int x=0;x<3;x++)for(int y=x+1;y<3;y++){int u=t[x],v=t[y];adj[u][v]=adj[v][u]=1;}
    }
    for(int u=0;u<N;u++)for(int v=0;v<N;v++)deg[u]+=adj[u][v];
    vector<Sig> out(B);
    for(int u=0;u<B;u++){
        vector<int> ideg;
        for(int v=B;v<N;v++)if(adj[u][v])ideg.push_back(deg[v]);
        sort(ideg.begin(),ideg.end());
        Sig q{}; q[0]=deg[u];q[1]=(int)ideg.size();
        for(int k=0;k<4;k++)q[2+k]=k<(int)ideg.size()?ideg[k]:-1;
        q[6]=ftype[u][0];q[7]=ftype[u][1];q[8]=ftype[u][2];q[9]=ftype[u][3-1];
        out[u]=q;
    }
    return out;
}

State transform(State const&s,array<int,N> const&m){
    State z;
    for(int k=0;k<FN;k++){auto t=dec3(s.f[k]);z.f[k]=enc3(m[t[0]],m[t[1]],m[t[2]]);}sort_state(z);return z;
}

State canonical(State const&s){
    auto sig=boundary_signature(s);
    vector<int> keep; vector<Sig> bestseq;
    for(int bi=0;bi<(int)boundary_maps.size();bi++){
        vector<Sig> seq(B); auto const&m=boundary_maps[bi];
        for(int old=0;old<B;old++)seq[m[old]]=sig[old];
        if(keep.empty()||seq<bestseq){bestseq=seq;keep.assign(1,bi);} else if(seq==bestseq)keep.push_back(bi);
    }
    bool have=false; State best;
    for(int bi:keep){
        auto bm=boundary_maps[bi];
        for(auto const&pp:iperms){auto m=bm;for(int k=0;k<I;k++)m[B+k]=B+pp[k];State z=transform(s,m);if(!have||z<best){best=z;have=true;}}
    }
    if(!have)throw runtime_error("canonical failure"); return best;
}

bool validate(State const&s,string*why=nullptr){
    auto fail=[&](string q){if(why)*why=q;return false;};
    if(!is_sorted(s.f.begin(),s.f.begin()+FN))return fail("faces unsorted");
    for(int k=0;k<FN;k++){
        auto t=dec3(s.f[k]); if(!(0<=t[0]&&t[0]<t[1]&&t[1]<t[2]&&t[2]<N))return fail("bad face");
        if(k&&s.f[k]==s.f[k-1])return fail("duplicate face");
    }
    array<int,N*N> ec{}; array<array<uint8_t,N>,N> adj{}; array<int,N> used{};
    vector<vector<int>> faceadj(FN);
    array<array<int,2>,N*N> ef{}; for(auto&x:ef)x={-1,-1};
    for(int k=0;k<FN;k++){auto t=dec3(s.f[k]);for(int u:t)used[u]++;
      for(int x=0;x<3;x++)for(int y=x+1;y<3;y++){int u=t[x],v=t[y];if(u>v)swap(u,v);int e=u*N+v; if(ec[e]<2)ef[e][ec[e]]=k; ec[e]++;adj[u][v]=adj[v][u]=1;}}
    for(int u=0;u<N;u++)if(!used[u])return fail("unused vertex");
    int E=0;
    for(int u=0;u<N;u++)for(int v=u+1;v<N;v++)if(ec[u*N+v]){
      E++;bool bd=(u<B&&v<B&&((v==u+1)||(u==0&&v==B-1)));
      if(ec[u*N+v]!=(bd?1:2))return fail("edge multiplicity");
      if(ec[u*N+v]==2){int a=ef[u*N+v][0],b=ef[u*N+v][1];faceadj[a].push_back(b);faceadj[b].push_back(a);}
    }
    if(N-E+FN!=1)return fail("Euler");
    vector<int> seen(FN);queue<int>q;q.push(0);seen[0]=1;while(!q.empty()){int u=q.front();q.pop();for(int v:faceadj[u])if(!seen[v])seen[v]=1,q.push(v);}if(accumulate(seen.begin(),seen.end(),0)!=FN)return fail("face disconnected");
    // Vertex links.
    for(int v=0;v<N;v++){
      array<array<uint8_t,N>,N> la{};array<int,N> ld{};array<uint8_t,N> lv{};
      for(int k=0;k<FN;k++){auto t=dec3(s.f[k]);int pos=-1;for(int z=0;z<3;z++)if(t[z]==v)pos=z;if(pos<0)continue;int a=t[(pos+1)%3],b=t[(pos+2)%3];la[a][b]=la[b][a]=1;lv[a]=lv[b]=1;}
      int cnt=0,start=-1,ones=0;for(int a=0;a<N;a++)if(lv[a]){cnt++;start=a;for(int b=0;b<N;b++)ld[a]+=la[a][b];if(ld[a]==1)ones++;else if(ld[a]!=2)return fail("link degree");}
      if(v<B){if(ones!=2)return fail("boundary link endpoints");}else{if(ones!=0)return fail("interior link cycle");}
      array<uint8_t,N> vs{};queue<int>qq;qq.push(start);vs[start]=1;int got=0;while(!qq.empty()){int a=qq.front();qq.pop();got++;for(int b=0;b<N;b++)if(la[a][b]&&!vs[b])vs[b]=1,qq.push(b);}if(got!=cnt)return fail("link disconnected");
    }
    return true;
}

vector<State> neighbors(State const&s){
    array<int,N*N> ec{}; array<array<int,2>,N*N> ef{};for(auto&x:ef)x={-1,-1};
    for(int k=0;k<FN;k++){auto t=dec3(s.f[k]);for(int x=0;x<3;x++)for(int y=x+1;y<3;y++){int u=t[x],v=t[y];if(u>v)swap(u,v);int e=u*N+v;if(ec[e]<2)ef[e][ec[e]]=k;ec[e]++;}}
    vector<State> out;
    for(int a=0;a<N;a++)for(int b=a+1;b<N;b++){
      int e=a*N+b;if(ec[e]!=2)continue;int i=ef[e][0],j=ef[e][1];auto x=dec3(s.f[i]),y=dec3(s.f[j]);int c=-1,d=-1;for(int v:x)if(v!=a&&v!=b)c=v;for(int v:y)if(v!=a&&v!=b)d=v;if(c<0||d<0||c==d)continue;int u=min(c,d),v=max(c,d);if(ec[u*N+v])continue;
      State z=s;z.f[i]=enc3(c,d,a);z.f[j]=enc3(c,d,b);sort_state(z);out.push_back(canonical(z));
    }
    sort(out.begin(),out.end());out.erase(unique(out.begin(),out.end()),out.end());return out;
}

State seed(){
    State s;int n=0;for(int i=1;i<B-1;i++)s.f[n++]=enc3(0,i,i+1);
    for(int v=B;v<N;v++){
      auto t=dec3(s.f[0]);s.f[0]=enc3(v,t[0],t[1]);s.f[n++]=enc3(v,t[1],t[2]);s.f[n++]=enc3(v,t[2],t[0]);sort(s.f.begin(),s.f.begin()+n);
    }
    if(n!=FN)throw runtime_error("seed face count");sort_state(s);return canonical(s);
}

uint64_t orbit_size(State const&s){unordered_set<State,HState> im;im.reserve(group_maps.size()*2);for(auto const&m:group_maps)im.insert(transform(s,m));return im.size();}

unsigned long long fact(int n){unsigned long long z=1;for(int i=2;i<=n;i++)z*=i;return z;}
cpp_int factorial128(int n){cpp_int z=1;for(int i=2;i<=n;i++)z*=i;return z;}
unsigned long long brown_labeled(int b,int in){
    // 2(2B-3)!(4I+2B-5)! I! / ((B-1)!(B-3)! I!?(unlabelled formula includes I! denominator) (3I+2B-3)!)
    // Formula used in prior proof already returns labeled count after cancellation of I!: 2(2B-3)!(4I+2B-5)! / ((B-1)!(B-3)!(3I+2B-3)!).
    cpp_int num=2*factorial128(2*b-3)*factorial128(4*in+2*b-5);
    cpp_int den=factorial128(b-1)*factorial128(b-3)*factorial128(3*in+2*b-3);
    cpp_int q=num/den;if(num%den)throw runtime_error("Brown noninteger");return (unsigned long long)q;
}

int main(int argc,char**argv){
    if(argc<4){cerr<<"usage: enum B I outprefix\n";return 2;}B=stoi(argv[1]);I=stoi(argv[2]);string outp=argv[3];
    if(B+I!=N||B<3||I<0||I>4)throw runtime_error("bad B,I");FN=B+2*I-2;gen_groups();
    auto t0=chrono::steady_clock::now();State s0=seed();string why;if(!validate(s0,&why))throw runtime_error("invalid seed "+why);
    unordered_set<State,HState> seen;seen.reserve(500000);vector<State> reps;reps.reserve(400000);queue<uint32_t> q;seen.insert(s0);reps.push_back(s0);q.push(0);
    uint64_t attempts=0;
    while(!q.empty()){
      auto id=q.front();q.pop();auto ns=neighbors(reps[id]);attempts+=ns.size();
      for(auto const&z:ns)if(seen.insert(z).second){reps.push_back(z);q.push((uint32_t)reps.size()-1);}
      if(id%10000==0){double sec=chrono::duration<double>(chrono::steady_clock::now()-t0).count();cerr<<"orbits "<<reps.size()<<" queue "<<q.size()<<" sec "<<sec<<"\n";}
    }
    {
    ofstream f(outp+".txt"); f<<B<<" "<<I<<" "<<FN<<" "<<reps.size()<<"\n";
    int tid[N][N][N];int z=0;
    for(int a=0;a<N;a++)for(int b=a+1;b<N;b++)for(int c=b+1;c<N;c++)tid[a][b][c]=z++;
    for(auto const&r:reps){for(int k=0;k<FN;k++){auto t=dec3(r.f[k]);if(k)f<<" ";f<<tid[t[0]][t[1]][t[2]];}f<<"\n";}
    }
    uint64_t sumorb=0,minorb=~0ULL,maxorb=0;array<uint64_t,500> hist{};
    for(size_t k=0;k<reps.size();k++){
      if(!validate(reps[k],&why))throw runtime_error("invalid rep "+to_string(k)+" "+why);
      auto z=orbit_size(reps[k]);sumorb+=z;minorb=min(minorb,z);maxorb=max(maxorb,z);if(z<hist.size())hist[z]++;
    }
    auto expect=brown_labeled(B,I);if(sumorb!=expect)throw runtime_error("orbit sum mismatch "+to_string(sumorb)+" != "+to_string(expect));
    ofstream bin(outp+".bin",ios::binary);uint32_t magic=0x31334E45, nb=B,ni=I,nf=FN,nr=reps.size();bin.write((char*)&magic,4);bin.write((char*)&nb,4);bin.write((char*)&ni,4);bin.write((char*)&nf,4);bin.write((char*)&nr,4);for(auto&s:reps)bin.write((char*)s.f.data(),FN*sizeof(uint16_t));bin.close();
    double sec=chrono::duration<double>(chrono::steady_clock::now()-t0).count();
    ofstream js(outp+".json");js<<"{\n  \"verified\": true,\n  \"B\": "<<B<<",\n  \"I\": "<<I<<",\n  \"faces\": "<<FN<<",\n  \"orbits\": "<<reps.size()<<",\n  \"orbit_size_sum\": "<<sumorb<<",\n  \"brown_labeled\": "<<expect<<",\n  \"group_size\": "<<group_maps.size()<<",\n  \"min_orbit\": "<<minorb<<",\n  \"max_orbit\": "<<maxorb<<",\n  \"neighbor_attempts\": "<<attempts<<",\n  \"seconds\": "<<sec<<",\n  \"orbit_histogram\": {";bool first=true;for(size_t i=0;i<hist.size();i++)if(hist[i]){if(!first)js<<",";first=false;js<<"\n    \""<<i<<"\": "<<hist[i];}js<<"\n  }\n}\n";js.close();
    cout<<"VERIFIED B="<<B<<" I="<<I<<" orbits="<<reps.size()<<" sum="<<sumorb<<" seconds="<<sec<<"\n";
}
