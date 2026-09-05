#pragma once
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
constexpr int MAXN=15,MAXF=18;
int N,B,I,F;
struct State {array<uint16_t,MAXF> t{};};
bool eq(const State&a,const State&b){return a.t==b.t;}
uint64_t hashstate(const State&s){uint64_t h=1469598103934665603ULL;for(int i=0;i<F;i++){h^=s.t[i];h*=1099511628211ULL;} h^=h>>33;h*=0xff51afd7ed558ccdULL;h^=h>>33;return h;}
struct Store{
 vector<State> states;vector<uint32_t> table;size_t mask;
 Store(size_t estimate){states.reserve(estimate);size_t sz=1;while(sz<estimate*1.6)sz*=2;table.assign(sz,0);mask=sz-1;}
 void rehash(){vector<uint32_t> nt(table.size()*2,0);size_t m=nt.size()-1;for(uint32_t i=0;i<states.size();i++){size_t h=hashstate(states[i])&m;while(nt[h])h=(h+1)&m;nt[h]=i+1;}table.swap(nt);mask=m;}
 bool insert(const State&s){if((states.size()+1)*10>table.size()*7)rehash();size_t h=hashstate(s)&mask;while(table[h]){if(eq(states[table[h]-1],s))return false;h=(h+1)&mask;}states.push_back(s);table[h]=states.size();return true;}
};
