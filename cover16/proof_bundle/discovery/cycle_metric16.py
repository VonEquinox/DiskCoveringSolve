"""Exact integer negative-cycle finder for cyclic interval-cap constraints.
Its output is independently checked as the old nonnegative Farkas covering vector.
"""
import numpy as np
from numba import njit
A=9974;C=21143;C6=37566;D=100000;S=10**12
@njit(cache=True)
def cycle_weights(B,key):
 n=2*B+len(key);us=np.zeros(n,np.int64);vs=np.zeros(n,np.int64);ws=np.zeros(n,np.int64);inds=np.full(n,-1,np.int64);wr=np.zeros(n,np.int64)
 for i in range(B):
  us[i]=i;vs[i]=(i+1)%B;wr[i]=int(i==B-1);ws[i]=A-D*wr[i];inds[i]=i
  us[B+i]=(i+1)%B;vs[B+i]=i;wr[B+i]=-int(i==B-1);ws[B+i]=-D*wr[B+i]
 for k,encoded in enumerate(key):
  m=encoded&((1<<B)-1);cap=C6 if encoded>>B else C
  i=0;j=0
  for v in range(B):
   if (m>>v&1) and not(m>>((v-1)%B)&1):i=v
   if not(m>>v&1) and (m>>((v-1)%B)&1):j=v
  e=2*B+k;us[e]=i;vs[e]=j;wr[e]=int(j<i);ws[e]=cap-D*wr[e];inds[e]=B+k
 dist=np.zeros(B,np.int64);pred=np.full(B,-1,np.int64);last=-1
 for it in range(B):
  last=-1
  for e in range(n):
   u=us[e];v=vs[e]
   if dist[v]>dist[u]+ws[e]:dist[v]=dist[u]+ws[e];pred[v]=e;last=v
  if last<0:return np.zeros(0,np.int64)
 v=last
 for i in range(B):v=us[pred[v]]
 start=v;out=np.zeros(B+len(key),np.int64);winding=0;cost=0
 for j in range(B+1):
  e=pred[v];winding+=wr[e];cost+=ws[e]
  if inds[e]>=0:out[inds[e]]+=1
  v=us[e]
  if v==start:break
 assert v==start and winding>0 and cost<0
 den=(S+winding-1)//winding
 out*=den
 # A separate direct cover/strict-cost check, repeated by the final verifier.
 for i in range(B):
  cover=out[i]
  for k,m in enumerate(key):
   if m>>i&1:cover+=out[B+k]
  assert cover>=S
 cost2=A*np.sum(out[:B])
 for k in range(len(key)):cost2+=(C6 if key[k]>>B else C)*out[B+k]
 assert cost2<D*S
 return out
