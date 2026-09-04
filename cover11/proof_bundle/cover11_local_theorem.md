# Strict local lower bound at the 11-disk candidate

Let \(G_{45}\) be the 45-edge core auxiliary graph of the candidate topology,
with nine circular anchors and twenty free vertices. The full auxiliary graph
\(G_{\mathrm{full}}\) additionally contains the two slack-face witness vertices
and their six incident edges. Restricting a full-graph realization to the core
gives, for every anchor-gap vector,
\[
 \rho_{G_{\mathrm{full}}}(g)^2\ge \rho_{G_{45}}(g)^2.
\]
Thus a lower bound proved for \(G_{45}\) is also a lower bound for the full
candidate topology.

For the local analysis use eight independent gap coordinates
\[
 u=(u_0,\ldots,u_7)\in\mathbb R^8,
 \qquad
 \Phi(u)=\left(u_0,\ldots,u_7,
 2\pi-\sum_{i=0}^7u_i\right).
\]
The ninth circular gap is therefore determined by the total angle. Write
\(g=\Phi(u)\), \(g_*=\Phi(u_*)\), and
\[
 \rho_{G_{45}}(\Phi(u))^2
 =\min_x\max_{vw\in E(G_{45})}|x_v-x_w|^2.
\]
All gradients, Hessians, orthogonal complements, and Euclidean norms below are
taken in \(\mathbb R^8\).

The full active KKT certificate gives an exact algebraic candidate
\((g_*,x_*,t_*)\). Its 45 edge multipliers are strictly positive, sum to
one, equilibrate every free vertex, and have zero tangential anchor
gradient. The equilibrium-plus-normalization matrix has rank 41, so its
affine stress space has dimension four.

Six nearby positive stresses are obtained by correcting rational proposal
stresses with
\[
 w_i^*=w_i^0+B_*^T(B_*B_*^T)^{-1}(b-B_*w_i^0).
\]
The rational inverse and root-box bounds prove
\[
 \|w_i^*-w_i^0\|_\infty<2\cdot10^{-37};
\]
hence every corrected stress remains positive. The graph-energy
Lipschitz lemma then bounds every effective-conductance perturbation by
\(10^{-34}\).

For each corrected stress define its Dirichlet lower-bound energy
\[
 E_i(g)=\min_x\sum_{vw}w_{i,vw}^*|x_v-x_w|^2,
 \qquad
 \widetilde E_i(u)=E_i(\Phi(u)).
\]
At the candidate all active edges have squared length \(t_*\), and the
candidate free vertices are harmonic, so
\[
 \widetilde E_i(u_*)=E_i(g_*)=t_*.
\]
The envelope theorem makes \(\nabla\widetilde E_i(u_*)\) linear in the stress.
Because the exact KKT stress has zero gradient and the affine stress space
is four-dimensional, the six gradients span a space
\[
 A_*:=\operatorname{span}\{\nabla\widetilde E_1(u_*),\ldots,
 \nabla\widetilde E_6(u_*)\}\subset\mathbb R^8
\]
of dimension at most four. The certified singular-value lower bound proves the
dimension is exactly four. Put
\[
 N_*:=A_*^\perp\subset\mathbb R^8.
\]
Consequently \(\dim A_*=\dim N_*=4\).

The exact rational central audit, together with the transfer bounds, gives
for all six energies
\[
\begin{aligned}
 \min_{|a|=1,a\in A_*}\max_i\nabla\widetilde E_i(u_*)\cdot a&\ge0.0013,\\
 n^T\nabla^2\widetilde E_i(u_*)n&\ge0.0026|n|^2,\\
 |a^T\nabla^2\widetilde E_i(u_*)n|&\le0.058|a||n|,\\
 |a^T\nabla^2\widetilde E_i(u_*)a|&\le0.57|a|^2,
\end{aligned}
\]
and the trilinear third-derivative norm is at most \(3.07\).

For \(h=u-u_*=a+n\), \(a\in A_*\), \(n\in N_*\), and \(\sigma=|h|\),
Taylor's theorem therefore yields
\[
 \rho_{G_{45}}(\Phi(u))^2-t_*
 \ge .0013|a|+.0013|n|^2-.058|a||n|
      -.285|a|^2-\frac{3.07}{6}\sigma^3.                 \tag{1}
\]
Using \(|a||n|\le |a|\sigma\), the right side at fixed \(\sigma\) is bounded
below by a concave quadratic in \(|a|\); its minimum is at an endpoint.
For \(0<\sigma\le 11/5000\), the two endpoint coefficients are
\[
 \frac{0.0026}{2}-\frac{3.07}{6}\sigma>0,
\]
and
\[
 0.0013-(0.058+0.57/2)\sigma-\frac{3.07}{6}\sigma^2>0.
\]
At \(\sigma=11/5000\) their exact rational margins are respectively
\(523/3000000\) and \(8143853/15000000000\). Consequently
\[
 \boxed{\rho_{G_{45}}(\Phi(u))^2>t_*\quad
  (0<|u-u_*|\le 0.0022).}
\]

The external branch-and-bound certificate measures the full nine-gap vector.
For \(h=u-u_*\),
\[
 \|\Phi(u)-\Phi(u_*)\|_2^2
 =\|h\|_2^2+\left(\sum_{i=0}^7h_i\right)^2
 \ge \|h\|_2^2.
\]
Therefore every certified local leaf contained in the radius-\(0.0022\)
nine-gap ball also lies in the eight-dimensional local ball required above.
Together with the exact external branch-and-bound certificate and
\(\rho_{G_{\mathrm{full}}}\ge\rho_{G_{45}}\), this proves that the candidate
topology has global minimum exactly \(t_*\).
