# The theory

Lie Groups and Differential Geomtry

The Lie Groups SO(3) and SE(3)

Riemann Geometry

Homotopy 


## Lie Groups and Differential Geometry

__Definition__ :
A _tangent vector_ at a point $p$ in a manifold $M$ is a _derivation_ at $p$.
The tangent vectors at $p$ form a vector space $T_pM$ called the _tangent space_ of $M$ at $p$.

__Definition__ :
A _Lie group_ is a subset $G$ of $\mathbb{R}^n$ such that $G$ is a group and a manifold in $\mathbb{R}^n$ and both the group operation $\cdot : G\rightarrow G$ and the inverse operation $.^{-1}: G\rightarrow G$ are smooth functions.

__Definition__ :
A _Lie algebra_ is an algebra $A$ togther with a _Lie bracket_ operator $\[.,.\] : A \times A \rightarrow A$.
An important fact is that the Lie algebra $A$ associated to a Lie group $M$ happens to be the tangent space
at the identity element $1$, that is: $A = T_1M$.

__Definition__ : _Lie bracket_



## SO(3) and SE(3) are Lie groups


## Riemannian Geometry
__Definition__ : 
Let $M$ be an $n$-dimensional manifold.
A _Riemannian metric_ on $M$ is a smoothly varying positive-definite inner-product $g_p$ on each tangent space $T_pM$.
We say that the pair $(M,g)$ is a _Riemannian manifold_.

We now define what is the analogue of "straight lines" in a Riemannian manifold $(M,g)$ :

__Definition__ : 
Let $\gamma : I \rightarrow M$ be a regular path.
Consider $D_t: \mathcal{X}(\gamma) \rightarrow \mathcal{X}(\gamma)$ the _covariant derivative along_ $\gamma$.
The velocity of $\gamma$ (its derivation $t\mapsto \gamma'(t)$) defines a vector field along $\gamma$, we call $D_t\gamma'$ the _acceleration_ of $\gamma$. 
When $D_t\gamma' = 0$, we say that $\gamma$ is a _geodesic_.

The idea is now to use geodesics to "explore" $M$, imagine sending probes with zero acceleration with all different velocities, after one second they report back their positions, giving you a "map" of $M$. Formally :

__Definition__ : 
Let $p\in M$ and $v\in T_pM$. We write $\gamma_v:I\rightarrow M$ for the (maximal) geodesic with $\gamma(0) = p$ and $\gamma'(0) = v$. We denote by $U_p \subset T_pM$ the set of those $v$ for which $\gamma_v(1)$ is well-defined. 
The _exponential map_ at $p$, written $exp_p:U_p\rightarrow M$, is defined by $exp_p(v) = \gamma_v(1)$.

__Remark__ : 
In the context of Lie groups, the exponential map from the Lie algebra (tangent space) to the corresponding Lie group (manifold).
There is an "inverse" computation that we call the _logarithm map_ which sends an element from the Lie group (manifold) to the Lie algebra (tangent space).


## Homotopy and Homology theory

__Definition__ : A _homotopy_ between two continuous functions $f$ and $g$ from $X$ to $Y$ is a continuous function $H : X \times [0,1]  \rightarrow Y$ such that $H(x,0) = f(x)$ and $H(x,1) = g(x)$ for all $x\in X$.

__Definition__ : Two paths $\gamma_1$ and $\gamma_2$ are said to be _homotopic_ if there exists a homotopy taking $\gamma_1$ to $\gamma_2$.

__Definition__ : The _fundamental group_ of a pointed topological space $(X,x)$ denoted $\pi_1(X,x)$ is the group of equivalence classes under homotopy of the loops based at $x$ in $X$.

__Example__ : The fundamental group of $\mathbb{R}^2$ is trivial, $\pi_1(\mathbb{R}^2 - \{p\}) = \mathbb{Z}$ and $\pi_1(\mathbb{R}^2 - \{p_1, p_2, ..., p_n\}) = F_n$ the free group on $n$ generators.


__Definition__ : Let $X$ be a topological space, _chain complexe_ $(C_*, d_*)$ with $d_n : C_n \rightarrow C_{n-1}, each $C_n$ is an abelian group and $d_{n-1} \circ d_n = 0 \forall n$.
$B_n := \im d_{n+1} and $Z_n := \ker d_n$. 
Finally $H_n(X) := Z_n /B_n ...


__Theorem (Hurewicz)__ : For $X$ connected, there exists a homomorphism $h : \pi_1(X) \rightarrow H_1(X)$ which is surjective and $\ker h = [\pi_1(X), \pi_1(X)]$ the commutator subgroup. This implies that the abelianisation of the first fundamental group is isomorphic to the first homology group, $\pi_1(X)^{ab} \cong H_1(X)$.