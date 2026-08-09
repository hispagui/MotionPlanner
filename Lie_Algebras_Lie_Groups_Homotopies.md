# The theory
layout :

Lie Groups and Differential Geomtry

The Lie Groups SO(3) and SE(3)

Riemann Geometry

Homotopy 


## Lie Groups and Differential Geometry
__Definition__ :
A group is a set $G$ together with a binary operation on ⁠$G$⁠, here denoted "$\cdot$", that combines any two elements $a$ and $b$ of $G$ to form an element of $G$ denoted $a \cdot b$⁠, such that the following three requirements, known as group axioms, are satisfied:

Associativity : $\forall a,b \in G, (a\cdot b) \cdot c = a \cdot (b\cdot c)$

An identity element : $\exists ! e \in G, \forall a \in G, e \cdot a = a$ and $a\cdot e = a$ called the _identity element_

Inverse elements : $\forall a\in G, \exists ! b \in G, a\cdot b = e$ and $b \cdot a = e$. Such a $b$ is called the _inverse_ of $a$.

__Definition__ :
An $n$-dimensional _manifold_ $M$ is a topological (usually Hausdorff and second countable), _locally Euclidean space_. A topological space $M$ is said to be locally Euclidean of dimension $n$ if every point $p\in M$ has a neighborhood $U$ such that there is a homeomorphisms $\phi : U \rightarrow V$ for $V\subset \mathbb{R}^n$ an open subset.

__Definition__ : 
An $m$-dimensional manifold is a _smooth manifold_ embedded in $\mathbb{R}^n$ (provided $n \ge m$) if every point $p\in M$ is contained by $U\subseteq M$, defined by some function $\varphi : \mathbb{R}^n \rightarrow M, V \mapsto U$, where $V$ is an open subset of $\mathbb{R}^n$ which contains the origin. Additionally, $\varphi$ must be : 

A homeomorphism : $\varphi$ and $\varphi^{-1}$ are continuous.

Smooth : $\varphi \in C^{\infty}$

...





__Definition__ :
A _tangent vector_ at a point $p$ in a manifold $M$ is a _derivation_ at $p$.
The tangent vectors at $p$ form a vector space $T_pM$ called the _tangent space_ of $M$ at $p$.

__Definition__ :
A _Lie group_ is a subset $G$ of $\mathbb{R}^n$ such that $G$ is a group and a manifold in $\mathbb{R}^n$ and both the group operation $\cdot : G\rightarrow G$ and the inverse operation $.^{-1}: G\rightarrow G$ are smooth functions.

__Definition__ : _Lie bracket_

__Definition__ :
A _Lie algebra_ is an algebra $A$ togther with a _Lie bracket_ operator $\[.,.\] : A \times A \rightarrow A$.
An important fact is that the Lie algebra $A$ associated to a Lie group $M$ happens to be the tangent space
at the identity element $1$, that is: $A = T_1M$.

## SO(3) and SE(3) are Lie groups


## Riemannian Geometry
__Definition__ : 
Let $M$ be an $n$-dimensional manifold.
A _Riemannian metric_ on $M$ is a smoothly varying positive-definite inner-product $g_p$ on each tangent space $T_pM$.
We say that the pair $(M,g)$ is a _Riemannian manifold_.

We now defin what is the analogue of "straight lines" in a Riemannian manifold $(M,g)$ :

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
In the context of Lie groups, the exponential map is actually maps elements from the Lie algebra (tangent space) to the corresponding Lie group (manifold).
There is an "inverse" computation that we call the logarithm map which sends an element from the Lie group (manifold) to the Lie algebra (tangent space).