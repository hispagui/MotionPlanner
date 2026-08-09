# Lie Groups

### Definition :
A group is a set $G$ together with a binary operation on ⁠$G$⁠, here denoted "$\cdot$", that combines any two elements $a$ and $b$ of $G$ to form an element of $G$ denoted $a \cdot b$⁠, such that the following three requirements, known as group axioms, are satisfied:

Associativity : $\forall a,b \in G, (a\cdot b) \cdot c = a \cdot (b\cdot c)$

An identity element : $\exists ! e \in G, \forall a \in G, e \cdot a = a$ and $a\cdot e = a$ called the _identity element_

Inverse elements : $\forall a\in G, \exists ! b \in G, a\cdot b = e$ and $b \cdot a = e$. Such a $b$ is called the _inverse_ of $a$.

### Definition :
An $n$-dimensional _manifold_ $M$ is a topological (usually Hausdorff and second countable), _locally Euclidean space_. A topological space $M$ is said to be locally Euclidean of dimension $n$ if every point $p\in M$ has a neighborhood $U$ such that there is a homeomorphisms $\phi : U \rightarrow V$ for $V\subset \mathbb{R}^n$ an open subset.

### Definition : 
An $m$-dimensional manifold is a _smooth manifold_ embedded in $\mathbb{R}^n$ (provided $n \ge m$) if every point $p\in M$ is contained by $U\subseteq M$, defined by some function $\varphi : \mathbb{R}^n \rightarrow M, V \mapsto U$, where $V$ is an open subset of $\mathbb{R}^n$ which contains the origin. Additionally, $\varphi$ must be : 

A homeomorphism : $\varphi$ and $\varphi^{-1}$ are continuous.

Smooth : $\varphi \in C^{\infty}$

...



### Definition :
A _tangent vector_ at a point $p$ in a manifold $M$ is a _derivation_ at $p$.
The tangent vectors at $p$ form a vector space $T_pM$ called the _tangent space_ of $M$ at $p$.

__Definition__ :
A _Lie group_ is a subset $G$ of $\mathbb{R}^n$ such that $G$ is a group and a manifold in $\mathbb{R}^n$ and both the group operation $\cdot : G\rightarrow G$ and the inverse operation $.^{-1}: G\rightarrow G$ are smooth functions.


__Definition__ :
A _Lie algebra_ is an algebra $A$ togther with a _Lie bracket_ operator $\[.,.\] : A \times A \rightarrow A$