# Rigid-body motions

In this project, we study the geometry and the topology of the configuration space : how orientations and poses live on curved manifolds,

how the free space around obstacles splits trajectories into distinct homotopy classes, and whether that topological structure actually helps a planner.

In this 



## Thread A - The topology of the configuration space itself (no obstacles)
    Here a configuration is an orientation $C = SO(3)$ or a full pose $C = SE(3)$.
    These spaces are curved and topologically non trivial, $SO(3) \cong \mathbb{RP}^3$ and $S^3$ sits above it as its universal cover, already $\pi_1(SO(3)) = \mathbb{Z}/2$.

    In this thread we represent interpolation between poses correctly *(see so3.py and se3.py)*

## Thread B - The topology of the free space (with obstacles)
    This is coming up next