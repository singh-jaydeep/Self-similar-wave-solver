# Self-similar-wave-solver
Numerical solver for a damped, linear wave equation on spherically symmetric, self-similar spacetimes. 

Solves the linear equation 
```math
\partial_s \partial_z (r\phi) - (1-k^2) |z| \partial_z^2 (r\phi) + c \partial_z (r\phi) + V(z) (r\phi) =0,
```
where $k >0$ is a small parameter, $c$ a damping parameter, and $V(z)$ a potential. The equation is defined on a domain 
```math
(s,z) \in [0,\infty) \times [-1,0]
```
with Dirichlet boundary conditions on $\{z=-1\}$, representing the center of symmetry. 

This equation governs spherically symmetric waves on a $3+1$-dimensional self-similar spacetime, after a transformation from standard double-null gauge to similarity coordinates via 
```math
(s,z) = (-\log|u|, v|u|^{-1+k^2}).
```
The coordinate $z$ parameterizes integral curves of the self-similar vector field, and $s$ is the natural time coordinate. If $k = c = 0$ and $V(z) = 0$, recover the wave equation on a subset of Minkowski spacetime. If $k \in (0,1)$, $c = 0$, and $V(z)$ is chosen appropriately, we model the $k$-self-similar spacetimes introduced in [this paper](https://www.jstor.org/stable/2118619).

