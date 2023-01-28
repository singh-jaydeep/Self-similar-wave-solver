# Self-similar-wave-solver
Numerical solver for a linear wave equation in similarity coordinates in C++

Solves the linear equation $ \partial_s \partial_z (r\phi) - (1-k) |z| \partial_z^2 (r\phi) + c \partial_z (r\phi) + V(z) (r\phi) =0, $ where $k >0$ is a small parameter, $c$ a damping parameter, and $V(z)$ a potential. The equation is defined on a domain $(s,z) \in $\mathbb{R}_+ \times [0,1]$ with Dirichlet boundary conditions at $z=0$. 

This equation governs spherically symmetric waves on a $3+1$-dimensional self-similar spacetime, after a transformation from standard double null gauge to similarity coordinates. $s$ becomes the time coordinate, $z$ the coordinate parameterizing integral curves of the self-similar vector field.
