# PencilVlasov.jl

Solving Vlasov-Poisson system in parallel in a 2D2V domain using [PencilArrays.jl](https://github.com/jipolanco/PencilArrays.jl)
and [PencilFFTs.jl](https://github.com/jipolanco/PencilFFTs.jl).

Does not work yet...


$$
\partial_t f + y \partial_x f +E \partial_y f = 0, f(t=0, x, y)= f_0(x, y), x\in [0, 4\pi], y\in \mathbb{R},
$$

where the electric field $E$ derives from a potential
$\phi(t, x)\in\mathbb{R}$ which satisfies a Poisson equation

$$
\partial_x^2 \phi = \int_{\mathbb{R}} f dy - 1.
$$

The initial condition is

```math
f_0(x, y, vx, vy)= \frac{1}{2\pi}e^{-(vx^2+vy^2)/2}(1 + 0.001\cos(2\pi kx x) \cos(2\pi ky y)).
```
