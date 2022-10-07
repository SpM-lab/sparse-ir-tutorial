# Discrete Lehmann representation

We explain the implementation of Discrete Lehmann Representation (DLR){cite:p}`DLR` in sparse-ir.
For both of fermions and bosons, we expand the Green's function as

$$
G(\tau) = - \sum_{p=1}^L K^\mathrm{L}(\tau, \bar{\omega}_p) c_p
$$

for $0 < \tau < \beta$.
The pole positions $\{\bar{\omega}_1, \cdots, \bar{\omega}_{L}\}$ are chosen as the extrema of $V'_{L-1}(\omega)$.
$\{K^\mathrm{L}(\tau, \bar{\omega}_p) \}$ forms a non-orthogonal basis set in $\tau$, which is common for fermions and bosons.


```{note}
The poles on the real-frequency axis selected for the DLR are based on a rank-revealing
decomposition, which offers guaranteed accuracy. Here, we instead select the pole locations
based on the zeros of the IR basis functions on the real axis {cite:p}`scipostreview`, which is a heuristic.  We do not expect that difference to matter, but please don't blame the DLR authors if we were wrong :-)
```

## Fermions
For fermions, this is equivalent to modeling the spectral function as

$$
    A(\omega) = \rho(\omega) = \sum_{p=1}^L c_p \delta(\omega - \bar{\omega}_p).
$$

The choice of the pole positions is heuristic but allows a numerically stable transform between $\rho_l$ and $c_p$ through the relation

$$
\rho_l = \sum_{p=1}^L \boldsymbol{V}_{lp} c_p,
$$

where the matrix $\boldsymbol{V}_{lp}~[\equiv V_l(\bar{\omega}_p)]$ is well-conditioned.
The Matsubara Green's function is expanded as

$$
G(\mathrm{i}\nu) = \int \mathrm{d}\omega \frac{A(\omega)}{\mathrm{i}\nu - \omega} = \sum_{p=1}^L \frac{c_p}{\mathrm{i}\nu - \bar{\omega}_p}.
$$


## Bosons

In the Matsubara-frequency space, The DLR for bosons is defined as follows:

$$
    A(\omega) = \sum_{p=1}^L c_p \tanh(\beta\bar{\omega}_p/2) \delta(\omega - \bar{\omega}_p),
$$

$$
    \rho(\omega) = \frac{A(\omega)}{\tanh(\beta\omega/2)} = \sum_{p=1}^L c_p \delta(\omega - \bar{\omega}_p),
$$

$$
G(\mathrm{i}\nu) = \int \mathrm{d}\omega \frac{A(\omega)}{\mathrm{i}\nu - \omega} = \sum_{p=1}^L \frac{c_p \tanh(\beta\bar{\omega}_p/2)}{\mathrm{i}\nu - \bar{\omega}_p}.
$$

When transforming data from IR to DLR, one can use the same transformation matrix for fermions and bosons.