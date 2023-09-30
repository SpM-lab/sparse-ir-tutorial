# Summation over Matsubara axis

In many cases, we want to perform the summation of a Greens-function-like object $A(\mathrm{i}\omega)$ over the Matsubara axis.

The Fourier transform of $A$ reads

$$
A(\tau) = \frac{1}{\beta} \sum_{\omega} A(\mathrm{i}\omega) e^{-\mathrm{i}\omega \tau}.
$$

This leads to the following the two formula:

$$
\begin{align}
 \sum_{\omega} A(\mathrm{i}\omega) e^{\mathrm{i}\omega 0^+} &= \beta A(\tau=0^-), \\
 \sum_{\omega} A(\mathrm{i}\omega) e^{\mathrm{i}\omega 0^-} &= \beta A(\tau=0^+). \\
\end{align}
$$

We now expand $A(\mathrm{i}\omega)$ at high frequencies as

$$
A(\mathrm{i}\omega) = \frac{c_1}{\mathrm{i}\omega} + \frac{c_2}{(\mathrm{i}\omega)^2} + \cdots.
$$

As discussed in [Sec. B3 of E. Gull's Ph. D thesis](https://www.research-collection.ethz.ch/handle/20.500.11850/104013),
$A(\tau=0^+) = A(\tau=0^-)$ if and only if $c_1 = 0$.
This condition is equivalent that $A(\mathrm{i}\omega)$ vanishes at high frequencies faster than $O(1/{\mathrm{i}\omega})$.
If $c_1 \neq 0$, the summation does NOT converge without a convergence factor and thus $ \sum_{\omega} A(\mathrm{i}\omega) e^{\mathrm{i}\omega 0^+} \neq  \sum_{\omega} A(\mathrm{i}\omega) e^{\mathrm{i}\omega 0^-}$.