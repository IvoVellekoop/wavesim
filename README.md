# Wavesim - A fast and accurate method for solving the Helmholtz and time-independent Maxwell's equation

When using this code, please refer to:

G. Osnabrugge, S. Leedumrongwatthanakun, I.M. Vellekoop - A convergent Born series for solving the inhomogeneous Helmholtz equation in arbitrarily large media, Journal of Computational Physics Volume 322, 1 October 2016, Pages 113–124, doi:10.1016/j.jcp.2016.06.034
Freely available at: http://www.sciencedirect.com/science/article/pii/S0021999116302595

G. Osnabrugge, M. Benedictus, I.M. Vellekoop - An ultra-thin boundary layer for high-accuracy simulations, Optics Express 29 (2), 11 January 2021, Pages 1649-1658, doi:10.1364/OE.412833

## wavesim.org

We are working to improve and accelerate Wavesim further. On 27 February 2024, we released a new version that uses a CUDA-based acceleration module (cumex) to provide around 2x speed up. Want to find out more? Want to participate in the forum for discussions, queries, and requests? Then please visit www.wavesim.org.

## wavesim_py

On 3 October 2024, we released a [Python version of wavesim](https://github.com/IvoVellekoop/wavesim_py) on GitHub. This is a Python implementation of the Modified Born Series (MBS) approach for solving the Helmholtz equation in arbitrarily large media through [domain decomposition](https://arxiv.org/abs/2410.02395). With this new framework, we simulated a complex 3D structure of a remarkable 315×315×315 wavelengths (3.1⋅10^7) in size in just 379 seconds by solving over two GPUs. This represents a factor of 1.93 increase over the largest possible simulation on a single GPU without domain decomposition.
