# Sequential-Search-Model-R-Conversion
A conversion of the MATLAB code for the GHK estimation of the homogenous Weitzmann Sequential Search Model from Ursu, Seiler, Honka (2024). The resulting code is in R, following an analogous structure and function.

While the MATLAB code replicates the estimation script for all 50 Monte Carlo seeds, R lacks a filename self-reference feature like mfilename, and so the code has been re-engineered such that the _same_ estimation script is ran 50 times, but with the seed passed as an argument. This achieves an identical result.

Note that exact replication of the MATLAB results is hindered by different RNG implementations between languages - both use a Mersenne-Twister algorithm but with different transformations and lower-level C implementations, so results will not be exactly the same. That said, the validity of the likelihood implementation has been verified using identical data and RNG inputs for various consumer examples.

---


