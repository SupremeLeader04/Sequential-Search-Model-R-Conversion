# Sequential-Search-Model-R-Conversion
A conversion of the MATLAB code for the GHK estimation of the homogenous Weitzmann Sequential Search Model from Ursu, Seiler, Honka (2024). The resulting code is in R, following an analogous structure and function.

While the MATLAB code replicates the estimation script for all 50 Monte Carlo seeds, R lacks a filename self-reference feature like mfilename, and so the code has been re-engineered such that the _same_ estimation script is ran 50 times, but with the seed passed as an argument. This achieves an identical result.

Note that exact replication of the MATLAB results is hindered by different RNG implementations between languages - both use a Mersenne-Twister algorithm but with different transformations and lower-level C implementations, so results will not be exactly the same. That said, the validity of the likelihood implementation has been verified using identical data and RNG inputs for various consumer examples.

Please note: much code is written by hand, other parts are LLM-assisted (especially the implementation of the third-party optimiser and the base likelihood function) but checked and edited for accuracy. For more precise knowledge of the authorship of each part of the code feel free to contact me on hayden.dyke23@imperial.ac.uk. If in doubt assume code is not mine.

---

## Structure

```
Sequential-Search-Model-R-Conversion
├── README.md # description of project
└── main # main container folder of R code
    ├── main.R # command script (runs estimation and averaging scripts over 50 seeds)
    ├── average.R # script to average results of estimation for each seed
    ├── cal_se_num_hessian.R # script to calculate standard errors of estimates
    ├── estWeitz_ghk_D100.R # main estimation script (D100 represents 100 error draws)
    ├── fminsearchcon.R # optimiser script (LLM copy of third-party MATLAB)
    ├── liklWeitz_ghk_1 # first part of likelihood function (prepares data and handles output)
    ├── liklWeitz_ghk_2 # second part of likelihood function (calculates likelihood values based on GHK specification)
    └── simWeitz.R # generates simulation data for 50 seeds
```
---
## How to Run

_Instructions to come post-upload._

---
