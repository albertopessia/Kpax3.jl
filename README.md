# Kpax3 - Bayesian cluster analysis of categorical data
## About
Kpax3 is a [Julia](http://julialang.org/) package for inferring the population structure of proteins. Output consists of a clustering of both the rows (statistical units) and columns (statistical variables) of the provided data matrix. It is an improved version of [K-Pax2](https://github.com/albertopessia/kpax2/), providing an MCMC algorithm for a proper Bayesian approach and a genetic algorithm for MAP estimation.

## Reference
To know more about the underlying statistical model, refer to the following publications (the first is the primary citation if you use the package):

Pessia, A. and Corander, J. (2017). "Bayesian cluster analysis of categorical data with supervised feature selection". _Submitted_

Pessia, A., Grad, Y., Cobey, S., Puranen, J. S., and Corander, J. (2015). "K-Pax2: Bayesian identification of cluster-defining amino acid positions in large sequence datasets". _Microbial Genomics_ **1**(1). [http://dx.doi.org/10.1099/mgen.0.000025](http://dx.doi.org/10.1099/mgen.0.000025)

## Installation
Kpax3 can be installed from within julia by typing

```julia
Pkg.clone("https://github.com/albertopessia/Kpax3.jl.git")
```

## Usage
The best way to learn how to use Kpax3 is by following the instructions in the [tutorial](tutorial/Kpax3_tutorial.jl). Full documentation of Kpax3 functions will be available in the near future.

## License
See [LICENSE.md](LICENSE.md)
