# This file is part of Kpax3. License is MIT.

# the following code has been modified from
# http://thirld.com/blog/2015/05/30/julia-profiling-cheat-sheet/

# do not forget to run julia with --track-allocation=user
function benchmark()
  include("src/boot.jl")

  fastafile = "test/data/mcmc_6.fasta"
  partition = "test/data/mcmc_6.csv"
  outfile = "build/benchmark.bin"
  T = 1000000
  burnin = 500000
  tstep = 1
  op = [0.6; 0.3; 0.1]
  α = 0.5
  θ = -0.25
  γ = [0.6; 0.35; 0.05]
  r = log(0.001) / log(0.95)
  λs1 = 1.0
  λs2 = 1.0
  parawm = 5.0
  maxclust = 1
  maxunit = 1
  verbose = false
  verbosestep = 100000

  x = AminoAcidData(fastafile)

  # Run once, to force compilation.
  println("-- First run --")
  srand(20150326)
  @time kpax3aa(x, partition, outfile, T, burnin=burnin, tstep=tstep, op=op,
                α=α, θ=θ, γ=γ, r=r, λs1=λs1, λs2=λs2, parawm=parawm,
                maxclust=maxclust, maxunit=maxunit, verbose=verbose,
                verbosestep=verbosestep)

  # Run a second time, with profiling.
  println("-- Second run --")
  srand(20150326)
  Profile.init(delay=0.01)
  Profile.clear()
  Profile.clear_malloc_data()
  @profile @time kpax3aa(x, partition, outfile, T, burnin=burnin, tstep=tstep,
                         op=op, α=α, θ=θ, γ=γ, r=r, λs1=λs1, λs2=λs2,
                         parawm=parawm, maxclust=maxclust, maxunit=maxunit,
                         verbose=verbose, verbosestep=verbosestep)

  # Write profile results
  r = Profile.retrieve()
  f = open("build/profile.bin", "w")
  serialize(f, r)
  close(f)

  nothing
end
