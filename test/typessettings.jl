# This file is part of Kpax3. License is MIT.

import StatsBase: WeightVec
import Distributions: Beta

T = 1000000
outfile = "outfile.bin"
burnin = 10000
t = 1
op = [1.0, 0.0, 0.0]
α = 0.0
θ = 1.0
γ = [0.6, 0.35, 0.05]
r = log(0.001) / log(0.95)
λs1 = 1.0
λs2 = 1.0
parawm = 5.0
maxunit = 500
maxclust = 500
verbose = true
verbosestep = 10000

settings = KSettings(outfile, T, burnin, t, op, α, θ, γ, r, λs1, λs2, parawm,
                     maxclust, maxunit, verbose, verbosestep)

@test settings.outfile == outfile
@test settings.T == T
@test settings.burnin == burnin
@test settings.t == t
@test isa(settings.op, WeightVec)
@test values(settings.op) == op
@test settings.α == α
@test settings.θ == θ
@test settings.γ == γ
@test settings.r == r
@test isa(settings.distws, Beta)
@test settings.distws.α == λs1
@test settings.distws.β == λs2
@test settings.parawm == parawm
@test settings.maxunit == maxunit
@test settings.maxclust == maxclust
@test settings.verbose == verbose
@test settings.verbosestep == verbosestep
