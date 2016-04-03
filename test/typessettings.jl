# This file is part of Kpax3. License is MIT.

import StatsBase: WeightVec
import Distributions: Beta

outfile = "../build/test.bin"
T = 1000000
burnin = 10000
tstep = 1
op = [0.7; 0.2; 0.1]
α = 0.0
θ = 1.0
γ = [0.6; 0.35; 0.05]
r = log(0.001) / log(0.95)
λs1 = 1.0
λs2 = 1.0
parawm = 5.0
maxclust = 500
maxunit = 500
verbose = true
verbosestep = 10000

settings = KSettings(outfile, T, burnin, tstep, op, α, θ, γ, r, λs1, λs2,
                     parawm, maxclust, maxunit, verbose, verbosestep)

@test settings.outfile == outfile
@test settings.T == T
@test settings.burnin == burnin
@test settings.tstep == tstep
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
@test settings.maxclust == maxclust
@test settings.maxunit == maxunit
@test settings.verbose == verbose
@test settings.verbosestep == verbosestep

@test_throws KDomainError KSettings(outfile, 0, burnin, tstep, op, α, θ, γ, r,
                                    λs1, λs2, parawm, maxclust, maxunit,
                                    verbose, verbosestep)
@test_throws KDomainError KSettings(outfile, T, -1, tstep, op, α, θ, γ, r,
                                    λs1, λs2, parawm, maxclust, maxunit,
                                    verbose, verbosestep)
@test_throws KDomainError KSettings(outfile, T, burnin, -1, op, α, θ, γ, r,
                                    λs1, λs2, parawm, maxclust, maxunit,
                                    verbose, verbosestep)
@test_throws KInputError KSettings(outfile, T, burnin, tstep, [1.0; 0.0], α,
                                    θ, γ, r, λs1, λs2, parawm, maxclust,
                                    maxunit, verbose, verbosestep)
@test_throws KDomainError KSettings(outfile, T, burnin, tstep, [1.0; 0.0; -1.0],
                                    α, θ, γ, r, λs1, λs2, parawm, maxclust,
                                    maxunit, verbose, verbosestep)
@test_throws KInputError KSettings(outfile, T, burnin, tstep, op, α, θ,
                                    [1.0; 0.0], r, λs1, λs2, parawm, maxclust,
                                    maxunit, verbose, verbosestep)
@test_throws KDomainError KSettings(outfile, T, burnin, tstep, op, α, θ,
                                    [1.0; 0.0; -1.0], r, λs1, λs2, parawm,
                                    maxclust, maxunit, verbose, verbosestep)
@test_throws KDomainError KSettings(outfile, T, burnin, tstep, op, α, θ, γ, 0.0,
                                    λs1, λs2, parawm, maxclust, maxunit,
                                    verbose, verbosestep)
@test_throws KDomainError KSettings(outfile, T, burnin, tstep, op, α, θ, γ, r,
                                    0.0, λs2, parawm, maxclust, maxunit,
                                    verbose, verbosestep)
@test_throws KDomainError KSettings(outfile, T, burnin, tstep, op, α, θ, γ, r,
                                    λs1, 0.0, parawm, maxclust, maxunit,
                                    verbose, verbosestep)
@test_throws KDomainError KSettings(outfile, T, burnin, tstep, op, α, θ, γ, r,
                                    λs1, λs2, 0.0, maxclust, maxunit,
                                    verbose, verbosestep)
@test_throws KDomainError KSettings(outfile, T, burnin, tstep, op, α, θ, γ, r,
                                    λs1, λs2, parawm, 0, maxunit,
                                    verbose, verbosestep)
@test_throws KDomainError KSettings(outfile, T, burnin, tstep, op, α, θ, γ, r,
                                    λs1, λs2, parawm, maxclust, 0,
                                    verbose, verbosestep)
@test_throws KDomainError KSettings(outfile, T, burnin, tstep, op, α, θ, γ, r,
                                    λs1, λs2, parawm, maxclust, maxunit,
                                    verbose, -1)
