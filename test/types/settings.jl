# This file is part of Kpax3. License is MIT.

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

settings = KSettings(outfile, T=T, burnin=burnin, tstep=tstep, op=op, α=α,
                     θ=θ, γ=γ, r=r, λs1=λs1, λs2=λs2, parawm=parawm,
                     maxclust=maxclust, maxunit=maxunit, verbose=verbose,
                     verbosestep=verbosestep)

@test settings.fpath == outfile
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

@test_throws KDomainError KSettings(outfile, T=0)
@test_throws KDomainError KSettings(outfile, burnin=-1)
@test_throws KDomainError KSettings(outfile, tstep=-1)
@test_throws KInputError  KSettings(outfile, op=[1.0; 0.0])
@test_throws KDomainError KSettings(outfile, op=[1.0; 0.0; -1.0])
@test_throws KInputError  KSettings(outfile, γ=[1.0; 0.0])
@test_throws KDomainError KSettings(outfile, γ=[1.0; 0.0; -1.0])
@test_throws KDomainError KSettings(outfile, r=0.0)
@test_throws KDomainError KSettings(outfile, r=-1.0)
@test_throws KDomainError KSettings(outfile, λs1=0.0)
@test_throws KDomainError KSettings(outfile, λs1=-1.0)
@test_throws KDomainError KSettings(outfile, λs2=0.0)
@test_throws KDomainError KSettings(outfile, λs2=-1.0)
@test_throws KDomainError KSettings(outfile, parawm=0.0)
@test_throws KDomainError KSettings(outfile, parawm=-1.0)
@test_throws KDomainError KSettings(outfile, maxclust=0)
@test_throws KDomainError KSettings(outfile, maxunit=0)
@test_throws KDomainError KSettings(outfile, verbosestep=-1)
