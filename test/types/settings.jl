# This file is part of Kpax3. License is MIT.

ifile = "data/proper_aa.fasta"
ofile = "../build/test.bin"
protein = false
miss = zeros(UInt8, 0)
l = 10
α = 0.0
θ = 1.0
γ = [0.3; 0.3; 0.4]
r = log(0.01) / log(0.9)
maxclust = 2
maxunit = 2
verbose = false
verbosestep = 1
popsize = 10
maxiter = 5
maxgap = 1
xrate = 0.8
mrate = 0.2
T = 10
burnin = 100
tstep = 2
op = [0.5; 0.3; 0.2]
λs1 = 1.5
λs2 = 1.5
parawm = 2.0

settings = KSettings(ifile, ofile, protein=protein, miss=miss, l=l, α=α, θ=θ,
                     γ=γ, r=r, maxclust=maxclust, maxunit=maxunit,
                     verbose=verbose, verbosestep=verbosestep, popsize=popsize,
                     maxiter=maxiter, maxgap=maxgap, xrate=xrate, mrate=mrate,
                     T=T, burnin=burnin, tstep=tstep, op=op, λs1=λs1, λs2=λs2,
                     parawm=parawm)

@test settings.ifile == ifile
@test settings.ofile == ofile
@test settings.l == l
@test settings.α == α
@test settings.θ == θ
@test settings.γ == γ
@test settings.r == r
@test settings.maxclust == maxclust
@test settings.maxunit == maxunit
@test settings.verbose == verbose
@test settings.verbosestep == verbosestep
@test settings.popsize == popsize
@test settings.maxiter == maxiter
@test settings.maxgap == maxgap
@test settings.xrate == xrate
@test settings.mrate == mrate
@test settings.T == T
@test settings.burnin == burnin
@test settings.tstep == tstep
@test isa(settings.op, WeightVec)
@test values(settings.op) == op
@test isa(settings.distws, Beta)
@test settings.distws.α == λs1
@test settings.distws.β == λs2
@test settings.parawm == parawm

settings = KSettings(ifile, ofile, protein=true, miss=zeros(UInt8, 0))

@test settings.protein
@test settings.miss == UInt8['?', '*', '#', '-', 'b', 'j', 'x', 'z']

settings = KSettings(ifile, ofile, protein=true,
                     miss=UInt8['?', '*', '#', 'b', 'j', 'x', 'z'])

@test settings.protein
@test settings.miss == UInt8['?', '*', '#', 'b', 'j', 'x', 'z']

settings = KSettings(ifile, ofile, protein=false, miss=zeros(UInt8, 0))

@test !settings.protein
@test settings.miss == UInt8['?', '*', '#', '-', 'b', 'd', 'h', 'k', 'm', 'n',
                             'r', 's', 'v', 'w', 'x', 'y', 'j', 'z']

settings = KSettings(ifile, ofile, protein=false,
                     miss=UInt8['?', '*', '#', 'b', 'd', 'h', 'k', 'm', 'n',
                                'r', 's', 'v', 'w', 'x', 'y', 'j', 'z'])

@test !settings.protein
@test settings.miss == UInt8['?', '*', '#', 'b', 'd', 'h', 'k', 'm', 'n', 'r',
                             's', 'v', 'w', 'x', 'y', 'j', 'z']

@test_throws KDomainError KSettings(ifile, ofile, miss=UInt8[63; 0])
@test_throws KDomainError KSettings(ifile, ofile, l=-1)
@test_throws KInputError  KSettings(ifile, ofile, γ=[1.0; 0.0])
@test_throws KDomainError KSettings(ifile, ofile, γ=[1.0; 0.0; -1.0])
@test_throws KDomainError KSettings(ifile, ofile, r=0.0)
@test_throws KDomainError KSettings(ifile, ofile, r=-1.0)
@test_throws KDomainError KSettings(ifile, ofile, maxclust=0)
@test_throws KDomainError KSettings(ifile, ofile, maxunit=0)
@test_throws KDomainError KSettings(ifile, ofile, verbosestep=-1)
@test_throws KDomainError KSettings(ifile, ofile, popsize=-1)
@test_throws KDomainError KSettings(ifile, ofile, popsize=0)
@test_throws KDomainError KSettings(ifile, ofile, popsize=1)
@test_throws KDomainError KSettings(ifile, ofile, maxiter=0)
@test_throws KDomainError KSettings(ifile, ofile, maxgap=-1)
@test_throws KDomainError KSettings(ifile, ofile, xrate=-1.0)
@test_throws KDomainError KSettings(ifile, ofile, xrate=2.0)
@test_throws KDomainError KSettings(ifile, ofile, mrate=-1.0)
@test_throws KDomainError KSettings(ifile, ofile, mrate=2.0)
@test_throws KDomainError KSettings(ifile, ofile, T=0)
@test_throws KDomainError KSettings(ifile, ofile, burnin=-1)
@test_throws KDomainError KSettings(ifile, ofile, tstep=-1)
@test_throws KInputError  KSettings(ifile, ofile, op=[1.0; 0.0])
@test_throws KDomainError KSettings(ifile, ofile, op=[1.0; 0.0; -1.0])
@test_throws KDomainError KSettings(ifile, ofile, λs1=0.0)
@test_throws KDomainError KSettings(ifile, ofile, λs1=-1.0)
@test_throws KDomainError KSettings(ifile, ofile, λs2=0.0)
@test_throws KDomainError KSettings(ifile, ofile, λs2=-1.0)
@test_throws KDomainError KSettings(ifile, ofile, parawm=0.0)
@test_throws KDomainError KSettings(ifile, ofile, parawm=-1.0)
