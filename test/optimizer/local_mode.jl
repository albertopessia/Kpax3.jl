# This file is part of Kpax3. License is MIT.

settings = KSettings("../build/test.bin", T=1, burnin=0, tstep=1,
                     op=[0.6; 0.3; 0.1], α=0.0, θ=1.0, γ=[0.6; 0.35; 0.05],
                     r=135.0, λs1=1.0, λs2=1.0, parawm=5.0, maxclust=6,
                     maxunit=6, verbose=true, verbosestep=1)

data = UInt8[1 1 1 1 0 0;
             0 0 1 1 0 0;
             1 0 0 0 1 1;
             1 0 1 0 1 0]

m, n = size(data)

#=
R = [1; 1; 2; 2; 3; 3]
k = length(unique(R))

priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, k, settings.γ, settings.r)

state = AminoAcidState(data, R, priorR, priorC, settings)

C = zeros(UInt8, n, m)

cs = ((UInt8[1], UInt8[2], UInt8[3], UInt8[4]),
      (UInt8[1; 1], UInt8[2; 2], UInt8[3; 3], UInt8[3; 4], UInt8[4; 3],
       UInt8[4; 4]),
      (UInt8[1; 1; 1], UInt8[2; 2; 2], UInt8[3; 3; 3], UInt8[3; 3; 4],
       UInt8[3; 4; 3], UInt8[4; 3; 3], UInt8[3; 4; 4], UInt8[4; 3; 4],
       UInt8[4; 4; 3], UInt8[4; 4; 4]),
      (UInt8[1; 1; 1; 1], UInt8[2; 2; 2; 2], UInt8[3; 3; 3; 3],
       UInt8[3; 3; 3; 4], UInt8[3; 3; 4; 3], UInt8[3; 4; 3; 3],
       UInt8[4; 3; 3; 3], UInt8[3; 3; 4; 4], UInt8[3; 4; 3; 4],
       UInt8[3; 4; 4; 3], UInt8[4; 3; 3; 4], UInt8[4; 3; 4; 3],
       UInt8[4; 4; 3; 3], UInt8[3; 4; 4; 4], UInt8[4; 3; 4; 4],
       UInt8[4; 4; 3; 4], UInt8[4; 4; 4; 3], UInt8[4; 4; 4; 4]),
      (UInt8[1; 1; 1; 1; 1], UInt8[2; 2; 2; 2; 2], UInt8[3; 3; 3; 3; 3],
       UInt8[3; 3; 3; 3; 4], UInt8[3; 3; 3; 4; 3], UInt8[3; 3; 4; 3; 3],
       UInt8[3; 4; 3; 3; 3], UInt8[4; 3; 3; 3; 3], UInt8[3; 3; 3; 4; 4],
       UInt8[3; 3; 4; 3; 4], UInt8[3; 3; 4; 4; 3], UInt8[3; 4; 3; 3; 4],
       UInt8[3; 4; 3; 4; 3], UInt8[3; 4; 4; 3; 3], UInt8[4; 3; 3; 3; 4],
       UInt8[4; 3; 3; 4; 3], UInt8[4; 3; 4; 3; 3], UInt8[4; 4; 3; 3; 3],
       UInt8[3; 3; 4; 4; 4], UInt8[3; 4; 3; 4; 4], UInt8[3; 4; 4; 3; 4],
       UInt8[3; 4; 4; 4; 3], UInt8[4; 3; 3; 4; 4], UInt8[4; 3; 4; 3; 4],
       UInt8[4; 4; 3; 3; 4], UInt8[4; 3; 4; 4; 3], UInt8[4; 4; 3; 4; 3],
       UInt8[4; 4; 4; 3; 3], UInt8[3; 4; 4; 4; 4], UInt8[4; 3; 4; 4; 4],
       UInt8[4; 4; 3; 4; 4], UInt8[4; 4; 4; 3; 4], UInt8[4; 4; 4; 4; 3],
       UInt8[4; 4; 4; 4; 4]),
      (UInt8[1; 1; 1; 1; 1; 1], UInt8[2; 2; 2; 2; 2; 2],
       UInt8[3; 3; 3; 3; 3; 3], UInt8[3; 3; 3; 3; 3; 4],
       UInt8[3; 3; 3; 3; 4; 3], UInt8[3; 3; 3; 4; 3; 3],
       UInt8[3; 3; 4; 3; 3; 3], UInt8[3; 4; 3; 3; 3; 3],
       UInt8[4; 3; 3; 3; 3; 3], UInt8[3; 3; 3; 3; 4; 4],
       UInt8[3; 3; 3; 4; 3; 4], UInt8[3; 3; 3; 4; 4; 3],
       UInt8[3; 3; 4; 3; 3; 4], UInt8[3; 3; 4; 3; 4; 3],
       UInt8[3; 3; 4; 4; 3; 3], UInt8[3; 4; 3; 3; 3; 4],
       UInt8[3; 4; 3; 3; 4; 3], UInt8[3; 4; 3; 4; 3; 3],
       UInt8[3; 4; 4; 3; 3; 3], UInt8[4; 3; 3; 3; 3; 4],
       UInt8[4; 3; 3; 3; 4; 3], UInt8[4; 3; 3; 4; 3; 3],
       UInt8[4; 3; 4; 3; 3; 3], UInt8[4; 4; 3; 3; 3; 3],
       UInt8[3; 3; 3; 4; 4; 4], UInt8[3; 3; 4; 3; 4; 4],
       UInt8[3; 3; 4; 4; 3; 4], UInt8[3; 3; 4; 4; 4; 3],
       UInt8[3; 4; 3; 3; 4; 4], UInt8[3; 4; 3; 4; 3; 4],
       UInt8[3; 4; 4; 3; 3; 4], UInt8[3; 4; 3; 4; 4; 3],
       UInt8[3; 4; 4; 3; 4; 3], UInt8[3; 4; 4; 4; 3; 3],
       UInt8[4; 3; 3; 3; 4; 4], UInt8[4; 3; 3; 4; 3; 4],
       UInt8[4; 3; 4; 3; 3; 4], UInt8[4; 4; 3; 3; 3; 4],
       UInt8[4; 3; 3; 4; 4; 3], UInt8[4; 3; 4; 3; 4; 3],
       UInt8[4; 4; 3; 3; 4; 3], UInt8[4; 3; 4; 4; 3; 3],
       UInt8[4; 4; 3; 4; 3; 3], UInt8[4; 4; 4; 3; 3; 3],
       UInt8[3; 3; 4; 4; 4; 4], UInt8[3; 4; 3; 4; 4; 4],
       UInt8[3; 4; 4; 3; 4; 4], UInt8[3; 4; 4; 4; 3; 4],
       UInt8[3; 4; 4; 4; 4; 3], UInt8[4; 3; 3; 4; 4; 4],
       UInt8[4; 3; 4; 3; 4; 4], UInt8[4; 3; 4; 4; 3; 4],
       UInt8[4; 3; 4; 4; 4; 3], UInt8[4; 4; 3; 3; 4; 4],
       UInt8[4; 4; 3; 4; 3; 4], UInt8[4; 4; 3; 4; 4; 3],
       UInt8[4; 4; 4; 3; 3; 4], UInt8[4; 4; 4; 3; 4; 3],
       UInt8[4; 4; 4; 4; 3; 3], UInt8[3; 4; 4; 4; 4; 4],
       UInt8[4; 3; 4; 4; 4; 4], UInt8[4; 4; 3; 4; 4; 4],
       UInt8[4; 4; 4; 3; 4; 4], UInt8[4; 4; 4; 4; 3; 4],
       UInt8[4; 4; 4; 4; 4; 3], UInt8[4; 4; 4; 4; 4; 4]))

logpnew = -Inf
logpold = -Inf
l = 0
for c1 in cs[state.k], c2 in cs[state.k], c3 in cs[state.k], c4 in cs[state.k]
  l += 1

  tmp = hcat(c1, c2, c3, c4)

  for l in 1:state.k
    state.C[state.cl[l], :] = tmp[l, :]
  end

  logpnew = logcondpostC(state.C, state.cl, state.k, state.v, state.n1s,
                         priorC.logω, priorC)

  if logpnew > logpold
    logpold = logpnew
    C[1:state.k, :] = copy(state.C[state.cl[1:state.k], :])
  end
end

C
=#

# k = 1
R = [1; 1; 1; 1; 1; 1]
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, k, settings.γ, settings.r)
C = zeros(UInt8, 6, 4)
C[1, :] = 0x01
cl = [1; 0; 0; 0; 0; 0]
k = 1
v = [6; 0; 0; 0; 0; 0]
n1s = zeros(Float64, 6, 4)
n1s[1, :] = sum(float(data), 2)'
logω = copy(priorC.logω)

s = AminoAcidState(data, R, priorR, priorC, settings)
computelocalmode!(s.v, s.n1s, s.C, s.cl, s.k, s.logpC, priorC)

@test s.C == C
@test_approx_eq_eps s.logpC[1] logpriorC(C, cl, k, priorC.logγ, logω) ε
@test_approx_eq_eps s.logpC[2] logcondpostC(C, cl, k, v, n1s, logω, priorC) ε

t = AminoAcidState(data, R, priorR, priorC, settings)
computelocalmode!(t, priorC)

@test t.C == C
@test_approx_eq_eps t.logpC[1] logpriorC(C, cl, k, priorC.logγ, logω) ε
@test_approx_eq_eps t.logpC[2] logcondpostC(C, cl, k, v, n1s, logω, priorC) ε

# k = 2
R = [1; 1; 1; 1; 2; 2]
k = 2
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, k, settings.γ, settings.r)
C = zeros(UInt8, 6, 4)
C[1:2, :] = UInt8[4 1 2 1; 3 1 2 1]
cl = [1; 2; 0; 0; 0; 0]
v = [4; 2; 0; 0; 0; 0]
n1s = zeros(Float64, 6, 4)
n1s[1:2, :] = hcat(sum(float(data[:, R .== 1]), 2),
                   sum(float(data[:, R .== 2]), 2))'
logω = copy(priorC.logω)

s = AminoAcidState(data, R, priorR, priorC, settings)
computelocalmode!(s.v, s.n1s, s.C, s.cl, s.k, s.logpC, priorC)

@test s.C == C
@test_approx_eq_eps s.logpC[1] logpriorC(C, cl, k, priorC.logγ, logω) ε
@test_approx_eq_eps s.logpC[2] logcondpostC(C, cl, k, v, n1s, logω, priorC) ε

t = AminoAcidState(data, R, priorR, priorC, settings)
computelocalmode!(t, priorC)

@test t.C == C
@test_approx_eq_eps t.logpC[1] logpriorC(C, cl, k, priorC.logγ, logω) ε
@test_approx_eq_eps t.logpC[2] logcondpostC(C, cl, k, v, n1s, logω, priorC) ε

# k = 3
R = [1; 1; 2; 2; 3; 3]
k = 3
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, k, settings.γ, settings.r)
C = zeros(UInt8, 6, 4)
C[1:3, :] = UInt8[4 3 2 1; 4 4 2 1; 3 3 2 1]
cl = [1; 2; 3; 0; 0; 0]
v = [2; 2; 2; 0; 0; 0]
n1s = zeros(Float64, 6, 4)
n1s[1:3, :] = hcat(sum(float(data[:, R .== 1]), 2),
                   sum(float(data[:, R .== 2]), 2),
                   sum(float(data[:, R .== 3]), 2))'
logω = copy(priorC.logω)

s = AminoAcidState(data, R, priorR, priorC, settings)
computelocalmode!(s.v, s.n1s, s.C, s.cl, s.k, s.logpC, priorC)

@test s.C == C
@test_approx_eq_eps s.logpC[1] logpriorC(C, cl, k, priorC.logγ, logω) ε
@test_approx_eq_eps s.logpC[2] logcondpostC(C, cl, k, v, n1s, logω, priorC) ε

t = AminoAcidState(data, R, priorR, priorC, settings)
computelocalmode!(t, priorC)

@test t.C == C
@test_approx_eq_eps t.logpC[1] logpriorC(C, cl, k, priorC.logγ, logω) ε
@test_approx_eq_eps t.logpC[2] logcondpostC(C, cl, k, v, n1s, logω, priorC) ε

# k = 4
R = [1; 2; 2; 2; 3; 4]
k = 4
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, k, settings.γ, settings.r)
C = zeros(UInt8, 6, 4)
C[1:4, :] = UInt8[2 3 2 1; 2 4 2 1; 2 3 2 1; 2 3 2 1]
cl = [1; 2; 3; 4; 0; 0]
v = [2; 2; 1; 1; 0; 0]
n1s = zeros(Float64, 6, 4)
n1s[1:4, :] = hcat(sum(float(data[:, R .== 1]), 2),
                   sum(float(data[:, R .== 2]), 2),
                   sum(float(data[:, R .== 3]), 2),
                   sum(float(data[:, R .== 4]), 2))'
logω = copy(priorC.logω)

s = AminoAcidState(data, R, priorR, priorC, settings)
computelocalmode!(s.v, s.n1s, s.C, s.cl, s.k, s.logpC, priorC)

@test s.C == C
@test_approx_eq_eps s.logpC[1] logpriorC(C, cl, k, priorC.logγ, logω) ε
@test_approx_eq_eps s.logpC[2] logcondpostC(C, cl, k, v, n1s, logω, priorC) ε

t = AminoAcidState(data, R, priorR, priorC, settings)
computelocalmode!(t, priorC)

@test t.C == C
@test_approx_eq_eps t.logpC[1] logpriorC(C, cl, k, priorC.logγ, logω) ε
@test_approx_eq_eps t.logpC[2] logcondpostC(C, cl, k, v, n1s, logω, priorC) ε

# k = 5
R = [1; 2; 3; 3; 4; 5]
k = 5
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, k, settings.γ, settings.r)
C = zeros(UInt8, 6, 4)
C[1:5, :] = UInt8[2 3 2 1; 2 3 2 1; 2 4 2 1; 2 3 2 1; 2 3 2 1]
cl = [1; 2; 3; 4; 5; 0]
v = [1; 1; 2; 1; 1; 0]
n1s = zeros(Float64, 6, 4)
n1s[1:5, :] = hcat(sum(float(data[:, R .== 1]), 2),
                   sum(float(data[:, R .== 2]), 2),
                   sum(float(data[:, R .== 3]), 2),
                   sum(float(data[:, R .== 4]), 2),
                   sum(float(data[:, R .== 5]), 2))'
logω = copy(priorC.logω)

s = AminoAcidState(data, R, priorR, priorC, settings)
computelocalmode!(s.v, s.n1s, s.C, s.cl, s.k, s.logpC, priorC)

@test s.C == C
@test_approx_eq_eps s.logpC[1] logpriorC(C, cl, k, priorC.logγ, logω) ε
@test_approx_eq_eps s.logpC[2] logcondpostC(C, cl, k, v, n1s, logω, priorC) ε

t = AminoAcidState(data, R, priorR, priorC, settings)
computelocalmode!(t, priorC)

@test t.C == C
@test_approx_eq_eps t.logpC[1] logpriorC(C, cl, k, priorC.logγ, logω) ε
@test_approx_eq_eps t.logpC[2] logcondpostC(C, cl, k, v, n1s, logω, priorC) ε

# k = 6
R = [1; 2; 3; 4; 5; 6]
k = 6
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, k, settings.γ, settings.r)
C = UInt8[1 1 2 2; 1 1 2 2 ; 1 1 2 2; 1 1 2 2; 1 1 2 2; 1 1 2 2]
cl = [1; 2; 3; 4; 5; 6]
v = [1; 1; 1; 1; 1; 1]
n1s = float(data)'
logω = copy(priorC.logω)

s = AminoAcidState(data, R, priorR, priorC, settings)
computelocalmode!(s.v, s.n1s, s.C, s.cl, s.k, s.logpC, priorC)

@test s.C == C
@test_approx_eq_eps s.logpC[1] logpriorC(C, cl, k, priorC.logγ, logω) ε
@test_approx_eq_eps s.logpC[2] logcondpostC(C, cl, k, v, n1s, logω, priorC) ε

t = AminoAcidState(data, R, priorR, priorC, settings)
computelocalmode!(t, priorC)

@test t.C == C
@test_approx_eq_eps t.logpC[1] logpriorC(C, cl, k, priorC.logγ, logω) ε
@test_approx_eq_eps t.logpC[2] logcondpostC(C, cl, k, v, n1s, logω, priorC) ε
