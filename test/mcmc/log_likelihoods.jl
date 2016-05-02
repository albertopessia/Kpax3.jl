# This file is part of Kpax3. License is MIT.

settings = KSettings("../build/test.bin", T=1, burnin=0, tstep=1,
                     op=[0.6; 0.3; 0.1], α=0.0, θ=1.0, γ=[0.6; 0.35; 0.05],
                     r=135.0, λs1=1.0, λs2=1.0, parawm=5.0, maxclust=6,
                     maxunit=6, verbose=true, verbosestep=1)

data = UInt8[0 0 0 0 0 1;
             1 1 1 1 1 0;
             0 0 1 0 1 1;
             1 1 0 1 0 0;
             1 1 0 0 0 0;
             0 0 0 1 1 0;
             1 1 1 0 0 0;
             0 0 0 1 1 1;
             0 0 1 0 0 0;
             1 0 0 1 0 1;
             0 1 0 0 1 0;
             0 0 0 0 0 1;
             1 1 1 0 0 0;
             0 0 0 1 1 0;
             1 1 0 0 1 1;
             0 0 1 1 0 0;
             1 1 0 1 0 0;
             0 0 1 0 1 1]

m, n = size(data)

priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, settings.γ, settings.r)

# merge
# [1; 1; 1; 2; 2; 3] => [1; 1; 1; 3; 3; 3]
state = AminoAcidState(data, [1; 1; 1; 2; 2; 3], priorR, priorC, settings)

mergesupport = KSupport(m, n, settings.maxclust, settings.maxunit)
mergesupport.C[1:2, :] = UInt8[3 1 2 1 1 1 4 3 2 1 2 1 1 1 1 1 2 1;
                               4 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1]
mergesupport.cl[1:2] = [1; 2]
mergesupport.k = 2

R = [1; 1; 1; 3; 3; 3]
hi = 3
hj = 2
n1s = hcat(sum(float(data[:, R .== 1]), 2), sum(float(data[:, R .== 3]), 2))'
v = [3; 3]

loglikmerge!(hi, hj, vec(n1s[2, :]), v[2], priorC, mergesupport, state)

correctloglik = 0.0
lidx = 0
for b in 1:m
  lidx = mergesupport.C[1, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[1, b], v[1], priorC.A[lidx], priorC.B[lidx])
  lidx = mergesupport.C[2, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[2, b], v[2], priorC.A[lidx], priorC.B[lidx])
end

@test_approx_eq_eps mergesupport.loglik correctloglik ε

# split
# [1; 1; 1; 3; 3; 3] => [1; 1; 1; 3; 3; 2]
state = AminoAcidState(data, [1; 1; 1; 3; 3; 3], priorR, priorC, settings)

splitsupport = KSupport(m, n, settings.maxclust, settings.maxunit)
splitsupport.C[1:3, :] = UInt8[3 1 2 1 1 1 4 3 2 1 2 1 1 1 1 1 2 1;
                               4 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1;
                               3 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1]
splitsupport.cl[1:3] = [1; 2; 3]
splitsupport.k = 3

R = [1; 1; 1; 3; 3; 2]
hi = 3
n1s = hcat(sum(float(data[:, R .== 1]), 2), sum(float(data[:, R .== 3]), 2),
           sum(float(data[:, R .== 2]), 2))'
v = [3; 2; 1]

splitsupport.ni = vec(n1s[3, :])
splitsupport.vi = v[3]
splitsupport.nj = vec(n1s[2, :])
splitsupport.vj = v[2]

logliksplit!(hi, priorC, splitsupport, state)

correctloglik = 0.0
lidx = 0
for b in 1:m
  lidx = splitsupport.C[1, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[1, b], v[1], priorC.A[lidx], priorC.B[lidx])
  lidx = splitsupport.C[2, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[3, b], v[3], priorC.A[lidx], priorC.B[lidx])
  lidx = splitsupport.C[3, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[2, b], v[2], priorC.A[lidx], priorC.B[lidx])
end

@test_approx_eq_eps splitsupport.loglik correctloglik ε

# biased random walk
# [1; 1; 1; 3; 3; 2] => [1; 1; 1; 3; 3; 3]
state = AminoAcidState(data, [1; 1; 1; 3; 3; 2], priorR, priorC, settings)

brwsupport = KSupport(m, n, settings.maxclust, settings.maxunit)
brwsupport.C[1:2, :] = UInt8[3 1 2 1 1 1 4 3 2 1 2 1 1 1 1 1 2 1;
                             4 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1]
brwsupport.cl[1] = 1
brwsupport.k = 2

R = [1; 1; 1; 3; 3; 3]
k = 2
hi = 2
hj = 3
n1s = hcat(sum(float(data[:, R .== 1]), 2),
           sum(float(data[:, R .== 3]), 2))'
v = [3; 3]

brwsupport.ni = vec(data[:, 6])

loglikbrw!(k, hi, hj, priorC, brwsupport, state)

correctloglik = 0.0
lidx = 0
for b in 1:m
  lidx = brwsupport.C[1, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[1, b], v[1], priorC.A[lidx], priorC.B[lidx])
  lidx = brwsupport.C[2, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[2, b], v[2], priorC.A[lidx], priorC.B[lidx])
end

@test_approx_eq_eps brwsupport.loglik correctloglik ε

# [1; 1; 1; 3; 3; 2] => [1; 1; 2; 3; 3; 2]
state = AminoAcidState(data, [1; 1; 1; 3; 3; 2], priorR, priorC, settings)

brwsupport = KSupport(m, n, settings.maxclust, settings.maxunit)
brwsupport.C[1:3, :] = UInt8[3 1 2 1 1 1 4 3 2 1 2 1 1 1 1 1 2 1;
                             4 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1;
                             3 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1]
brwsupport.cl[1] = 3
brwsupport.k = 3

R = [1; 1; 2; 3; 3; 2]
k = 3
hi = 1
hj = 2
n1s = hcat(sum(float(data[:, R .== 1]), 2), sum(float(data[:, R .== 2]), 2),
           sum(float(data[:, R .== 3]), 2))'
v = [2; 2; 2]

brwsupport.ni = vec(data[:, 3])

loglikbrw!(k, hi, hj, priorC, brwsupport, state)

correctloglik = 0.0
lidx = 0
for b in 1:m
  lidx = brwsupport.C[1, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[3, b], v[3], priorC.A[lidx], priorC.B[lidx])
  lidx = brwsupport.C[2, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[1, b], v[1], priorC.A[lidx], priorC.B[lidx])
  lidx = brwsupport.C[3, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[2, b], v[2], priorC.A[lidx], priorC.B[lidx])
end

@test_approx_eq_eps brwsupport.loglik correctloglik ε

# [3; 3; 3; 2; 2; 4] => [3; 3; 1; 2; 2; 4]
state = AminoAcidState(data, [2; 2; 2; 1; 1; 3], priorR, priorC, settings)

state.R = [3; 3; 3; 2; 2; 4]
state.C[2:4, :] = copy(state.C[1:3, :])

state.emptycluster[1] = true
state.emptycluster[2:4] = false
state.cl[1:3] = [2; 3; 4]
state.k = 3

state.v[2:4] = copy(state.v[1:3])
state.n1s[2:4, :] = copy(state.n1s[1:3, :])
state.unit[4] = copy(state.unit[3])
state.unit[3] = copy(state.unit[2])
state.unit[2] = copy(state.unit[1])

brwsupport = KSupport(m, n, settings.maxclust, settings.maxunit)
brwsupport.C[1:4, :] = UInt8[3 1 2 1 1 1 4 3 2 1 2 1 1 1 1 1 2 1;
                             4 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1;
                             3 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1;
                             3 1 2 1 1 1 4 3 2 1 2 1 1 1 1 1 2 1]
brwsupport.cl[1:2] = [2; 4]
brwsupport.k = 4

R = [3; 3; 1; 2; 2; 4]
k = 4
hi = 3
hj = 1
n1s = hcat(sum(float(data[:, R .== 1]), 2),
           sum(float(data[:, R .== 2]), 2),
           sum(float(data[:, R .== 3]), 2),
           sum(float(data[:, R .== 4]), 2))'
v = [1; 2; 2; 1]

brwsupport.ni = vec(data[:, 3])

loglikbrw!(k, hi, hj, priorC, brwsupport, state)

correctloglik = 0.0
lidx = 0
for b in 1:m
  lidx = brwsupport.C[1, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[2, b], v[2], priorC.A[lidx], priorC.B[lidx])
  lidx = brwsupport.C[2, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[4, b], v[4], priorC.A[lidx], priorC.B[lidx])
  lidx = brwsupport.C[3, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[3, b], v[3], priorC.A[lidx], priorC.B[lidx])
  lidx = brwsupport.C[4, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[1, b], v[1], priorC.A[lidx], priorC.B[lidx])
end

@test_approx_eq_eps brwsupport.loglik correctloglik ε
