# This file is part of Kpax3. License is MIT.

ifile = "data/proper_aa.fasta"
ofile = "../build/test.bin"

settings = KSettings(ifile, ofile, maxclust=15, maxunit=1, op=[0.0; 1.0; 0.0])

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

(m, n) = size(data)

priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, settings.γ, settings.r)

# [2; 2; 2; 1; 1; 3] => [2; 2; 2; 1; 1; 1]
support = KSupport(m, n, settings.maxclust, settings.maxunit)
state = AminoAcidState(data, [2; 2; 2; 1; 1; 3], priorR, priorC, settings)

k = 2
i = 6
hi = state.R[i]
hj = state.R[5]
cl = [1; 2]
v = [3; 3]
n1s = hcat(sum(float(data)[:, 1:3], 2), sum(float(data)[:, 4:6], 2))'

initsupportbrw!(k, i, state.v[hi], data, priorC, settings, support)

support.lograR = logratiopriorrowbrwmerge(k, state.v[hj], priorR)

support.C[1:2, :] = UInt8[4 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1;
                          3 1 2 1 1 1 4 3 2 1 2 1 1 1 1 1 2 1]
support.cl[1:2] = [2; 1]
support.k = 2

support.logpC = [logpriorC(support.C, cl, k, priorC);
                 logcondpostC(support.C, cl, k, v, n1s, priorC)]

loglikbrw!(k, hi, hj, priorC, support, state)

performbrw!(i, hi, hj, k, priorC, settings, support, state)

@test state.R == [2; 2; 2; 1; 1; 1]
@test state.k == 2

@test state.C[1, :] == support.C[2, :]
@test !state.emptycluster[1]
@test state.cl[1] == 1
@test state.v[1] == 3
@test state.n1s[1, :] == n1s[2, :]
@test state.unit[1][1:3] == [4; 5; 6]

@test state.C[2, :] == support.C[1, :]
@test !state.emptycluster[2]
@test state.cl[2] == 2
@test state.v[2] == 3
@test state.n1s[2, :] == n1s[1, :]
@test state.unit[2][1:3] == [1; 2; 3]

@test state.emptycluster[3]

@test_approx_eq_eps state.logpR logdPriorRow(n, k, v, priorR) ε

@test_approx_eq_eps state.logpC[1] logpriorC(state.C, state.cl, state.k,
                                             priorC) ε
@test_approx_eq_eps state.logpC[2] logcondpostC(state.C, state.cl,
                                                state.k, state.v,
                                                state.n1s, priorC) ε

correctloglik = 0.0
lidx = 0
for b in 1:m
  lidx = support.C[1, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[1, b], v[1], priorC.A[lidx], priorC.B[lidx])
  lidx = support.C[2, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[2, b], v[2], priorC.A[lidx], priorC.B[lidx])
end

@test_approx_eq_eps state.loglik correctloglik ε

@test_approx_eq_eps state.logpp (logdPriorRow(n, k, v, priorR) +
                                 logpriorC(state.C, state.cl, state.k, priorC) +
                                 correctloglik) ε

# [2; 2; 2; 1; 1; 3] => [2; 2; 3; 1; 1; 3]
support = KSupport(m, n, settings.maxclust, settings.maxunit)
state = AminoAcidState(data, [2; 2; 2; 1; 1; 3], priorR, priorC, settings)

k = 3
i = 3
hi = state.R[i]
hj = state.R[6]
cl = [1; 2; 3]
v = [2; 2; 2]
n1s = hcat(sum(float(data)[:, 4:5], 2), sum(float(data)[:, 1:2], 2),
           sum(float(data)[:, [3; 6]], 2))'

initsupportbrw!(k, i, state.v[hi], data, priorC, settings, support)

support.lograR = logratiopriorrowbrwmove(state.v[hi], state.v[hj], priorR)

support.C[1:3, :] = UInt8[4 1 2 1 1 1 3 4 2 1 2 3 1 1 1 1 2 1;
                          3 1 2 1 1 1 4 3 2 1 2 3 1 1 1 1 2 1;
                          4 1 2 1 1 1 3 3 2 1 2 4 1 1 1 1 2 1]
support.cl = [1; 2; 3]
support.k = 3

support.logpC = [logpriorC(support.C, cl, k, priorC);
                 logcondpostC(support.C, cl, k, v, n1s, priorC)]

loglikbrw!(k, hi, hj, priorC, support, state)

performbrw!(i, hi, hj, k, priorC, settings, support, state)

@test state.R == [2; 2; 3; 1; 1; 3]
@test state.k == 3

@test state.C[1, :] == support.C[1, :]
@test !state.emptycluster[1]
@test state.cl[1] == 1
@test state.v[1] == 2
@test state.n1s[1, :] == n1s[1, :]
@test state.unit[1][1:2] == [4; 5]

@test state.C[2, :] == support.C[2, :]
@test !state.emptycluster[2]
@test state.cl[2] == 2
@test state.v[2] == 2
@test state.n1s[2, :] == n1s[2, :]
@test state.unit[2][1:2] == [1; 2]

@test state.C[3, :] == support.C[3, :]
@test !state.emptycluster[3]
@test state.cl[3] == 3
@test state.v[3] == 2
@test state.n1s[3, :] == n1s[3, :]
@test state.unit[3][1:2] == [6; 3]

@test_approx_eq_eps state.logpR logdPriorRow(n, k, v, priorR) ε

@test_approx_eq_eps state.logpC[1] logpriorC(state.C, state.cl, state.k,
                                             priorC) ε
@test_approx_eq_eps state.logpC[2] logcondpostC(state.C, state.cl,
                                                state.k, state.v,
                                                state.n1s, priorC) ε

correctloglik = 0.0
lidx = 0
for b in 1:m
  lidx = support.C[1, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[1, b], v[1], priorC.A[lidx], priorC.B[lidx])
  lidx = support.C[2, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[2, b], v[2], priorC.A[lidx], priorC.B[lidx])
  lidx = support.C[3, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[3, b], v[3], priorC.A[lidx], priorC.B[lidx])
end

@test_approx_eq_eps state.loglik correctloglik ε

@test_approx_eq_eps state.logpp (logdPriorRow(n, k, v, priorR) +
                                 logpriorC(state.C, state.cl, state.k, priorC) +
                                 correctloglik) ε

# [3; 3; 3; 2; 2; 4] => [3; 3; 1; 2; 2; 4]
support = KSupport(m, n, settings.maxclust, settings.maxunit)
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

k = 4
i = 3
hi = state.R[i]
hj = 1
cl = [1; 2; 3; 4]
v = [2; 1; 2; 1]
n1s = hcat(sum(float(data)[:, 4:5], 2), float(data)[:, 6],
           sum(float(data)[:, 1:2], 2), float(data)[:, 3])'

initsupportbrw!(k, i, state.v[hi], data, priorC, settings, support)

support.lograR = logratiopriorrowbrwsplit(k, state.v[hi], priorR)

support.C[1:4, :] = UInt8[4 1 2 1 1 1 3 4 2 1 2 3 1 1 1 3 2 1;
                          3 1 2 1 1 1 4 3 2 1 2 3 1 1 1 3 2 1;
                          4 1 2 1 1 1 3 3 2 1 2 4 1 1 1 3 2 1;
                          3 1 2 1 1 1 3 3 2 1 2 3 1 1 1 4 2 1]
support.cl[1:4] = [2; 4; 3; 1]
support.k = 4

support.logpC = [logpriorC(support.C, cl, k, priorC);
                 logcondpostC(support.C, cl, k, v, n1s, priorC)]

loglikbrw!(k, hi, hj, priorC, support, state)

performbrw!(i, hi, hj, k, priorC, settings, support, state)

@test state.R == [3; 3; 1; 2; 2; 4]
@test state.k == 4

@test state.C[1, :] == support.C[4, :]
@test !state.emptycluster[1]
@test state.cl[1] == 1
@test state.v[1] == 1
@test state.n1s[1, :] == n1s[4, :]
@test state.unit[1][1] == 3

@test state.C[2, :] == support.C[1, :]
@test !state.emptycluster[2]
@test state.cl[2] == 2
@test state.v[2] == 2
@test state.n1s[2, :] == n1s[1, :]
@test state.unit[2][1:2] == [4; 5]

@test state.C[3, :] == support.C[3, :]
@test !state.emptycluster[3]
@test state.cl[3] == 3
@test state.v[3] == 2
@test state.n1s[3, :] == n1s[3, :]
@test state.unit[3][1:2] == [1; 2]

@test state.C[4, :] == support.C[2, :]
@test !state.emptycluster[4]
@test state.cl[4] == 4
@test state.v[4] == 1
@test state.n1s[4, :] == n1s[2, :]
@test state.unit[4][1] == 6

@test_approx_eq_eps state.logpR logdPriorRow(n, k, v, priorR) ε

@test_approx_eq_eps state.logpC[1] logpriorC(state.C, state.cl, state.k,
                                             priorC) ε
@test_approx_eq_eps state.logpC[2] logcondpostC(state.C, state.cl,
                                                state.k, state.v,
                                                state.n1s, priorC) ε

correctloglik = 0.0
lidx = 0
for b in 1:m
  lidx = support.C[1, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[1, b], v[1], priorC.A[lidx], priorC.B[lidx])
  lidx = support.C[2, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[2, b], v[2], priorC.A[lidx], priorC.B[lidx])
  lidx = support.C[3, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[3, b], v[3], priorC.A[lidx], priorC.B[lidx])
  lidx = support.C[4, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[4, b], v[4], priorC.A[lidx], priorC.B[lidx])
end

@test_approx_eq_eps state.loglik correctloglik ε

@test_approx_eq_eps state.logpp (logdPriorRow(n, k, v, priorR) +
                                 logpriorC(state.C, state.cl, state.k, priorC) +
                                 correctloglik) ε

# [2; 2; 2; 1; 1; 3] => [2; 2; 4; 1; 1; 3]
# allocate new resources
settings = KSettings(ifile, ofile, maxclust=3, maxunit=1, op=[0.0; 1.0; 0.0])

priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, settings.γ, settings.r,
                           maxclust=settings.maxclust)

support = KSupport(m, n, settings.maxclust, settings.maxunit)
state = AminoAcidState(data, [2; 2; 2; 1; 1; 3], priorR, priorC, settings)

k = 4
i = 3
hi = state.R[i]
hj = 0
cl = [1; 2; 3; 4]
v = [2; 1; 2; 1]
n1s = hcat(sum(float(data)[:, 4:5], 2), float(data)[:, 6],
           sum(float(data)[:, 1:2], 2), float(data)[:, 3])'

initsupportbrw!(k, i, state.v[hi], data, priorC, settings, support)

support.lograR = logratiopriorrowbrwsplit(k, state.v[hi], priorR)

support.C[1:4, :] = UInt8[4 1 2 1 1 1 3 4 2 1 2 3 1 1 1 3 2 1;
                          3 1 2 1 1 1 4 3 2 1 2 3 1 1 1 3 2 1;
                          4 1 2 1 1 1 3 3 2 1 2 4 1 1 1 3 2 1;
                          3 1 2 1 1 1 3 3 2 1 2 3 1 1 1 4 2 1]
support.cl[1:4] = [1; 3; 2; 4]
support.k = 4

support.logpC = [logpriorC(support.C, cl, k, priorC);
                 logcondpostC(support.C, cl, k, v, n1s, priorC)]

loglikbrw!(k, hi, hj, priorC, support, state)

performbrw!(i, hi, hj, k, priorC, settings, support, state)

@test state.R == [2; 2; 4; 1; 1; 3]
@test state.k == 4

@test state.C[1, :] == support.C[1, :]
@test !state.emptycluster[1]
@test state.cl[1] == 1
@test state.v[1] == 2
@test state.n1s[1, :] == n1s[1, :]
@test state.unit[1][1:2] == [4; 5]

@test state.C[2, :] == support.C[3, :]
@test !state.emptycluster[2]
@test state.cl[2] == 2
@test state.v[2] == 2
@test state.n1s[2, :] == n1s[3, :]
@test state.unit[2][1:2] == [1; 2]

@test state.C[3, :] == support.C[2, :]
@test !state.emptycluster[3]
@test state.cl[3] == 3
@test state.v[3] == 1
@test state.n1s[3, :] == n1s[2, :]
@test state.unit[3][1] == 6

@test state.C[4, :] == support.C[4, :]
@test !state.emptycluster[4]
@test state.cl[4] == 4
@test state.v[4] == 1
@test state.n1s[4, :] == n1s[4, :]
@test state.unit[4][1] == 3

@test_approx_eq_eps state.logpR logdPriorRow(n, k, v, priorR) ε

@test_approx_eq_eps state.logpC[1] logpriorC(state.C, state.cl, state.k,
                                             priorC) ε
@test_approx_eq_eps state.logpC[2] logcondpostC(state.C, state.cl,
                                                state.k, state.v,
                                                state.n1s, priorC) ε

correctloglik = 0.0
lidx = 0
for b in 1:m
  lidx = support.C[1, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[1, b], v[1], priorC.A[lidx], priorC.B[lidx])
  lidx = support.C[2, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[2, b], v[2], priorC.A[lidx], priorC.B[lidx])
  lidx = support.C[3, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[3, b], v[3], priorC.A[lidx], priorC.B[lidx])
  lidx = support.C[4, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[4, b], v[4], priorC.A[lidx], priorC.B[lidx])
end

@test_approx_eq_eps state.loglik correctloglik ε

@test_approx_eq_eps state.logpp (logdPriorRow(n, k, v, priorR) +
                                 logpriorC(state.C, state.cl, state.k, priorC) +
                                 correctloglik) ε
