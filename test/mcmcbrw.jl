# This file is part of Kpax3. License is MIT.

ε = eps()

settings = KSettings("typesmcmc.bin", 1, 0, 1, [0.0; 1.0; 0.0], 0.0, 1.0,
                     [0.6; 0.35; 0.05], 135.0, 1.0, 1.0, 5.0, 15, 1, true, 1)

data = UInt8[0x00 0x00 0x00 0x00 0x00 0x01;
             0x01 0x01 0x01 0x01 0x01 0x00;
             0x00 0x00 0x01 0x00 0x01 0x01;
             0x01 0x01 0x00 0x01 0x00 0x00;
             0x01 0x01 0x00 0x00 0x00 0x00;
             0x00 0x00 0x00 0x01 0x01 0x00;
             0x01 0x01 0x01 0x00 0x00 0x00;
             0x00 0x00 0x00 0x01 0x01 0x01;
             0x00 0x00 0x01 0x00 0x00 0x00;
             0x01 0x00 0x00 0x01 0x00 0x01;
             0x00 0x01 0x00 0x00 0x01 0x00;
             0x00 0x00 0x00 0x00 0x00 0x01;
             0x01 0x01 0x01 0x00 0x00 0x00;
             0x00 0x00 0x00 0x01 0x01 0x00;
             0x01 0x01 0x00 0x00 0x01 0x01;
             0x00 0x00 0x01 0x01 0x00 0x00;
             0x01 0x01 0x00 0x01 0x00 0x00;
             0x00 0x00 0x01 0x00 0x01 0x01]

m, n = size(data)

# [7; 7; 7; 5; 5; 9] => [7; 7; 7; 5; 5; 5]
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, 3, settings.γ, settings.r)
support = KSupport(m, n, settings.maxclust, settings.maxunit)
mcmcobj = AminoAcidMCMC(data, [7; 7; 7; 5; 5; 9], priorR, priorC, settings)

k = 2
i = 6
hi = mcmcobj.R[i]
hj = mcmcobj.R[5]
cl = [1; 2]
v = [3; 3]
n1s = hcat(sum(float(data)[:, 4:6], 2), sum(float(data)[:, 1:3], 2))'

initsupportbrw!(k, i, mcmcobj.v[hi], data, support)

support.lograR = logratiopriorrowbrwmerge(k, mcmcobj.v[hj], priorR)

support.C[1:2, :] = UInt8[4 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1;
                          3 1 2 1 1 1 4 3 2 1 2 1 1 1 1 1 2 1]

support.logpC = [logpriorC(support.C, cl, priorC.logγ, support.logω);
                 logcondpostC(support.C, cl, v, n1s, support.logω, priorC)]

loglikbrw!(k, hi, hj, priorC, support, mcmcobj)

performbrw!(i, hi, hj, k, priorC, settings, support, mcmcobj)

@test mcmcobj.R == [7; 7; 7; 5; 5; 5]

@test mcmcobj.C[5, :] == support.C[1, :]
@test mcmcobj.filledcluster[5]
@test mcmcobj.cl[1] == 5
@test mcmcobj.v[5] == 3
@test mcmcobj.n1s[5, :] == n1s[1, :]
@test mcmcobj.unit[5] == [4; 5; 6]

@test mcmcobj.C[7, :] == support.C[2, :]
@test mcmcobj.filledcluster[7]
@test mcmcobj.cl[2] == 7
@test mcmcobj.v[7] == 3
@test mcmcobj.n1s[7, :] == n1s[2, :]
@test mcmcobj.unit[7] == [1; 2; 3]

@test !mcmcobj.filledcluster[9]

@test_approx_eq_eps priorC.logω [0.0; 0.0; log(k - 1.0) - log(k); -log(k)] ε

@test_approx_eq_eps mcmcobj.logpR logdPriorRow(n, k, v, priorR) ε

@test_approx_eq_eps mcmcobj.logpC[1] logpriorC(mcmcobj.C, mcmcobj.cl,
                                               priorC.logγ, priorC.logω) ε
@test_approx_eq_eps mcmcobj.logpC[2] logcondpostC(mcmcobj.C, mcmcobj.cl,
                                                  mcmcobj.v, mcmcobj.n1s,
                                                  priorC.logω, priorC) ε

correctloglik = 0.0
lidx = 0
for b in 1:m
  lidx = support.C[1, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[1, b], v[1], priorC.A[lidx], priorC.B[lidx])
  lidx = support.C[2, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[2, b], v[2], priorC.A[lidx], priorC.B[lidx])
end

@test_approx_eq_eps mcmcobj.loglik correctloglik ε

# [7; 7; 7; 5; 5; 9] => [7; 7; 9; 5; 5; 9]
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, 3, settings.γ, settings.r)
support = KSupport(m, n, settings.maxclust, settings.maxunit)
mcmcobj = AminoAcidMCMC(data, [7; 7; 7; 5; 5; 9], priorR, priorC, settings)

k = 3
i = 3
hi = mcmcobj.R[i]
hj = mcmcobj.R[6]
cl = [1; 2; 3]
v = [2; 2; 2]
n1s = hcat(sum(float(data)[:, 4:5], 2), sum(float(data)[:, 1:2], 2),
           sum(float(data)[:, [3; 6]], 2))'

initsupportbrw!(k, i, mcmcobj.v[hi], data, support)

support.lograR = logratiopriorrowbrwmove(mcmcobj.v[hi], mcmcobj.v[hj], priorR)

support.C[1:3, :] = UInt8[4 1 2 1 1 1 3 4 2 1 2 3 1 1 1 1 2 1;
                          3 1 2 1 1 1 4 3 2 1 2 3 1 1 1 1 2 1;
                          4 1 2 1 1 1 3 3 2 1 2 4 1 1 1 1 2 1]

support.logpC = [logpriorC(support.C, cl, priorC.logγ, support.logω);
                 logcondpostC(support.C, cl, v, n1s, support.logω, priorC)]

loglikbrw!(k, hi, hj, priorC, support, mcmcobj)

performbrw!(i, hi, hj, k, priorC, settings, support, mcmcobj)

@test mcmcobj.R == [7; 7; 9; 5; 5; 9]

@test mcmcobj.C[5, :] == support.C[1, :]
@test mcmcobj.filledcluster[5]
@test mcmcobj.cl[1] == 5
@test mcmcobj.v[5] == 2
@test mcmcobj.n1s[5, :] == n1s[1, :]
@test mcmcobj.unit[5] == [4; 5]

@test mcmcobj.C[7, :] == support.C[2, :]
@test mcmcobj.filledcluster[7]
@test mcmcobj.cl[2] == 7
@test mcmcobj.v[7] == 2
@test mcmcobj.n1s[7, :] == n1s[2, :]
@test mcmcobj.unit[7] == [1; 2]

@test mcmcobj.C[9, :] == support.C[3, :]
@test mcmcobj.filledcluster[9]
@test mcmcobj.cl[3] == 9
@test mcmcobj.v[9] == 2
@test mcmcobj.n1s[9, :] == n1s[3, :]
@test mcmcobj.unit[9] == [6; 3]

@test_approx_eq_eps priorC.logω [0.0; 0.0; log(k - 1.0) - log(k); -log(k)] ε

@test_approx_eq_eps mcmcobj.logpR logdPriorRow(n, k, v, priorR) ε

@test_approx_eq_eps mcmcobj.logpC[1] logpriorC(mcmcobj.C, mcmcobj.cl,
                                               priorC.logγ, priorC.logω) ε
@test_approx_eq_eps mcmcobj.logpC[2] logcondpostC(mcmcobj.C, mcmcobj.cl,
                                                  mcmcobj.v, mcmcobj.n1s,
                                                  priorC.logω, priorC) ε

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

@test_approx_eq_eps mcmcobj.loglik correctloglik ε

# [7; 7; 7; 5; 5; 9] => [7; 7; 2; 5; 5; 9]
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, 3, settings.γ, settings.r)
support = KSupport(m, n, settings.maxclust, settings.maxunit)
mcmcobj = AminoAcidMCMC(data, [7; 7; 7; 5; 5; 9], priorR, priorC, settings)

k = 4
i = 3
hi = mcmcobj.R[i]
hj = 2
cl = [1; 2; 3; 4]
v = [2; 2; 1; 1]
n1s = hcat(sum(float(data)[:, 4:5], 2), sum(float(data)[:, 1:2], 2),
           float(data)[:, 6], float(data)[:, 3])'

initsupportbrw!(k, i, mcmcobj.v[hi], data, support)

support.lograR = logratiopriorrowbrwsplit(k, mcmcobj.v[hi], priorR)

support.C[1:4, :] = UInt8[4 1 2 1 1 1 3 4 2 1 2 3 1 1 1 3 2 1;
                          3 1 2 1 1 1 4 3 2 1 2 3 1 1 1 3 2 1;
                          4 1 2 1 1 1 3 3 2 1 2 4 1 1 1 3 2 1;
                          3 1 2 1 1 1 3 3 2 1 2 3 1 1 1 4 2 1]

support.logpC = [logpriorC(support.C, cl, priorC.logγ, support.logω);
                 logcondpostC(support.C, cl, v, n1s, support.logω, priorC)]

loglikbrw!(k, hi, hj, priorC, support, mcmcobj)

performbrw!(i, hi, hj, k, priorC, settings, support, mcmcobj)

@test mcmcobj.R == [7; 7; 2; 5; 5; 9]

@test mcmcobj.C[2, :] == support.C[4, :]
@test mcmcobj.filledcluster[2]
@test mcmcobj.cl[1] == 2
@test mcmcobj.v[2] == 1
@test mcmcobj.n1s[2, :] == n1s[4, :]
@test mcmcobj.unit[2] == [3]

@test mcmcobj.C[5, :] == support.C[1, :]
@test mcmcobj.filledcluster[5]
@test mcmcobj.cl[2] == 5
@test mcmcobj.v[5] == 2
@test mcmcobj.n1s[5, :] == n1s[1, :]
@test mcmcobj.unit[5] == [4; 5]

@test mcmcobj.C[7, :] == support.C[2, :]
@test mcmcobj.filledcluster[7]
@test mcmcobj.cl[3] == 7
@test mcmcobj.v[7] == 2
@test mcmcobj.n1s[7, :] == n1s[2, :]
@test mcmcobj.unit[7] == [1; 2]

@test mcmcobj.C[9, :] == support.C[3, :]
@test mcmcobj.filledcluster[9]
@test mcmcobj.cl[4] == 9
@test mcmcobj.v[9] == 1
@test mcmcobj.n1s[9, :] == n1s[3, :]
@test mcmcobj.unit[9] == [6]

@test_approx_eq_eps priorC.logω [0.0; 0.0; log(k - 1.0) - log(k); -log(k)] ε

@test_approx_eq_eps mcmcobj.logpR logdPriorRow(n, k, v, priorR) ε

@test_approx_eq_eps mcmcobj.logpC[1] logpriorC(mcmcobj.C, mcmcobj.cl,
                                               priorC.logγ, priorC.logω) ε
@test_approx_eq_eps mcmcobj.logpC[2] logcondpostC(mcmcobj.C, mcmcobj.cl,
                                                  mcmcobj.v, mcmcobj.n1s,
                                                  priorC.logω, priorC) ε

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

@test_approx_eq_eps mcmcobj.loglik correctloglik ε

# [2; 2; 2; 1; 1; 3] => [2; 2; 4; 1; 1; 3]
# allocate new resources
settings = KSettings("typesmcmc.bin", 1, 0, 1, [0.0; 1.0; 0.0], 0.0, 1.0,
                     [0.6; 0.35; 0.05], 135.0, 1.0, 1.0, 5.0, 3, 1, true, 1)
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, 3, settings.γ, settings.r)
support = KSupport(m, n, settings.maxclust, settings.maxunit)
mcmcobj = AminoAcidMCMC(data, [2; 2; 2; 1; 1; 3], priorR, priorC, settings)

k = 4
i = 3
hi = mcmcobj.R[i]
hj = 0
cl = [1; 2; 3; 4]
v = [2; 2; 1; 1]
n1s = hcat(sum(float(data)[:, 4:5], 2), sum(float(data)[:, 1:2], 2),
           float(data)[:, 6], float(data)[:, 3])'

initsupportbrw!(k, i, mcmcobj.v[hi], data, support)

support.lograR = logratiopriorrowbrwsplit(k, mcmcobj.v[hi], priorR)

support.C[1:4, :] = UInt8[4 1 2 1 1 1 3 4 2 1 2 3 1 1 1 3 2 1;
                          3 1 2 1 1 1 4 3 2 1 2 3 1 1 1 3 2 1;
                          4 1 2 1 1 1 3 3 2 1 2 4 1 1 1 3 2 1;
                          3 1 2 1 1 1 3 3 2 1 2 3 1 1 1 4 2 1]

support.logpC = [logpriorC(support.C, cl, priorC.logγ, support.logω);
                 logcondpostC(support.C, cl, v, n1s, support.logω, priorC)]

loglikbrw!(k, hi, hj, priorC, support, mcmcobj)

performbrw!(i, hi, hj, k, priorC, settings, support, mcmcobj)

@test mcmcobj.R == [2; 2; 4; 1; 1; 3]

@test mcmcobj.C[1, :] == support.C[1, :]
@test mcmcobj.filledcluster[1]
@test mcmcobj.cl[1] == 1
@test mcmcobj.v[1] == 2
@test mcmcobj.n1s[1, :] == n1s[1, :]
@test mcmcobj.unit[1] == [4; 5]

@test mcmcobj.C[2, :] == support.C[2, :]
@test mcmcobj.filledcluster[2]
@test mcmcobj.cl[2] == 2
@test mcmcobj.v[2] == 2
@test mcmcobj.n1s[2, :] == n1s[2, :]
@test mcmcobj.unit[2] == [1; 2]

@test mcmcobj.C[3, :] == support.C[3, :]
@test mcmcobj.filledcluster[3]
@test mcmcobj.cl[3] == 3
@test mcmcobj.v[3] == 1
@test mcmcobj.n1s[3, :] == n1s[3, :]
@test mcmcobj.unit[3] == [6]

@test mcmcobj.C[4, :] == support.C[4, :]
@test mcmcobj.filledcluster[4]
@test mcmcobj.cl[4] == 4
@test mcmcobj.v[4] == 1
@test mcmcobj.n1s[4, :] == n1s[4, :]
@test mcmcobj.unit[4] == [3]

@test_approx_eq_eps priorC.logω [0.0; 0.0; log(k - 1.0) - log(k); -log(k)] ε

@test_approx_eq_eps mcmcobj.logpR logdPriorRow(n, k, v, priorR) ε

@test_approx_eq_eps mcmcobj.logpC[1] logpriorC(mcmcobj.C, mcmcobj.cl,
                                               priorC.logγ, priorC.logω) ε
@test_approx_eq_eps mcmcobj.logpC[2] logcondpostC(mcmcobj.C, mcmcobj.cl,
                                                  mcmcobj.v, mcmcobj.n1s,
                                                  priorC.logω, priorC) ε

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

@test_approx_eq_eps mcmcobj.loglik correctloglik ε
