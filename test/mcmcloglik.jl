# This file is part of Kpax3. License is MIT.

ε = eps()

settings = KSettings("outfile.bin", 1, 0, 1, [1.0; 0.0; 0.0], 0.0, 1.0,
                     [0.6; 0.35; 0.05], 135.0, 1.0, 1.0, 5.0, 6, 6, true, 1)

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

# merge
# [1; 1; 1; 2; 2; 3] => [1; 1; 1; 3; 3; 3]
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, 3, settings.γ, settings.r)
mcmcobj = AminoAcidMCMC(data, [1; 1; 1; 2; 2; 3], priorR, priorC, settings)

mergesupport = KSupport(m, n, settings.maxclust, settings.maxunit)
mergesupport.C[1:2, :] = UInt8[3 1 2 1 1 1 4 3 2 1 2 1 1 1 1 1 2 1;
                               4 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1]

R = [1; 1; 1; 3; 3; 3]
hi = 3
hj = 2
n1s = hcat(sum(float(data[:, R .== 1]), 2), sum(float(data[:, R .== 3]), 2))'
v = [3; 3]

loglikmerge!(hi, hj, vec(n1s[2, :]), v[2], priorC, mergesupport, mcmcobj)

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
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, 2, settings.γ, settings.r)
mcmcobj = AminoAcidMCMC(data, [1; 1; 1; 3; 3; 3], priorR, priorC, settings)

splitsupport = KSupport(m, n, settings.maxclust, settings.maxunit)
splitsupport.C[1:3, :] = UInt8[3 1 2 1 1 1 4 3 2 1 2 1 1 1 1 1 2 1;
                               4 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1;
                               3 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1]

R = [1; 1; 1; 3; 3; 2]
hi = 3
n1s = hcat(sum(float(data[:, R .== 1]), 2), sum(float(data[:, R .== 2]), 2),
           sum(float(data[:, R .== 3]), 2))'
v = [3; 1; 2]

splitsupport.ni = vec(n1s[3, :])
splitsupport.vi = v[3]
splitsupport.nj = vec(n1s[2, :])
splitsupport.vj = v[2]

logliksplit!(hi, priorC, splitsupport, mcmcobj)

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
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, 3, settings.γ, settings.r)
mcmcobj = AminoAcidMCMC(data, [1; 1; 1; 3; 3; 2], priorR, priorC, settings)

brwsupport = KSupport(m, n, settings.maxclust, settings.maxunit)
brwsupport.C[1:2, :] = UInt8[3 1 2 1 1 1 4 3 2 1 2 1 1 1 1 1 2 1;
                               4 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1]

R = [1; 1; 1; 3; 3; 3]
k = 2
hi = 2
hj = 3
n1s = hcat(sum(float(data[:, R .== 1]), 2), sum(float(data[:, R .== 3]), 2))'
v = [3; 3]

brwsupport.ni = vec(data[:, 6])

loglikbrw!(k, hi, hj, priorC, brwsupport, mcmcobj)

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
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, 3, settings.γ, settings.r)
mcmcobj = AminoAcidMCMC(data, [1; 1; 1; 3; 3; 2], priorR, priorC, settings)

brwsupport = KSupport(m, n, settings.maxclust, settings.maxunit)
brwsupport.C[1:3, :] = UInt8[3 1 2 1 1 1 4 3 2 1 2 1 1 1 1 1 2 1;
                               4 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1;
                               3 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1]

R = [1; 1; 2; 3; 3; 2]
k = 3
hi = 1
hj = 2
n1s = hcat(sum(float(data[:, R .== 1]), 2), sum(float(data[:, R .== 2]), 2),
           sum(float(data[:, R .== 3]), 2))'
v = [2; 2; 2]

brwsupport.ni = vec(data[:, 3])

loglikbrw!(k, hi, hj, priorC, brwsupport, mcmcobj)

correctloglik = 0.0
lidx = 0
for b in 1:m
  lidx = brwsupport.C[1, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[1, b], v[1], priorC.A[lidx], priorC.B[lidx])
  lidx = brwsupport.C[2, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[2, b], v[2], priorC.A[lidx], priorC.B[lidx])
  lidx = brwsupport.C[3, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[3, b], v[3], priorC.A[lidx], priorC.B[lidx])
end

@test_approx_eq_eps brwsupport.loglik correctloglik ε

# [1; 1; 1; 4; 4; 2] => [1; 1; 3; 4; 4; 2]
priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, 3, settings.γ, settings.r)
mcmcobj = AminoAcidMCMC(data, [1; 1; 1; 4; 4; 2], priorR, priorC, settings)

brwsupport = KSupport(m, n, settings.maxclust, settings.maxunit)
brwsupport.C[1:4, :] = UInt8[3 1 2 1 1 1 4 3 2 1 2 1 1 1 1 1 2 1;
                               4 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1;
                               3 1 2 1 1 1 3 4 2 1 2 1 1 1 1 1 2 1;
                               3 1 2 1 1 1 4 3 2 1 2 1 1 1 1 1 2 1]

R = [1; 1; 3; 4; 4; 2]
k = 4
hi = 1
hj = 3
n1s = hcat(sum(float(data[:, R .== 1]), 2), sum(float(data[:, R .== 2]), 2),
           sum(float(data[:, R .== 3]), 2), sum(float(data[:, R .== 4]), 2))'
v = [2; 1; 1; 2]

brwsupport.ni = vec(data[:, 3])

loglikbrw!(k, hi, hj, priorC, brwsupport, mcmcobj)

correctloglik = 0.0
lidx = 0
for b in 1:m
  lidx = brwsupport.C[1, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[1, b], v[1], priorC.A[lidx], priorC.B[lidx])
  lidx = brwsupport.C[2, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[2, b], v[2], priorC.A[lidx], priorC.B[lidx])
  lidx = brwsupport.C[3, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[4, b], v[4], priorC.A[lidx], priorC.B[lidx])
  lidx = brwsupport.C[4, b] + 4 * (b - 1)
  correctloglik += logmarglik(n1s[3, b], v[3], priorC.A[lidx], priorC.B[lidx])
end

@test_approx_eq_eps brwsupport.loglik correctloglik ε
