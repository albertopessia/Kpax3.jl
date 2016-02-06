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
R = [1; 1; 1; 2; 2; 3]
k = length(unique(R)) - 1

priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, k + 1, settings.γ, settings.r)

mcmcobj = AminoAcidMCMC(data, R, priorR, priorC, settings)

ij = [4; 6]
S = 1

hi = mcmcobj.R[ij[1]]
hj = mcmcobj.R[ij[2]]

vi = mcmcobj.v[hi] + mcmcobj.v[hj]
ni = zeros(Float64, size(data, 1))
for b in 1:size(data, 1)
  ni[b] = mcmcobj.n1s[hi, b] + mcmcobj.n1s[hj, b]
end

# test merge functions by merging cluster 2 and cluster 3
mergesupport = KSupport(m, n, settings.maxclust, settings.maxunit)

initsupport!(ij, S, k, data, priorC, mergesupport)

# simulate a merge by first simulating its inverse split
lcp = zeros(Float64, 2)
for b in 1:m
  lcp[1] += computeclusteriseqprobs!(data[b, 5], b, priorC, mergesupport)
  lcp[2] += computeclusterjseqprobs!(data[b, 5], b, priorC, mergesupport)
end
updateclusteri!(5, data, mergesupport)

# simulate new matrix C
simcmerge!(k, hi, hj, vi, ni, priorC, mergesupport, mcmcobj)

loglikmerge!(hi, hj, ni, vi, priorC, mergesupport, mcmcobj)

correctloglik = 0.0
lidx = 0
for b in 1:m
  lidx = mergesupport.C[1, b] + 4 * (b - 1)
  correctloglik += logmarglik(mcmcobj.n1s[1, b], mcmcobj.v[1], priorC.A[lidx],
                              priorC.B[lidx])
  lidx = mergesupport.C[2, b] + 4 * (b - 1)
  correctloglik += logmarglik(ni[b], vi, priorC.A[lidx], priorC.B[lidx])
end

@test_approx_eq_eps correctloglik mergesupport.loglik ε

# split
R = [1; 1; 1; 1; 1; 2]
k = length(unique(R)) + 1

priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, k - 1, settings.γ, settings.r)

mcmcobj = AminoAcidMCMC(data, R, priorR, priorC, settings)

ij = [1; 5]
S = 3

hi = mcmcobj.R[ij[1]]

# test split functions by splitting cluster 1
splitsupport = KSupport(m, n, settings.maxclust, settings.maxunit)

initsupport!(ij, S, k, data, priorC, splitsupport)

# simulate split
# put unit 2 into cluster 1
lcp = zeros(Float64, 2)
for b in 1:m
  lcp[1] += computeclusteriseqprobs!(data[b, 2], b, priorC, splitsupport)
  lcp[2] += computeclusterjseqprobs!(data[b, 2], b, priorC, splitsupport)
end
updateclusteri!(2, data, splitsupport)

# put unit 3 into cluster 1
lcp[1] = lcp[2] = 0.0
for b in 1:m
  lcp[1] += computeclusteriseqprobs!(data[b, 3], b, priorC, splitsupport)
  lcp[2] += computeclusterjseqprobs!(data[b, 3], b, priorC, splitsupport)
end
updateclusteri!(3, data, splitsupport)

# put unit 4 into cluster 2
lcp[1] = lcp[2] = 0.0
for b in 1:m
  lcp[1] += computeclusteriseqprobs!(data[b, 4], b, priorC, splitsupport)
  lcp[2] += computeclusterjseqprobs!(data[b, 4], b, priorC, splitsupport)
end
updateclusterj!(4, data, splitsupport)

# simulate new matrix C
simcsplit!(k, 1, priorC, splitsupport, mcmcobj)

logliksplit!(hi, priorC, splitsupport, mcmcobj)

correctloglik = 0.0
lidx = 0
for b in 1:m
  lidx = splitsupport.C[1, b] + 4 * (b - 1)
  correctloglik += logmarglik(splitsupport.ni[b], splitsupport.vi,
                              priorC.A[lidx], priorC.B[lidx])
  lidx = splitsupport.C[2, b] + 4 * (b - 1)
  correctloglik += logmarglik(mcmcobj.n1s[2, b], mcmcobj.v[2], priorC.A[lidx],
                              priorC.B[lidx])
  lidx = splitsupport.C[3, b] + 4 * (b - 1)
  correctloglik += logmarglik(splitsupport.nj[b], splitsupport.vj,
                              priorC.A[lidx], priorC.B[lidx])
end

@test_approx_eq_eps correctloglik splitsupport.loglik ε
