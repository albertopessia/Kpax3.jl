# This file is part of Kpax3. License is MIT.

ifile = "data/proper_aa.fasta"
ofile = "../build/test.bin"

settings = KSettings(ifile, ofile, maxclust=3, maxunit=1)

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
priorC = AminoAcidPriorCol(data, settings.γ, settings.r,
                           maxclust=settings.maxclust)

# merge
R = [13; 13; 13; 42; 42; 76]
k = length(unique(R)) - 1

ij = [4; 6]
S = 1
u = 5

# test merge functions by merging cluster 2 and cluster 3
# test initialization
mergesupport = KSupport(m, n, 1, 1)

initsupportsplitmerge!(ij, S, k, data, priorC, settings, mergesupport)

wi = zeros(Float64, 4, m)
for col in 1:m
  wi[1, col] = priorC.logγ[1] +
               logmarglik(data[col, ij[1]], 1, priorC.A[1, col],
                          priorC.B[1, col])

  wi[2, col] = priorC.logγ[2] +
               logmarglik(data[col, ij[1]], 1, priorC.A[2, col],
                          priorC.B[2, col])

  wi[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
               logmarglik(data[col, ij[1]], 1, priorC.A[3, col],
                          priorC.B[3, col])

  wi[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
               logmarglik(data[col, ij[1]], 1, priorC.A[4, col],
                          priorC.B[4, col])
end

ci = [log(sum(exp(mergesupport.wi.w[:, b])))::Float64 for b in 1:m]

wj = zeros(Float64, 4, m)
for col in 1:m
  wj[1, col] = priorC.logγ[1] +
               logmarglik(data[col, ij[2]], 1, priorC.A[1, col],
                          priorC.B[1, col])

  wj[2, col] = priorC.logγ[2] +
               logmarglik(data[col, ij[2]], 1, priorC.A[2, col],
                          priorC.B[2, col])

  wj[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
               logmarglik(data[col, ij[2]], 1, priorC.A[3, col],
                          priorC.B[3, col])

  wj[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
               logmarglik(data[col, ij[2]], 1, priorC.A[4, col],
                          priorC.B[4, col])
end

cj = [log(sum(exp(mergesupport.wj.w[:, b])))::Float64 for b in 1:m]

@test mergesupport.vi == 1
@test mergesupport.ni == float(data[:, ij[1]])
@test mergesupport.ui == [ij[1]; 0]

@test maximum(abs(mergesupport.wi.w - wi)) <= ε
@test maximum(abs(mergesupport.wi.c - ci)) <= ε
@test mergesupport.wi.z == zeros(Float64, 4, m)

@test mergesupport.vj == 1
@test mergesupport.nj == float(data[:, ij[2]])
@test mergesupport.uj == [ij[2]; 0]

@test maximum(abs(mergesupport.wj.w - wj)) <= ε
@test maximum(abs(mergesupport.wj.c - cj)) <= ε
@test mergesupport.wj.z == zeros(Float64, 4, m)

# test the move of the first unit (u) to cluster 2 (test inverse split operator)
mergesupport = KSupport(m, n, 1, 1)

initsupportsplitmerge!(ij, S, k, data, priorC, settings, mergesupport)

lcp = zeros(Float64, 2)
for b in 1:m
  lcp[1] += computeclusteriseqprobs!(data[b, u], b, priorC, mergesupport)
  lcp[2] += computeclusterjseqprobs!(data[b, u], b, priorC, mergesupport)
end
updateclusteri!(u, data, mergesupport)

wi = zeros(Float64, 4, m)
for col in 1:m
  wi[1, col] = priorC.logγ[1] +
               logmarglik(data[col, ij[1]], 1, priorC.A[1, col],
                          priorC.B[1, col]) +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[1, col], priorC.B[1, col])

  wi[2, col] = priorC.logγ[2] +
               logmarglik(data[col, ij[1]], 1, priorC.A[2, col],
                          priorC.B[2, col]) +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[2, col], priorC.B[2, col])

  wi[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
               logmarglik(data[col, ij[1]], 1, priorC.A[3, col],
                          priorC.B[3, col]) +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[3, col], priorC.B[3, col])

  wi[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
               logmarglik(data[col, ij[1]], 1, priorC.A[4, col],
                          priorC.B[4, col]) +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[4, col], priorC.B[4, col])
end

zi = copy(wi)

ci = [log(sum(exp(mergesupport.wi.w[:, b])))::Float64 for b in 1:m]

wj = zeros(Float64, 4, m)
zj = zeros(Float64, 4, m)
for col in 1:m
  wj[1, col] = priorC.logγ[1] +
               logmarglik(data[col, ij[2]], 1, priorC.A[1, col],
                          priorC.B[1, col])

  wj[2, col] = priorC.logγ[2] +
               logmarglik(data[col, ij[2]], 1, priorC.A[2, col],
                          priorC.B[2, col])

  wj[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
               logmarglik(data[col, ij[2]], 1, priorC.A[3, col],
                          priorC.B[3, col])

  wj[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
               logmarglik(data[col, ij[2]], 1, priorC.A[4, col],
                          priorC.B[4, col])

  zj[1, col] = wj[1, col] +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[1, col], priorC.B[1, col])

  zj[2, col] = wj[2, col] +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[2, col], priorC.B[2, col])

  zj[3, col] = wj[3, col] +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[3, col], priorC.B[3, col])

  zj[4, col] = wj[4, col] +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[4, col], priorC.B[4, col])

end

cj = [log(sum(exp(mergesupport.wj.w[:, b])))::Float64 for b in 1:m]

@test mergesupport.vi == 2
@test mergesupport.ni == float(data[:, ij[1]]) + float(data[:, u])
@test mergesupport.ui == [ij[1]; u]

@test maximum(abs(mergesupport.wi.w - wi)) <= ε
@test maximum(abs(mergesupport.wi.c - ci)) <= ε
@test maximum(abs(mergesupport.wi.z - zi)) <= ε

@test mergesupport.vj == 1
@test mergesupport.nj == float(data[:, ij[2]])
@test mergesupport.uj == [ij[2]; 0]

@test maximum(abs(mergesupport.wj.w - wj)) <= ε
@test maximum(abs(mergesupport.wj.c - cj)) <= ε
@test maximum(abs(mergesupport.wj.z - zj)) <= ε

# test the move of the first unit (u) to cluster 3 (test inverse split operator)
mergesupport = KSupport(m, n, 1, 1)

initsupportsplitmerge!(ij, S, k, data, priorC, settings, mergesupport)

lcp = zeros(Float64, 2)
for b in 1:m
  lcp[1] += computeclusteriseqprobs!(data[b, u], b, priorC, mergesupport)
  lcp[2] += computeclusterjseqprobs!(data[b, u], b, priorC, mergesupport)
end
updateclusterj!(u, data, mergesupport)

wi = zeros(Float64, 4, m)
zi = zeros(Float64, 4, m)
for col in 1:m
  wi[1, col] = priorC.logγ[1] +
               logmarglik(data[col, ij[1]], 1, priorC.A[1, col],
                          priorC.B[1, col])

  wi[2, col] = priorC.logγ[2] +
               logmarglik(data[col, ij[1]], 1, priorC.A[2, col],
                          priorC.B[2, col])

  wi[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
               logmarglik(data[col, ij[1]], 1, priorC.A[3, col],
                          priorC.B[3, col])

  wi[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
               logmarglik(data[col, ij[1]], 1, priorC.A[4, col],
                          priorC.B[4, col])

  zi[1, col] = wi[1, col] +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[1, col], priorC.B[1, col])

  zi[2, col] = wi[2, col] +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[2, col], priorC.B[2, col])

  zi[3, col] = wi[3, col] +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[3, col], priorC.B[3, col])

  zi[4, col] = wi[4, col] +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[4, col], priorC.B[4, col])
end

ci = [log(sum(exp(mergesupport.wi.w[:, b])))::Float64 for b in 1:m]

wj = zeros(Float64, 4, m)
for col in 1:m
  wj[1, col] = priorC.logγ[1] +
               logmarglik(data[col, ij[2]], 1, priorC.A[1, col],
                          priorC.B[1, col]) +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[1, col], priorC.B[1, col])

  wj[2, col] = priorC.logγ[2] +
               logmarglik(data[col, ij[2]], 1, priorC.A[2, col],
                          priorC.B[2, col]) +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[2, col], priorC.B[2, col])

  wj[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
               logmarglik(data[col, ij[2]], 1, priorC.A[3, col],
                          priorC.B[3, col]) +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[3, col], priorC.B[3, col])

  wj[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
               logmarglik(data[col, ij[2]], 1, priorC.A[4, col],
                          priorC.B[4, col]) +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[4, col], priorC.B[4, col])
end
zj = copy(wj)

cj = [log(sum(exp(mergesupport.wj.w[:, b])))::Float64 for b in 1:m]

@test mergesupport.vi == 1
@test mergesupport.ni == float(data[:, ij[1]])
@test mergesupport.ui == [ij[1]; 0]

@test maximum(abs(mergesupport.wi.w - wi)) <= ε
@test maximum(abs(mergesupport.wi.c - ci)) <= ε
@test maximum(abs(mergesupport.wi.z - zi)) <= ε

@test mergesupport.vj == 2
@test mergesupport.nj == float(data[:, ij[2]]) + float(data[:, u])
@test mergesupport.uj == [ij[2]; u]

@test maximum(abs(mergesupport.wj.w - wj)) <= ε
@test maximum(abs(mergesupport.wj.c - cj)) <= ε
@test maximum(abs(mergesupport.wj.z - zj)) <= ε

# split
R = [13; 13; 13; 13; 13; 76]
k = length(unique(R)) + 1

priorC = AminoAcidPriorCol(data, settings.γ, settings.r, maxclust=2)

ij = [1; 5]
S = 3
u = 4

# test split functions by splitting cluster 1
# test initialization
splitsupport = KSupport(m, n, 1, 1)

initsupportsplitmerge!(ij, S, k, data, priorC, settings, splitsupport)

len = k + settings.maxclust - 1
logω = Vector{Float64}[[log(k - 1) - log(k); -log(k)] for k in 1:len]
@test priorC.logω == logω

wi = zeros(Float64, 4, m)
for col in 1:m
  wi[1, col] = priorC.logγ[1] +
               logmarglik(data[col, ij[1]], 1, priorC.A[1, col],
                          priorC.B[1, col])

  wi[2, col] = priorC.logγ[2] +
               logmarglik(data[col, ij[1]], 1, priorC.A[2, col],
                          priorC.B[2, col])

  wi[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
               logmarglik(data[col, ij[1]], 1, priorC.A[3, col],
                          priorC.B[3, col])

  wi[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
               logmarglik(data[col, ij[1]], 1, priorC.A[4, col],
                          priorC.B[4, col])
end

ci = [log(sum(exp(splitsupport.wi.w[:, b])))::Float64 for b in 1:m]

wj = zeros(Float64, 4, m)
for col in 1:m
  wj[1, col] = priorC.logγ[1] +
               logmarglik(data[col, ij[2]], 1, priorC.A[1, col],
                          priorC.B[1, col])

  wj[2, col] = priorC.logγ[2] +
               logmarglik(data[col, ij[2]], 1, priorC.A[2, col],
                          priorC.B[2, col])

  wj[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
               logmarglik(data[col, ij[2]], 1, priorC.A[3, col],
                          priorC.B[3, col])

  wj[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
               logmarglik(data[col, ij[2]], 1, priorC.A[4, col],
                          priorC.B[4, col])
end

cj = [log(sum(exp(splitsupport.wj.w[:, b])))::Float64 for b in 1:m]

@test splitsupport.vi == 1
@test splitsupport.ni == float(data[:, ij[1]])
@test splitsupport.ui == [ij[1]; zeros(Int, S)]

@test maximum(abs(splitsupport.wi.w - wi)) <= ε
@test maximum(abs(splitsupport.wi.c - ci)) <= ε
@test splitsupport.wi.z == zeros(Float64, 4, m)

@test splitsupport.vj == 1
@test splitsupport.nj == float(data[:, ij[2]])
@test splitsupport.uj == [ij[2]; zeros(Int, S)]

@test maximum(abs(splitsupport.wj.w - wj)) <= ε
@test maximum(abs(splitsupport.wj.c - cj)) <= ε
@test splitsupport.wj.z == zeros(Float64, 4, m)

# test the move of the first unit (u) to cluster 1
priorC = AminoAcidPriorCol(data, settings.γ, settings.r, maxclust=2)

splitsupport = KSupport(m, n, 1, 1)

initsupportsplitmerge!(ij, S, k, data, priorC, settings, splitsupport)

len = k + settings.maxclust - 1
logω = Vector{Float64}[[log(k - 1) - log(k); -log(k)] for k in 1:len]
@test priorC.logω == logω

lcp = zeros(Float64, 2)
for b in 1:m
  lcp[1] += computeclusteriseqprobs!(data[b, u], b, priorC, splitsupport)
  lcp[2] += computeclusterjseqprobs!(data[b, u], b, priorC, splitsupport)
end
updateclusteri!(u, data, splitsupport)

wi = zeros(Float64, 4, m)
for col in 1:m
  wi[1, col] = priorC.logγ[1] +
               logmarglik(data[col, ij[1]], 1, priorC.A[1, col],
                          priorC.B[1, col]) +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[1, col], priorC.B[1, col])

  wi[2, col] = priorC.logγ[2] +
               logmarglik(data[col, ij[1]], 1, priorC.A[2, col],
                          priorC.B[2, col]) +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[2, col], priorC.B[2, col])

  wi[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
               logmarglik(data[col, ij[1]], 1, priorC.A[3, col],
                          priorC.B[3, col]) +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[3, col], priorC.B[3, col])

  wi[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
               logmarglik(data[col, ij[1]], 1, priorC.A[4, col],
                          priorC.B[4, col]) +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[4, col], priorC.B[4, col])
end
zi = copy(wi)

ci = [log(sum(exp(splitsupport.wi.w[:, b])))::Float64 for b in 1:m]

wj = zeros(Float64, 4, m)
zj = zeros(Float64, 4, m)
for col in 1:m
  wj[1, col] = priorC.logγ[1] +
               logmarglik(data[col, ij[2]], 1, priorC.A[1, col],
                          priorC.B[1, col])

  wj[2, col] = priorC.logγ[2] +
               logmarglik(data[col, ij[2]], 1, priorC.A[2, col],
                          priorC.B[2, col])

  wj[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
               logmarglik(data[col, ij[2]], 1, priorC.A[3, col],
                          priorC.B[3, col])

  wj[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
               logmarglik(data[col, ij[2]], 1, priorC.A[4, col],
                          priorC.B[4, col])

  zj[1, col] = wj[1, col] +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[1, col], priorC.B[1, col])

  zj[2, col] = wj[2, col] +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[2, col], priorC.B[2, col])

  zj[3, col] = wj[3, col] +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[3, col], priorC.B[3, col])

  zj[4, col] = wj[4, col] +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[4, col], priorC.B[4, col])
end

cj = [log(sum(exp(splitsupport.wj.w[:, b])))::Float64 for b in 1:m]

@test splitsupport.vi == 2
@test splitsupport.ni == float(data[:, ij[1]]) + float(data[:, u])
@test splitsupport.ui == [ij[1]; u; 0; 0]

@test maximum(abs(splitsupport.wi.w - wi)) <= ε
@test maximum(abs(splitsupport.wi.c - ci)) <= ε
@test maximum(abs(splitsupport.wi.z - zi)) <= ε

@test splitsupport.vj == 1
@test splitsupport.nj == float(data[:, ij[2]])
@test splitsupport.uj == [ij[2]; 0; 0; 0]

@test maximum(abs(splitsupport.wj.w - wj)) <= ε
@test maximum(abs(splitsupport.wj.c - cj)) <= ε
@test maximum(abs(splitsupport.wj.z - zj)) <= ε

# test the move of the first unit (u) to cluster 2
priorC = AminoAcidPriorCol(data, settings.γ, settings.r, maxclust=2)

splitsupport = KSupport(m, n, 1, 1)

initsupportsplitmerge!(ij, S, k, data, priorC, settings, splitsupport)

len = k + settings.maxclust - 1
logω = Vector{Float64}[[log(k - 1) - log(k); -log(k)] for k in 1:len]
@test priorC.logω == logω

lcp = zeros(Float64, 2)
for b in 1:m
  lcp[1] += computeclusteriseqprobs!(data[b, u], b, priorC, splitsupport)
  lcp[2] += computeclusterjseqprobs!(data[b, u], b, priorC, splitsupport)
end
updateclusterj!(u, data, splitsupport)

wi = zeros(Float64, 4, m)
zi = zeros(Float64, 4, m)
for col in 1:m
  wi[1, col] = priorC.logγ[1] +
               logmarglik(data[col, ij[1]], 1, priorC.A[1, col],
                          priorC.B[1, col])

  wi[2, col] = priorC.logγ[2] +
               logmarglik(data[col, ij[1]], 1, priorC.A[2, col],
                          priorC.B[2, col])

  wi[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
               logmarglik(data[col, ij[1]], 1, priorC.A[3, col],
                          priorC.B[3, col])

  wi[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
               logmarglik(data[col, ij[1]], 1, priorC.A[4, col],
                          priorC.B[4, col])

  zi[1, col] = wi[1, col] +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[1, col], priorC.B[1, col])

  zi[2, col] = wi[2, col] +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[2, col], priorC.B[2, col])

  zi[3, col] = wi[3, col] +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[3, col], priorC.B[3, col])

  zi[4, col] = wi[4, col] +
               logcondmarglik(data[col, u], data[col, ij[1]], 1,
                              priorC.A[4, col], priorC.B[4, col])
end

ci = [log(sum(exp(splitsupport.wi.w[:, b])))::Float64 for b in 1:m]

wj = zeros(Float64, 4, m)
for col in 1:m
  wj[1, col] = priorC.logγ[1] +
               logmarglik(data[col, ij[2]], 1, priorC.A[1, col],
                          priorC.B[1, col]) +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[1, col], priorC.B[1, col])

  wj[2, col] = priorC.logγ[2] +
               logmarglik(data[col, ij[2]], 1, priorC.A[2, col],
                          priorC.B[2, col]) +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[2, col], priorC.B[2, col])

  wj[3, col] = priorC.logγ[3] + priorC.logω[k][1] +
               logmarglik(data[col, ij[2]], 1, priorC.A[3, col],
                          priorC.B[3, col]) +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[3, col], priorC.B[3, col])

  wj[4, col] = priorC.logγ[3] + priorC.logω[k][2] +
               logmarglik(data[col, ij[2]], 1, priorC.A[4, col],
                          priorC.B[4, col]) +
               logcondmarglik(data[col, u], data[col, ij[2]], 1,
                              priorC.A[4, col], priorC.B[4, col])
end
zj = copy(wj)

cj = [log(sum(exp(splitsupport.wj.w[:, b])))::Float64 for b in 1:m]

@test splitsupport.vi == 1
@test splitsupport.ni == float(data[:, ij[1]])
@test splitsupport.ui == [ij[1]; 0; 0; 0]

@test maximum(abs(splitsupport.wi.w - wi)) <= ε
@test maximum(abs(splitsupport.wi.c - ci)) <= ε
@test maximum(abs(splitsupport.wi.z - zi)) <= ε

@test splitsupport.vj == 2
@test splitsupport.nj == float(data[:, ij[2]]) + float(data[:, u])
@test splitsupport.uj == [ij[2]; u; 0; 0]

@test maximum(abs(splitsupport.wj.w - wj)) <= ε
@test maximum(abs(splitsupport.wj.c - cj)) <= ε
@test maximum(abs(splitsupport.wj.z - zj)) <= ε

# biased random walk
# move unit 4 into cluster 1
priorC = AminoAcidPriorCol(data, settings.γ, settings.r, maxclust=3)
brwsupport = KSupport(m, n, 1, 1)

# new_k = 2; i = 4; vi = 1
initsupportbrw!(2, 4, 1, data, priorC, settings, brwsupport)

@test brwsupport.ni == float(data[:, 4])
@test brwsupport.vi == 1

# move unit 3 into its own cluster
priorC = AminoAcidPriorCol(data, settings.γ, settings.r, maxclust=3)
brwsupport = KSupport(m, n, 1, 1)

# new_k = 3; i = 3; vi = 3
initsupportbrw!(3, 3, 3, data, priorC, settings, brwsupport)

@test brwsupport.ni == float(data[:, 3])
@test brwsupport.vi == 3

# move unit 3 into its own cluster
priorC = AminoAcidPriorCol(data, settings.γ, settings.r, maxclust=2)
brwsupport = KSupport(m, n, 1, 1)

# new_k = 3; i = 3; vi = 3
initsupportbrw!(3, 3, 3, data, priorC, settings, brwsupport)

len = 3 + settings.maxclust - 1
logω = Vector{Float64}[[log(k - 1) - log(k); -log(k)] for k in 1:len]
@test priorC.logω == logω

@test brwsupport.ni == float(data[:, 3])
@test brwsupport.vi == 3
