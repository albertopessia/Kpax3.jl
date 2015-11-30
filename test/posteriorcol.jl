# This file is part of Kpax3. License is MIT.

ε = eps()

data = UInt8[0x01 0x01 0x01 0x01 0x00 0x00;
             0x00 0x00 0x00 0x00 0x01 0x01;
             0x01 0x01 0x00 0x00 0x01 0x00;
             0x01 0x01 0x00 0x00 0x00 0x00]

m, n = size(data)

maxclust = 100

α = 0.0
θ = 1.0
priorR = EwensPitman(α, θ)

R = [1, 1, 1, 1, 2, 2]

k = maximum(R)

g = [0.6, 0.35, 0.05]
r = log(0.001) / log(0.95)
priorC = AminoAcidPriorCol(data, k, g, r)

C = zeros(UInt8, m, maxclust)

emptycluster = trues(maxclust)
emptycluster[13] = false
emptycluster[76] = false

cluster = [Cluster(0, zeros(Int, 1), zeros(Float64, 1), 0.0)
           for g in 1:maxclust]

cluster[13].v = sum(R .== 1)
cluster[13].unit = find(R .== 1)
cluster[13].n1s = zeros(Float64, m)
for i in cluster[13].unit
  cluster[13].n1s += data[:, i]
end

cluster[76].v = sum(R .== 2)
cluster[76].unit = find(R .== 2)
cluster[76].n1s = zeros(Float64, m)
for i in cluster[76].unit
  cluster[76].n1s += data[:, i]
end

p = zero(Float64)

cs = ([0x01 0x01], [0x02 0x02], [0x03 0x03],
      [0x03 0x04], [0x04 0x03], [0x04 0x04])

for c1 in cs, c2 in cs, c3 in cs, c4 in cs
  tmp = [c1; c2; c3; c4]

  C[:, 13] = tmp[:, 1]
  C[:, 76] = tmp[:, 2]

  p += exp(logcondpostC(C, priorC, cluster, emptycluster))
end

@test_approx_eq_eps p 1.0 ε

p = zero(Float64)

ss = [0x01, 0x02, 0x03]

for s1 in ss, s2 in ss, s3 in ss, s4 in ss
  S = [s1, s2, s3, s4]
  p += exp(logcondpostS(S, priorC, cluster, emptycluster))
end

@test_approx_eq_eps p 1.0 ε

N = 1000000
for i in 1:N
  rcolpartition!(C, priorC, cluster, emptycluster)
end

Ctest = [0x04 0x03; 0x03 0x04; 0x01 0x01; 0x02 0x02]
