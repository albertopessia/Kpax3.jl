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

cl = find(!emptycluster)

cluster = [Cluster(0, zeros(Int, 1), zeros(Float64, 1), 0.0)
           for g in 1:maxclust]

cluster[cl[1]].v = sum(R .== 1)
cluster[cl[1]].unit = find(R .== 1)
cluster[cl[1]].n1s = zeros(Float64, m)
for i in cluster[cl[1]].unit
  cluster[cl[1]].n1s += data[:, i]
end

cluster[cl[2]].v = sum(R .== 2)
cluster[cl[2]].unit = find(R .== 2)
cluster[cl[2]].n1s = zeros(Float64, m)
for i in cluster[cl[2]].unit
  cluster[cl[2]].n1s += data[:, i]
end

p = zero(Float64)

cs = ([0x01 0x01], [0x02 0x02], [0x03 0x03],
      [0x03 0x04], [0x04 0x03], [0x04 0x04])

for c1 in cs, c2 in cs, c3 in cs, c4 in cs
  tmp = [c1; c2; c3; c4]

  C[:, cl[1]] = tmp[:, 1]
  C[:, cl[2]] = tmp[:, 2]

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

# configuration with the highest posterior probability
Ctest = [0x02 0x02; 0x02 0x02; 0x01 0x01; 0x01 0x01]
N = 10000000
Sp = zeros(Float64, m, 3)
Cp = zero(Float64)

trueSp = hcat([0.388699761141992228274943954602349549531936645507812500,
               0.404110835846893112766053945961175486445426940917968750,
               0.207189403011114547936699636920820921659469604492187500],
              [0.388699761141992228274943954602349549531936645507812500
               0.404110835846893112766053945961175486445426940917968750
               0.207189403011114547936699636920820921659469604492187500],
              [0.768636427030134350424361855402821674942970275878906250
               0.231359112662804122795279226920683868229389190673828125
               0.000004460307061158063528846090539659030582697596400977],
              [0.745422155055375235122028243495151400566101074218750000
               0.254411468571633325730374508566455915570259094238281250
               0.000166376372991255944535152200280947454302804544568062])'

for i in 1:N
  rcolpartition!(C, priorC, cluster, emptycluster)

  for j in 1:m
    if C[j, cl[1]] == 0x01
      @test C[j, cl[2]] == 0x01

      Sp[j, 1] += 1.0
    elseif C[j, cl[1]] == 0x02
      @test C[j, cl[2]] == 0x02

      Sp[j, 2] += 1.0
    else
      Sp[j, 3] += 1.0
    end
  end

  if C[:, cl] == Ctest
    Cp += 1.0
  end
end

Sp /= N
Cp /= N

C[:, cl] = Ctest
p = exp(logcondpostC(C, priorC, cluster, emptycluster))

@test maximum(abs(Sp - trueSp)) < 0.001
@test_approx_eq_eps Cp p 0.001
