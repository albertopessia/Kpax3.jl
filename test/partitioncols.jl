# This file is part of Kpax3. License is MIT.

ε = 1.0e-15

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
n1s = Float64[1; 5; 3; 3; 2; 2; 3; 3; 1; 3; 2; 1; 3; 2; 4; 2; 3; 3]

r1 = 2.0
r2 = 100.0

A1 = zeros(Float64, 4, m)
A1[1, :] = (r1 + 1.0) * (n1s + 0.5) / (n + 1)
A1[2, :] = 1.0
A1[3, :] = 1.0
A1[4, :] = r1

B1 = zeros(Float64, 4, m)
B1[1, :] = (r1 + 1.0) - A1[1, :]
B1[2, :] = 1.0
B1[3, :] = r1
B1[4, :] = 1.0

A2 = zeros(Float64, 4, m)
A2[1, :] = n1s + 0.5
A2[2, :] = 1.0
A2[3, :] = 1.0
A2[4, :] = r2

B2 = zeros(Float64, 4, m)
B2[1, :] = n - n1s + 0.5
B2[2, :] = 1.0
B2[3, :] = r2
B2[4, :] = 1.0

for k in 1:n
  ω = [1.0; 1.0; 1.0 - 1.0 / k; 1.0 / k]

  for γ in ([1.0; 0.0; 0.0], [0.0; 1.0; 0.0], [0.0; 0.0; 1.0],
            [0.4; 0.3; 0.3], [0.5; 0.3; 0.2], [0.7; 0.2; 0.1],
            [0.1; 0.1; 0.1], [0.3; 0.1; 0.1], [0.0; 0.2; 0.1])
    x1 = AminoAcidPriorCol(data, k, γ, r1)
    x2 = AminoAcidPriorCol(data, k, γ, r2)

    for s in 1:3
      @test_approx_eq_eps x1.γ[s] (γ[s] / sum(γ)) ε
      @test_approx_eq_eps x1.logγ[s] log(γ[s] / sum(γ)) ε
      @test_approx_eq_eps x1.ω[s] ω[s] ε
      @test_approx_eq_eps x1.logω[s] log(ω[s]) ε

      @test_approx_eq_eps x2.γ[s] (γ[s] / sum(γ)) ε
      @test_approx_eq_eps x2.logγ[s] log(γ[s] / sum(γ)) ε
      @test_approx_eq_eps x2.ω[s] ω[s] ε
      @test_approx_eq_eps x2.logω[s] log(ω[s]) ε
    end

    for b in 1:m, s in 1:4
      @test_approx_eq_eps x1.A[s, b] A1[s, b] ε
      @test_approx_eq_eps x1.B[s, b] B1[s, b] ε

      @test_approx_eq_eps x2.A[s, b] A2[s, b] ε
      @test_approx_eq_eps x2.B[s, b] B2[s, b] ε
    end
  end
end

data = UInt8[0x01 0x01 0x01 0x01 0x00 0x00;
             0x00 0x00 0x00 0x00 0x01 0x01;
             0x01 0x01 0x00 0x00 0x01 0x00;
             0x01 0x01 0x00 0x00 0x00 0x00]

m, n = size(data)

maxclust = 100

α = 0.0
θ = 1.0
priorR = EwensPitman(α, θ)

R = [1; 1; 1; 1; 2; 2]

k = maximum(R)

g = [0.6; 0.35; 0.05]
r = log(0.001) / log(0.95)
priorC = AminoAcidPriorCol(data, k, g, r)

C = zeros(UInt8, maxclust, m)

emptycluster = trues(maxclust)
emptycluster[13] = false
emptycluster[76] = false

cl = find(!emptycluster)

cluster = [KCluster(0, zeros(Int, 1), zeros(Float64, 1))
           for g in 1:maxclust]

cluster[cl[1]].v = sum(R .== 1)
cluster[cl[1]].unit = find(R .== 1)
cluster[cl[1]].n1s = zeros(Float64, m)
for a in cluster[cl[1]].unit
  cluster[cl[1]].n1s += data[:, a]
end

cluster[cl[2]].v = sum(R .== 2)
cluster[cl[2]].unit = find(R .== 2)
cluster[cl[2]].n1s = zeros(Float64, m)
for a in cluster[cl[2]].unit
  cluster[cl[2]].n1s += data[:, a]
end

cs = ([0x01; 0x01], [0x02; 0x02], [0x03; 0x03],
      [0x03; 0x04], [0x04; 0x03], [0x04; 0x04])

logp1 = zeros(Float64, (2 + 2^k)^m)
logp2 = zeros(Float64, (2 + 2^k)^m)

l = 0

for c1 in cs, c2 in cs, c3 in cs, c4 in cs
  l += 1

  tmp = hcat(c1, c2, c3, c4)

  C[cl[1], :] = tmp[1, :]
  C[cl[2], :] = tmp[2, :]

  logp1[l] = logpriorC(C, emptycluster, priorC.logγ, priorC.logω)
  logp2[l] = logcondpostC(C, cluster, emptycluster, priorC.logγ, priorC.logω,
                          priorC.A, priorC.B)
end

M = maximum(logp1)
p1 = exp(M + log(sum(exp(logp1 - M))))

M = maximum(logp2)
p2 = exp(M + log(sum(exp(logp2 - M))))

@test_approx_eq_eps p1 1.0 ε
@test_approx_eq_eps p2 1.0 ε

ss = [0x01; 0x02; 0x03]
logp3 = zeros(Float64, 3^m)
l = 0

for s1 in ss, s2 in ss, s3 in ss, s4 in ss
  l += 1
  S = [s1; s2; s3; s4]
  logp3[l] = logcondpostS(S, cluster, emptycluster, priorC.logγ, priorC.logω,
                          priorC.A, priorC.B)
end

M = maximum(logp3)
p3 = exp(M + log(sum(exp(logp3 - M))))

@test_approx_eq_eps p3 1.0 ε

# configuration with the highest posterior probability
Ctest = [0x02 0x02 0x01 0x01; 0x02 0x02 0x01 0x01]
N = 10000000
Sp = zeros(Float64, 3, m)
Cp = zero(Float64)

trueSp = hcat([0.388699761141992228274943954602349549531936645507812500
               0.404110835846893112766053945961175486445426940917968750
               0.207189403011114547936699636920820921659469604492187500],
              [0.388699761141992228274943954602349549531936645507812500
               0.404110835846893112766053945961175486445426940917968750
               0.207189403011114547936699636920820921659469604492187500],
              [0.768636427030134350424361855402821674942970275878906250
               0.231359112662804122795279226920683868229389190673828125
               0.000004460307061158063528846090539659030582697596400977],
              [0.745422155055375235122028243495151400566101074218750000
               0.254411468571633325730374508566455915570259094238281250
               0.000166376372991255944535152200280947454302804544568062])

for t in 1:N
  rpostpartitioncols!(C, cluster, emptycluster, priorC.logγ, priorC.logω,
                      priorC.A, priorC.B)

  for b in 1:m
    if C[cl[1], b] == 0x01
      @test C[cl[2], b] == 0x01

      Sp[1, b] += 1.0
    elseif C[cl[1], b] == 0x02
      @test C[cl[2], b] == 0x02

      Sp[2, b] += 1.0
    else
      Sp[3, b] += 1.0
    end
  end

  if C[cl, :] == Ctest
    Cp += 1.0
  end
end

Sp /= N
Cp /= N

C[cl, :] = Ctest
p4 = exp(logcondpostC(C, cluster, emptycluster, priorC.logγ, priorC.logω,
                      priorC.A, priorC.B))

@test maximum(abs(Sp - trueSp)) < 0.001
@test_approx_eq_eps Cp p4 0.001
