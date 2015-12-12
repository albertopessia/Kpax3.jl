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

R = [1; 1; 1; 2; 2; 3]

k = maximum(R)

g = [0.6; 0.35; 0.05]
r = log(0.001) / log(0.95)
priorC = AminoAcidPriorCol(data, k, g, r)

C = zeros(UInt8, maxclust, m)

cl = [13; 42; 76]

v = zeros(Int, maxclust)
n1s = zeros(Float64, maxclust, m)

v[cl[1]] = sum(R .== 1)
n1s[cl[1], :] = copy(sum(float(data[:, R .== 1]), 2))

v[cl[2]] = sum(R .== 2)
n1s[cl[2], :] = copy(sum(float(data[:, R .== 2]), 2))

v[cl[3]] = sum(R .== 3)
n1s[cl[3], :] = copy(sum(float(data[:, R .== 3]), 2))

cs = ([0x01; 0x01; 0x01], [0x02; 0x02; 0x02], [0x03; 0x03; 0x03],
      [0x03; 0x03; 0x04], [0x03; 0x04; 0x03], [0x04; 0x03; 0x03],
      [0x03; 0x04; 0x04], [0x04; 0x03; 0x04], [0x04; 0x04; 0x03],
      [0x04; 0x04; 0x04])

logp1 = zeros(Float64, (2 + 2^k)^m)
logp2 = zeros(Float64, (2 + 2^k)^m)

l = 0

for c1 in cs, c2 in cs, c3 in cs, c4 in cs
  l += 1

  tmp = hcat(c1, c2, c3, c4)

  C[cl[1], :] = tmp[1, :]
  C[cl[2], :] = tmp[2, :]
  C[cl[3], :] = tmp[3, :]

  logp1[l] = logpriorC(C, cl, priorC.logγ, priorC.logω)
  logp2[l] = logcondpostC(C, cl, v, n1s, priorC.logγ, priorC.logω, priorC.A,
                          priorC.B)
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
  logp3[l] = logcondpostS(S, cl, v, n1s, priorC.logγ, priorC.logω, priorC.A,
                          priorC.B)
end

M = maximum(logp3)
p3 = exp(M + log(sum(exp(logp3 - M))))

@test_approx_eq_eps p3 1.0 ε

# configuration with the highest posterior probability
Ctest = [0x01 0x01 0x01 0x01; 0x01 0x01 0x01 0x01; 0x01 0x01 0x01 0x01]
N = 10000000
Sp = zeros(Float64, 3, m)
Cp = zero(Float64)

trueSp = hcat([0.650881788845628306283686015376588329672813415527343750000000
               0.345392727967450685611083827097900211811065673828125000000000
               0.003725483186920867158947734409935037547256797552108764648438],
              [0.650854890249376039079720612789969891309738159179687500000000
               0.345378454132021173172972794418456032872200012207031250000000
               0.003766655618603001551975006933048462087754160165786743164063],
              [0.746989515332794340451982861850410699844360351562500000000000
               0.252948830588882567216302277302020229399204254150390625000000
               0.000061654078322953418311511142313463551545282825827598571777],
              [0.730447714956305782507683943549636751413345336914062500000000
               0.266613264143007511197680514669627882540225982666015625000000
               0.002939020900686826857917122168828427675180137157440185546875])

for t in 1:N
  rpostpartitioncols!(C, cl, v, n1s, priorC.logγ, priorC.logω, priorC.A,
                      priorC.B)

  for b in 1:m
    if C[cl[1], b] == 0x01
      @test C[cl[2], b] == 0x01
      @test C[cl[3], b] == 0x01

      Sp[1, b] += 1.0
    elseif C[cl[1], b] == 0x02
      @test C[cl[2], b] == 0x02
      @test C[cl[3], b] == 0x02

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
p4 = exp(logcondpostC(C, cl, v, n1s, priorC.logγ, priorC.logω, priorC.A,
                      priorC.B))

@test maximum(abs(Sp - trueSp)) < 0.001
@test_approx_eq_eps Cp p4 0.001

logpr, logpo = rpostpartitioncols!(C, cl, v, n1s, priorC.logγ, priorC.logω,
                                   priorC.A, priorC.B)

@test_approx_eq_eps logpr logpriorC(C, cl, priorC.logγ, priorC.logω) ε
@test_approx_eq_eps logpo logcondpostC(C, cl, v, n1s, priorC.logγ, priorC.logω,
                                       priorC.A, priorC.B) ε
