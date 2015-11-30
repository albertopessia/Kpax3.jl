# This file is part of Kpax3. License is MIT.

ε = eps()

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
n1s = Float64[1, 5, 3, 3, 2, 2, 3, 3, 1, 3, 2, 1, 3, 2, 4, 2, 3, 3]

r1 = 2.0
r2 = 100.0

A1 = zeros(Float64, m, 4)
A1[:, 1] = (r1 + 1.0) * (n1s + 0.5) / (n + 1)
A1[:, 2] = 1.0
A1[:, 3] = 1.0
A1[:, 4] = r1

B1 = zeros(Float64, m, 4)
B1[:, 1] = (r1 + 1.0) - A1[:, 1]
B1[:, 2] = 1.0
B1[:, 3] = r1
B1[:, 4] = 1.0

A2 = zeros(Float64, m, 4)
A2[:, 1] = n1s + 0.5
A2[:, 2] = 1.0
A2[:, 3] = 1.0
A2[:, 4] = r2

B2 = zeros(Float64, m, 4)
B2[:, 1] = n - n1s + 0.5
B2[:, 2] = 1.0
B2[:, 3] = r2
B2[:, 4] = 1.0

for k in 1:n
  ω = [1.0, 1.0, 1.0 - 1.0 / k, 1.0 / k]

  for γ in ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
            [0.4, 0.3, 0.3], [0.5, 0.3, 0.2], [0.7, 0.2, 0.1],
            [0.1, 0.1, 0.1], [0.3, 0.1, 0.1], [0.0, 0.2, 0.1])
    x1 = AminoAcidPriorCol(data, k, γ, r1)
    x2 = AminoAcidPriorCol(data, k, γ, r2)

    for i in 1:3
      @test_approx_eq_eps x1.γ[i] (γ[i] / sum(γ)) ε
      @test_approx_eq_eps x1.logγ[i] log(γ[i] / sum(γ)) ε
      @test_approx_eq_eps x1.ω[i] ω[i] ε
      @test_approx_eq_eps x1.logω[i] log(ω[i]) ε

      @test_approx_eq_eps x2.γ[i] (γ[i] / sum(γ)) ε
      @test_approx_eq_eps x2.logγ[i] log(γ[i] / sum(γ)) ε
      @test_approx_eq_eps x2.ω[i] ω[i] ε
      @test_approx_eq_eps x2.logω[i] log(ω[i]) ε
    end

    for i in 1:4, j in 1:m
      @test_approx_eq_eps x1.A[j, i] A1[j, i] ε
      @test_approx_eq_eps x1.B[j, i] B1[j, i] ε

      @test_approx_eq_eps x2.A[j, i] A2[j, i] ε
      @test_approx_eq_eps x2.B[j, i] B2[j, i] ε
    end
  end
end
