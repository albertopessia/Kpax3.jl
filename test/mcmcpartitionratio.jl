# This file is part of Kpax3. License is MIT.

ε = 1.0e-13

# suppose we split a cluster with 15 units into two clusters of 8 and 7 units
# respectively, moving from k = 5 to k = 6
n = 50
v = [22; 15; 7; 5; 1]
z = [22;  8; 7; 5; 1; 7]

ka = 5
kb = 6

vi = 8
vj = 7

for (α, θ) in ((0.4, -0.3), (0.4, 0.0), (0.4, 2.1), (0.0, 2.1), (-2.4, 10))
  ep = EwensPitman(α, θ)
  logratio = logdPriorRow(n, kb, z, ep) - logdPriorRow(n, ka, v, ep)
  @test_approx_eq_eps logratiopriorrowsplit(kb, vi, vj, ep) logratio ε
  @test_approx_eq_eps logratiopriorrowmerge(ka, vi, vj, ep) -logratio ε
end

# biased random walk
n = 50
v = [22; 15; 7; 5; 1]
z = [22; 14; 7; 5; 1; 1]

ka = 5
kb = 6

vi = 14
vj = 1

for (α, θ) in ((0.4, -0.3), (0.4, 0.0), (0.4, 2.1), (0.0, 2.1), (-2.4, 10))
  ep = EwensPitman(α, θ)
  logratio = logdPriorRow(n, kb, z, ep) - logdPriorRow(n, ka, v, ep)
  @test_approx_eq_eps logratiopriorrowbrwsplit(kb, vi + vj, ep) logratio ε
  @test_approx_eq_eps logratiopriorrowbrwmerge(ka, vi, ep) -logratio ε
end
