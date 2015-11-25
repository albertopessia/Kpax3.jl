# This file is part of Kpax3. License is MIT.

ε = 1.0e-13

for (α, β) in ([0.1, 0.1], [0.5, 0.5], [1.0, 1.0], [10.0, 10.0], [100.0, 100.0],
               [0.2, 0.1], [1.0, 0.5], [2.0, 1.0], [20.0, 10.0], [200.0, 100.0],
               [0.1, 0.2], [0.5, 1.0], [1.0, 2.0], [10.0, 20.0], [100.0, 200.0])
  for (n, x) in ([  1.0, 0.0], [  1.0, 1.0],
                 [  5.0, 0.0], [  5.0, 1.0], [  5.0,  3.0], [   5.0,  5.0],
                 [ 10.0, 0.0], [ 10.0, 1.0], [ 10.0,  5.0], [  10.0, 10.0],
                 [100.0, 0.0], [100.0, 1.0], [100.0, 10.0], [100.0, 100.0])
    logp = lgamma(α + x) + lgamma(β + n - x) - lgamma(α + β + n) +
           lgamma(α + β) - lgamma(α) - lgamma(β)

    @test_approx_eq_eps marglik(x, n, α, β) exp(logp) ε
    @test_approx_eq_eps logmarglik(x, n, α, β) logp ε
  end
end

for (α, β) in ([0.1, 0.1], [0.5, 0.5], [1.0, 1.0], [10.0, 10.0], [100.0, 100.0],
               [0.2, 0.1], [1.0, 0.5], [2.0, 1.0], [20.0, 10.0], [200.0, 100.0],
               [0.1, 0.2], [0.5, 1.0], [1.0, 2.0], [10.0, 20.0], [100.0, 200.0])
  for (n, x) in ((  1.0, [  0.0, 0.0,  1.0,  0.0, 1.0]),
                 ( 10.0, [  2.0, 8.0,  5.0, 10.0, 0.0]),
                 (100.0, [100.0, 2.0, 72.0, 34.0, 0.0]))
    m = length(x)

    a = zeros(Float64, m)
    b = zeros(Float64, m)

    a[:] = α
    b[:] = β

    q = marglik(x, n, a, b)
    logq = logmarglik(x, n, a, b)

    logp = lgamma(α + x) + lgamma(β + n - x) - lgamma(α + β + n) +
           lgamma(α + β) - lgamma(α) - lgamma(β)

    for i in 1:m
      @test_approx_eq_eps q[i] exp(logp[i]) ε
      @test_approx_eq_eps logq[i] logp[i] ε
    end
  end
end
