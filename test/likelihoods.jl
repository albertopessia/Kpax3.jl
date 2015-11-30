# This file is part of Kpax3. License is MIT.

ε = 1.0e-13

for (α, β) in ([0.1, 0.1], [0.5, 0.5], [1.0, 1.0], [10.0, 10.0], [100.0, 100.0],
               [0.2, 0.1], [1.0, 0.5], [2.0, 1.0], [20.0, 10.0], [200.0, 100.0],
               [0.1, 0.2], [0.5, 1.0], [1.0, 2.0], [10.0, 20.0], [100.0, 200.0])
  for (n, y) in ([  1.0, 0.0], [  1.0, 1.0],
                 [  5.0, 0.0], [  5.0, 1.0], [  5.0,  3.0], [   5.0,  5.0],
                 [ 10.0, 0.0], [ 10.0, 1.0], [ 10.0,  5.0], [  10.0, 10.0],
                 [100.0, 0.0], [100.0, 1.0], [100.0, 10.0], [100.0, 100.0])
    logp = lgamma(α + y) + lgamma(β + n - y) - lgamma(α + β + n) +
           lgamma(α + β) - lgamma(α) - lgamma(β)

    logcp0 = log(β + n - y) - log(α + β + n)
    logcp1 = log(α + y) - log(α + β + n)

    @test_approx_eq_eps marglik(y, n, α, β) exp(logp) ε
    @test_approx_eq_eps logmarglik(y, n, α, β) logp ε

    @test_approx_eq_eps condmarglik(0x00, y, n, α, β) exp(logcp0) ε
    @test_approx_eq_eps condmarglik(0x01, y, n, α, β) exp(logcp1) ε

    @test_approx_eq_eps logcondmarglik(0x00, y, n, α, β) logcp0 ε
    @test_approx_eq_eps logcondmarglik(0x01, y, n, α, β) logcp1 ε
  end
end

for (α, β) in ([0.1, 0.1], [0.5, 0.5], [1.0, 1.0], [10.0, 10.0], [100.0, 100.0],
               [0.2, 0.1], [1.0, 0.5], [2.0, 1.0], [20.0, 10.0], [200.0, 100.0],
               [0.1, 0.2], [0.5, 1.0], [1.0, 2.0], [10.0, 20.0], [100.0, 200.0])
  for (n, y) in ((  1.0, [  0.0, 0.0,  1.0,  0.0, 1.0]),
                 ( 10.0, [  2.0, 8.0,  5.0, 10.0, 0.0]),
                 (100.0, [100.0, 2.0, 72.0, 34.0, 0.0]))
    m = length(y)
    n1s = cld(m, 2)

    a = fill(α, m)
    b = fill(β, m)

    q = marglik(y, n, a, b)
    logq = logmarglik(y, n, a, b)

    q0 = condmarglik(fill(0x00, m), y, n, a, b)
    logq0 = logcondmarglik(fill(0x00, m), y, n, a, b)

    q1 = condmarglik(fill(0x01, m), y, n, a, b)
    logq1 = logcondmarglik(fill(0x01, m), y, n, a, b)

    qx = condmarglik([fill(0x01, n1s); fill(0x00, m - n1s)], y, n, a, b)
    logqx = logcondmarglik([fill(0x01, n1s); fill(0x00, m - n1s)], y, n, a, b)

    logp = lgamma(a + y) + lgamma(b + n - y) - lgamma(a + b + n) +
           lgamma(a + b) - lgamma(a) - lgamma(b)

    logcp0 = log(b + n - y) - log(a + b + n)
    logcp1 = log(a + y) - log(a + b + n)

    logcpx = log([a[1:n1s] + y[1:n1s]; b[(n1s + 1):m] + n - y[(n1s + 1):m]]) -
             log(a + b + n)

    for i in 1:m
      @test_approx_eq_eps q[i] exp(logp[i]) ε
      @test_approx_eq_eps logq[i] logp[i] ε

      @test_approx_eq_eps q0[i] exp(logcp0[i]) ε
      @test_approx_eq_eps logq0[i] logcp0[i] ε

      @test_approx_eq_eps q1[i] exp(logcp1[i]) ε
      @test_approx_eq_eps logq1[i] logcp1[i] ε

      @test_approx_eq_eps qx[i] exp(logcpx[i]) ε
      @test_approx_eq_eps logqx[i] logcpx[i] ε
    end
  end
end
