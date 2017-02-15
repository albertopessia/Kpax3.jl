# This file is part of Kpax3. License is MIT.

function test_mcmc_clusterweight()
  n = 10
  k = 2

  priorR = Kpax3.EwensPitman(0.5, -0.1)
  @test_approx_eq_eps Kpax3.clusterweight(n, priorR) log(n - 0.5) ε

  priorR = Kpax3.EwensPitman(0.5, 0.0)
  @test_approx_eq_eps Kpax3.clusterweight(n, priorR) log(n - 0.5) ε

  priorR = Kpax3.EwensPitman(0, 2.0)
  @test_approx_eq_eps Kpax3.clusterweight(n, priorR) log(n) ε

  priorR = Kpax3.EwensPitman(-2, 5)
  @test_approx_eq_eps Kpax3.clusterweight(n, priorR) log(n + 2.0) ε

  priorR = Kpax3.EwensPitman(0.5, -0.1)
  @test_approx_eq_eps Kpax3.clusterweight(n, k, priorR) log(0.5 * k - 0.1) ε

  priorR = Kpax3.EwensPitman(0.5, 0.0)
  @test_approx_eq_eps Kpax3.clusterweight(n, k, priorR) log(0.5 * k) ε

  priorR = Kpax3.EwensPitman(0, 2.0)
  @test_approx_eq_eps Kpax3.clusterweight(n, k, priorR) log(2) ε

  priorR = Kpax3.EwensPitman(-2, 5)
  @test_approx_eq_eps Kpax3.clusterweight(n, k, priorR) log(-2.0 * (k - 5)) ε

  nothing
end

test_mcmc_clusterweight()
