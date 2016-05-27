# This file is part of Kpax3. License is MIT.

function test_mcmc_splitmergeweight()
  n = 10

  priorR = EwensPitman(0.5, -0.1)
  @test_approx_eq_eps splitmergeweight(n, priorR) log(n - 0.5) ε

  priorR = EwensPitman(0.5, 0.0)
  @test_approx_eq_eps splitmergeweight(n, priorR) log(n - 0.5) ε

  priorR = EwensPitman(0, 2.0)
  @test_approx_eq_eps splitmergeweight(n, priorR) log(n) ε

  priorR = EwensPitman(-2, 5)
  @test_approx_eq_eps splitmergeweight(n, priorR) log(n + 2.0) ε

  nothing
end

test_mcmc_splitmergeweight()
