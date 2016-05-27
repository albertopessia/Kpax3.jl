# This file is part of Kpax3. License is MIT.

function test_prior_row_exceptions()
  @test_throws KDomainError EwensPitman(1.0, 0.0)
  @test_throws KDomainError EwensPitman(2.0, 0.0)

  @test_throws KDomainError EwensPitman(0.5, -1.0)
  @test_throws KDomainError EwensPitman(0.5, -0.5)

  @test_throws KDomainError EwensPitman(0.0, -1.0)
  @test_throws KDomainError EwensPitman(0.0, 0.0)

  @test_throws KDomainError EwensPitman(0.5, 0)
  @test_throws KDomainError EwensPitman(0.0, 0)
  @test_throws KDomainError EwensPitman(0.5, 1)
  @test_throws KDomainError EwensPitman(0.0, 1)

  @test_throws KDomainError EwensPitman(-1.0, 0)
  @test_throws KDomainError EwensPitman(-1.0, -1)
  @test_throws KDomainError EwensPitman(-1.0, 1.0)

  nothing
end

test_prior_row_exceptions()

function test_prior_row_constructor()
  priorR = EwensPitman(0.5, -0.1)
  @test isa(priorR, EwensPitmanPAUT)
  @test priorR.α == 0.5
  @test priorR.θ == -0.1

  priorR = EwensPitman(0.5, 0.0)
  @test isa(priorR, EwensPitmanPAZT)
  @test priorR.α == 0.5

  priorR = EwensPitman(0, 2.0)
  @test isa(priorR, EwensPitmanZAPT)
  @test priorR.θ == 2.0

  priorR = EwensPitman(-2, 5)
  @test isa(priorR, EwensPitmanNAPT)
  @test priorR.α == -2.0
  @test priorR.L == 5

  nothing
end

test_prior_row_constructor()
