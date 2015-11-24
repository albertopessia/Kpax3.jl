# This file is part of Kpax3. License is MIT.

include("data/partitions.jl")

@test_throws Kpax3DomainError EwensPitman(1.0, 0.0)
@test_throws Kpax3DomainError EwensPitman(2.0, 0.0)

@test_throws Kpax3DomainError EwensPitman(0.5, -1.0)
@test_throws Kpax3DomainError EwensPitman(0.5, -0.5)

@test_throws Kpax3DomainError EwensPitman(0.0, -1.0)
@test_throws Kpax3DomainError EwensPitman(0.0, 0.0)

@test_throws Kpax3DomainError EwensPitman(0.5, 0)
@test_throws Kpax3DomainError EwensPitman(0.0, 0)
@test_throws Kpax3DomainError EwensPitman(0.5, 1)
@test_throws Kpax3DomainError EwensPitman(0.0, 1)

@test_throws Kpax3DomainError EwensPitman(-1.0, 0)
@test_throws Kpax3DomainError EwensPitman(-1.0, -1)
@test_throws Kpax3DomainError EwensPitman(-1.0, 1.0)

for (α, θ) in ((0.4, -0.3), (0.4, 0.0), (0.4, 2.1), (0.0, 2.1), (-2.4, 3))
  ep = EwensPitman(α, θ)

  for i in 1:6
    po = TestPartition(i)
    pr = 0.0

    for j in 1:po.B
      pr += dEwensPitman(ep, po.partition[:, j])
    end

    @test_approx_eq_eps pr 1.0 3.0e-15

    pr = 0.0

    for j in 1:(po.C - 1)
      pr += (po.index[j + 1] - po.index[j]) * dEwensPitman(ep, i, po.k[j],
                                                           po.blocksize[:, j])
    end

    pr += dEwensPitman(ep, i, po.k[po.C], po.blocksize[:, po.C])

    @test_approx_eq_eps pr 1.0 3.0e-15
  end
end
