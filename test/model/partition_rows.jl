# This file is part of Kpax3. License is MIT.

function test_partition_rows_functions()
  include("../data/partitions.jl")

  for (α, θ) in ((0.4, -0.3), (0.4, 0.0), (0.4, 2.1), (0.0, 2.1), (-2.4, 3))
    ep = EwensPitman(α, θ)

    for i in 1:6
      po = TestPartition(i)

      pr = 0.0
      qr = 0.0

      for j in 1:po.B
        pr += dPriorRow(po.partition[:, j], ep)
        qr += exp(logdPriorRow(po.partition[:, j], ep))
      end

      @test_approx_eq_eps pr 1.0 ε
      @test_approx_eq_eps qr 1.0 ε

      pr = 0.0
      qr = 0.0

      for j in 1:(po.C - 1)
        pr += (po.index[j + 1] - po.index[j]) *
               dPriorRow(i, po.k[j], po.blocksize[:, j], ep)
        qr += (po.index[j + 1] - po.index[j]) *
               exp(logdPriorRow(i, po.k[j], po.blocksize[:, j], ep))
      end

      pr += dPriorRow(i, po.k[po.C], po.blocksize[:, po.C], ep)
      qr += exp(logdPriorRow(i, po.k[po.C], po.blocksize[:, po.C], ep))

      @test_approx_eq_eps pr 1.0 ε
      @test_approx_eq_eps qr 1.0 ε
    end
  end

  nothing
end

test_partition_rows_functions()
