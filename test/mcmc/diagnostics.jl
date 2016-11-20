# This file is part of Kpax3. License is MIT.

# these tests must be run after "mcmc/posterior.jl"
function test_traceR()
  fileroot = "../build/mcmc_6"
  maxlag = 20

  (entropy, avgd) = traceR(fileroot, maxlag=maxlag)

  fpR = open(string(fileroot, "_row_partition.bin"), "r")

  tmp = zeros(Int, 1)

  read!(fpR, tmp)
  n = tmp[1]

  read!(fpR, tmp)
  read!(fpR, tmp)
  N = tmp[1]

  k = zeros(Int, N)
  R = zeros(Int, n, N)
  v = zeros(Int, n, N)

  e = 0.0

  k0 = zeros(Int, 1)
  R0 = zeros(Int, n)

  T = 1
  while !eof(fpR)
    read!(fpR, k0)
    read!(fpR, R0)

    k[T] = k0[1]
    copy!(R, 1 + n * (T - 1), normalizepartition(R0, n), 1, n)

    for i in 1:n
      v[R[i, T], T] += 1
    end

    e = 0.0
    for g in 1:k[T]
      e -= v[g, T] * (log(v[g, T]) - log(n)) / n
    end

    @test_approx_eq_eps entropy[T] e ε

    T += 1
  end

  close(fpR)

  # test lag 1
  l = 1
  z = zeros(Float64, N - l)
  for t in 1:(N - l)
    z[t] = distsj(R[:, t], v[:, t], k[t], R[:, t+l], v[:, t+l], k[t+l], n)
  end
  @test_approx_eq_eps avgd[l] mean(z) ε

  # test lag 5
  l = 5
  z = zeros(Float64, N - l)
  for t in 1:(N - l)
    z[t] = distsj(R[:, t], v[:, t], k[t], R[:, t+l], v[:, t+l], k[t+l], n)
  end
  @test_approx_eq_eps avgd[l] mean(z) ε

  # test lag 13
  l = 13
  z = zeros(Float64, N - l)
  for t in 1:(N - l)
    z[t] = distsj(R[:, t], v[:, t], k[t], R[:, t+l], v[:, t+l], k[t+l], n)
  end
  @test_approx_eq_eps avgd[l] mean(z) ε

  nothing
end

test_traceR()

function test_traceC()
  fileroot = "../build/mcmc_6"
  maxlag = 20

  (entropy, avgd) = traceC(fileroot, maxlag=maxlag)

  fpC = open(string(fileroot, "_col_partition.bin"), "r")

  tmp = zeros(Int, 1)

  read!(fpC, tmp)
  n = tmp[1]

  read!(fpC, tmp)
  m = tmp[1]

  read!(fpC, tmp)
  N = tmp[1]

  C = zeros(UInt8, m, N)
  e = 0.0

  C0 = zeros(UInt8, m)
  v0 = zeros(Float64, 3)

  T = 1
  while !eof(fpC)
    readbytes!(fpC, C0, m)
    copy!(C, 1 + m * (T - 1), C0, 1, m)

    fill!(v0, 0.0)
    for b in 1:m
      v0[C[b, T]] += 1
    end
    v0 /= m

    e = 0.0
    if v0[1] > 0.0
      e -= v0[1] * log(v0[1])
    end

    if v0[2] > 0.0
      e -= v0[2] * log(v0[2])
    end

    if v0[3] > 0.0
      e -= v0[3] * log(v0[3])
    end

    @test_approx_eq_eps entropy[T] e ε

    T += 1
  end

  close(fpC)

  # test lag 1
  l = 1
  z = zeros(Float64, N - l)
  for t in 1:(N - l)
    z[t] = Distances.hamming(C[:, t], C[:, t + l]) / m
  end
  @test_approx_eq_eps avgd[l] mean(z) ε

  # test lag 5
  l = 5
  z = zeros(Float64, N - l)
  for t in 1:(N - l)
    z[t] = Distances.hamming(C[:, t], C[:, t + l]) / m
  end

  @test_approx_eq_eps avgd[l] mean(z) ε

  # test lag 13
  l = 13
  z = zeros(Float64, N - l)
  for t in 1:(N - l)
    z[t] = Distances.hamming(C[:, t], C[:, t + l]) / m
  end

  @test_approx_eq_eps avgd[l] mean(z) ε

  nothing
end

test_traceC()
