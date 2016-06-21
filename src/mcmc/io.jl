# This file is part of Kpax3. License is MIT.

function processchain(ifile::AbstractString)
  n = zeros(Int, 1)
  m = zeros(Int, 1)
  N = zeros(Int, 1)

  k = zeros(Int, 1)

  h = 0

  fpR = open(string(ifile, "_row_partition.bin"), "r")

  read!(fpR, n)
  read!(fpR, m)
  read!(fpR, N)

  R = zeros(Int, n[1])

  absfreqk = zeros(Float64, n[1])
  absfreqR = zeros(Float64, div(n[1] * (n[1] - 1), 2))

  T = 0
  while !eof(fpR)
    read!(fpR, k)
    read!(fpR, R)

    absfreqk[k[1]] += 1

    h = 0
    for i in 1:(n[1] - 1), j in (i + 1):n[1]
      h += 1
      absfreqR[h] += float(R[i] == R[j])
    end

    T += 1
  end

  close(fpR)

  if T != N[1]
    throw(KInputError(string("Expecting ", N[1], " simulations but found ", T)))
  end

  fpC = open(string(ifile, "_col_partition.bin"), "r")

  read!(fpC, n)
  read!(fpC, m)
  read!(fpC, N)

  C = zeros(UInt8, m[1])

  absfreqS = zeros(Float64, 3, m[1])

  T = 0
  while !eof(fpC)
    read!(fpC, k)
    readbytes!(fpC, C, m[1])

    for b in 1:m[1]
      absfreqS[C[b], b] += 1
    end

    T += 1
  end

  close(fpC)

  if T != N[1]
    throw(KInputError(string("Expecting ", N[1], " simulations but found ", T)))
  end

  (absfreqk / T, absfreqR / T, absfreqS / T)
end

function savestate!(fpR::IOStream,
                    fpC::IOStream,
                    state::AminoAcidState)
  write(fpR, state.k)
  write(fpR, state.R)

  write(fpC, state.k)
  for b in 1:size(state.C, 2)
    write(fpC, state.C[state.cl[1], b] < 0x03 ? state.C[state.cl[1], b] : 0x03)
  end

  nothing
end
