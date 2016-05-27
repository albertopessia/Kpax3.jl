# This file is part of Kpax3. License is MIT.

function processchain(ifile::AbstractString)
  fp = open(ifile, "r")

  # discard prior hyperparameters
  tmp = zeros(Float64, 6)
  read!(fp, tmp)

  n = zeros(Int, 1)
  m = zeros(Int, 1)

  read!(fp, n)
  read!(fp, m)

  k = zeros(Int, 1)
  R = zeros(Int, n[1])
  C = zeros(UInt8, m[1])

  absfreqk = zeros(Float64, n[1])
  absfreqR = zeros(Float64, div(n[1] * (n[1] - 1), 2))
  absfreqS = zeros(Float64, 3, m[1])

  T = 0
  while !eof(fp)
    read!(fp, k)
    read!(fp, R)
    readbytes!(fp, C, m[1])

    absfreqk[k[1]] += 1

    idx = 0
    for i in 1:(n[1] - 1), j in (i + 1):n[1]
      idx += 1
      absfreqR[idx] += float(R[i] == R[j])
    end

    for b in 1:m[1]
      absfreqS[C[b], b] += 1
    end

    T += 1
  end

  close(fp)

  (absfreqk / T, absfreqR / T, absfreqS / T)
end

function savestate!(fp::IOStream,
                    state::AminoAcidState)
  write(fp, state.k)
  write(fp, state.R)
  for b in 1:size(state.C, 2)
    write(fp, state.C[state.cl[1], b] < 0x03 ? state.C[state.cl[1], b] : 0x03)
  end
  nothing
end
