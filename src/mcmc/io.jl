# This file is part of Kpax3. License is MIT.

function readresults(infile::AbstractString)
  fp = open(infile, "r")

  n = zeros(Int, 1)
  m = zeros(Int, 1)

  read!(fp, n)
  read!(fp, m)

  k = zeros(Int, 1)
  R = zeros(Int, n[1])
  C = zeros(UInt8, n[1] * m[1])

  absfreqk = zeros(Float64, n[1])
  absfreqR = zeros(Float64, div(n[1] * (n[1] - 1), 2))
  absfreqS = zeros(Float64, 4, m[1])

  T = 0
  while !eof(fp)
    read!(fp, k)
    read!(fp, R)
    readbytes!(fp, C, k[1] * m[1])

    absfreqk[k[1]] += 1

    idx = 0
    for i in 1:(n[1] - 1), j in (i + 1):n[1]
      absfreqR[idx += 1] += float(R[i] == R[j])
    end

    for b in 1:m[1]
      absfreqS[C[1 + k[1] * (b - 1)], b] += 1
    end

    T += 1
  end

  close(fp)

  for b in 1:m[1]
    absfreqS[3, b] += absfreqS[4, b]
  end

  (absfreqk / T, absfreqR / T, absfreqS[1:3, :] / T)
end

function saveresults!(fp::IOStream,
                      mcmcobj::AminoAcidMCMC)
  write(fp, mcmcobj.k)
  write(fp, mcmcobj.R)
  write(fp, mcmcobj.C[mcmcobj.cl[1:mcmcobj.k], :])
  nothing
end
