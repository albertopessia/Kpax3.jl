# This file is part of Kpax3. License is MIT.

function readresults(fp::IOStream)
  n = zeros(Int, 1)
  m = zeros(Int, 1)

  read!(fp, n)
  read!(fp, m)

  k = zeros(Int, 1)
  R = zeros(Int, n[1])

  absfreqk = zeros(Float64, n[1])
  absfreqR = zeros(Float64, div(n[1] * (n[1] - 1), 2))

  T = 0
  while !eof(fp)
    read!(fp, k)
    read!(fp, R)

    absfreqk[k[1]] += 1

    idx = 0
    for i in 1:(n[1] - 1), j in (i + 1):n[1]
      absfreqR[idx += 1] += Float64(R[i] == R[j])
    end

    T += 1
  end

  (absfreqk / T, absfreqR / T)
end

function saveresults!(fp::IOStream,
                      mcmcobj::AminoAcidMCMC)
  write(fp, length(mcmcobj.cl))
  write(fp, mcmcobj.R)
  #write(fp, vec(mcmcobj.C[!mcmcobj.emptycluster, :][1, :]))

  nothing
end
