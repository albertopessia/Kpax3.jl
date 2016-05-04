# This file is part of Kpax3. License is MIT.

abstract StateList

type AminoAcidStateList <: StateList
  state::Vector{AminoAcidState}
  logpp::Vector{Float64}
  rank::Vector{Int}
end

function AminoAcidStateList(data::Matrix{UInt8},
                            D::Matrix{Float64},
                            kset::UnitRange{Int},
                            priorR::PriorRowPartition,
                            priorC::AminoAcidPriorCol,
                            settings::KSettings)
  (m, n) = size(data)

  if settings.verbose
    @printf("Initializing state...\n")
  end

  (s, slp) = initializestate(data, D, kset, priorR, priorC, settings)

  if settings.verbose
    @printf("Initialization done. Creating solution set... ")
  end

  state = Array{AminoAcidState}(settings.popsize)
  state[1] = s

  logpp = zeros(Float64, settings.popsize)
  logpp[1] = slp

  R = zeros(Int, n)
  for i in 2:settings.popsize
    copy!(R, s.R)
    modifypartition!(R, s.k)
    state[i] = AminoAcidState(data, R, priorR, priorC, settings)
    logpp[i] = state[i].logpR + state[i].logpC[1] + state[i].loglik
  end

  rank = Int[i for i in 1:settings.popsize]
  sortperm!(rank, logpp, rev=true, initialized=true)

  if settings.verbose
    @printf("done.\n")
  end

  AminoAcidStateList(state, logpp, rank)
end

function AminoAcidStateList(data::Matrix{UInt8},
                            R::Vector{Int},
                            priorR::PriorRowPartition,
                            priorC::AminoAcidPriorCol,
                            settings::KSettings)
  state = Array{AminoAcidState}(settings.popsize)
  state[1] = AminoAcidState(data, R, priorR, priorC, settings)

  logpp = zeros(Float64, settings.popsize)
  logpp[1] = state[1].logpR + state[1].logpC[1] + state[1].loglik

  if settings.verbose
    @printf("Creating solution set... ")
  end

  R = zeros(Int, size(data, 2))
  for i in 2:settings.popsize
    copy!(R, state[1].R)
    modifypartition!(R, state[1].k)
    state[i] = AminoAcidState(data, R, priorR, priorC, settings)
    logpp[i] = state[i].logpR + state[i].logpC[1] + state[i].loglik
  end

  rank = Int[i for i in 1:settings.popsize]
  sortperm!(rank, logpp, rev=true, initialized=true)

  if settings.verbose
    @printf("done.\n")
  end

  AminoAcidStateList(state, logpp, rank)
end

function modifypartition!(R::Vector{Int},
                          k::Int)
  n = length(R)
  q = k

  if (k > 0) && (k < n + 1)
    q = rand(max(1, k - 10):min(n, k + 10))

    if q < k
      modifymerge!(R, k, q)
    elseif q > k
      modifysplit!(R, k, q)
    else
      modifyscramble!(R, k)
    end
  end

  q
end

function modifymerge!(R::Vector{Int},
                      k::Int,
                      q::Int)
  n = length(R)

  cset = zeros(Int, k)
  empty = trues(n)
  l = 1
  for a in 1:n
    if empty[R[a]]
      empty[R[a]] = false
      cset[l] = R[a]
      l += 1
    end
  end

  c = zeros(Int, 2)

  while k != q
    StatsBase.sample!(cset, c, replace=false, ordered=false)

    for a in 1:n
      if R[a] == c[2]
        R[a] = c[1]
      end
    end

    l = 1
    while cset[l] != c[2]
      l += 1
    end

    deleteat!(cset, l)
    k -= 1
  end

  nothing
end

function modifysplit!(R::Vector{Int},
                      k::Int,
                      q::Int)
  n = length(R)

  t = zeros(Int, n)
  for a in 1:n
    t[R[a]] += 1
  end

  g = 0

  while k != q
    k += 1

    w = StatsBase.WeightVec(Float64[t[a] > 0 ? t[a] - 1 : 0 for a in 1:n])
    g = StatsBase.sample(w)

    for a in 1:n
      if (R[a] == g) && ((t[k] == 0) || (rand() <= 0.25))
        R[a] = k
        t[g] -= 1
        t[k] += 1
      end
    end
  end

  nothing
end

function modifyscramble!(R::Vector{Int},
                         k::Int)
  n = length(R)

  t = zeros(Int, n)
  for a in 1:n
    t[R[a]] += 1
  end

  v = zeros(Int, n)
  moved = false

  l = 1
  g = 1
  h = 1
  while g < k
    if t[l] > 0
      t[l] = 0
      v[l] = 0

      for a in (l + 1):n
        v[a] = t[a] > 0 ? t[a] - 1 : 0
      end

      w = StatsBase.WeightVec(v)
      h = StatsBase.sample(w)

      keepgoing = true
      moved = false
      a = 1
      while keepgoing
        if R[a] == h
          if moved
            if rand() <= 0.05
              R[a] = g
              t[h] -= 1

              if t[h] == 1
                keepgoing = false
              end
            end
          else
            moved = true
            R[a] = g
            t[h] -= 1

            if t[h] == 1
              keepgoing = false
            end
          end
        end

        a += 1
        if a > n
          keepgoing = false
        end
      end

      g += 1
    end

    l += 1
  end

  nothing
end
