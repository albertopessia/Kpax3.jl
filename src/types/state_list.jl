# This file is part of Kpax3. License is MIT.

abstract StateList

type AminoAcidStateList <: StateList
  state::Array{AminoAcidState}
  logpp::Array{Float64}
  rank::Array{Int}
end

function AminoAcidStateList(x::AminoAcidData,
                            d::Vector{Float64},
                            N::Int,
                            settings::KSettings;
                            kset::UnitRange{Int}=1:0)
  m, n = size(x.data)

  if length(d) != div(n * (n - 1), 2)
    throw(KInputError(string("Argument 'd' does not have the correct length. ",
                             "Expecting a length of ", div(n * (n - 1), 2),
                             " but found ", length(d), " instead.")))
  end

  if N < 2
    throw(KDomainError("Argument 'N' is less than 2."))
  end

  if length(kset) == 0
    kset = 2:max(ceil(Int, sqrt(n)), 100)
  elseif kset[1] > 0
    if kset[end] > n
      kset = (kset[1] == 1) ? 2:n : kset[1]:n
    elseif kset[1] == 1
      kset = 2:kset[end]
    end
  else
    throw(KDomainError("First element of 'kset' is less than one."))
  end

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(x.data, 1, settings.γ, settings.r)

  R = ones(Int, n)

  s = AminoAcidState(x.data, R, priorR, priorC, settings)
  slp = s.logpR + s.logpC[1] + s.loglik

  D = zeros(Float64, n, n)
  idx = 1
  for j in 1:(n - 1), i in (j + 1):n
    D[i, j] = D[j, i] = d[idx]
    idx += 1
  end

  t1 = copystate(s)
  tlp1 = slp

  t2 = copystate(s)
  tlp2 = slp

  niter = 0
  for k in kset
    updateprior!(priorC, k)

    copy!(R, kmedoids(D, k).assignments)
    updatestate!(t1, x.data, R, priorR, priorC, settings)
    tlp1 = t1.logpR + t1.logpC[1] + t1.loglik

    niter = 0
    while niter < 10
      copy!(R, kmedoids(D, k).assignments)
      updatestate!(t2, x.data, R, priorR, priorC, settings)
      tlp2 = t2.logpR + t2.logpC[1] + t2.loglik

      if tlp2 > tlp1
        copystate!(t1, t2)
        tlp1 = tlp2
      end

      niter += 1
    end

    if tlp1 > slp
      copystate!(s, t1)
      slp = tlp1
    end
  end

  state = Array{AminoAcidState}(N)
  state[1] = s

  logpp = zeros(Float64, N)
  logpp[1] = slp

  for i in 2:N
    copy!(R, s.R)
    q = modifypartition!(R, s.k)
    updateprior!(priorC, q)
    state[i] = AminoAcidState(x.data, R, priorR, priorC, settings)
    logpp[i] = state[i].logpR + state[i].logpC[1] + state[i].loglik
  end

  rank = Int[i for i in 1:N]
  sortperm!(rank, logpp, rev=true, initialized=true)

  AminoAcidStateList(state, logpp, rank)
end

function modifypartition!(R::Vector{Int},
                          k::Int)
  n = length(R)
  q = k

  if (k > 1) && (k < n)
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
