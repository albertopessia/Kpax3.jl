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

  s = initializestate(data, D, kset, priorR, priorC, settings)

  if settings.verbose
    @printf("Initialization done. Creating solution set... ")
  end

  state = Array{AminoAcidState}(settings.popsize)
  logpp = zeros(Float64, settings.popsize)
  rank = Int[i for i in 1:settings.popsize]

  state[1] = AminoAcidState(data, s.R, priorR, priorC, settings)
  logpp[1] = state[1].logpp

  R = zeros(Int, n)
  for i in 2:settings.popsize
    copy!(R, s.R)
    modifypartition!(R, s.k)
    state[i] = AminoAcidState(data, R, priorR, priorC, settings)
    logpp[i] = state[i].logpp
  end

  sortperm!(rank, logpp, rev=true, initialized=true)

  if settings.verbose
    @printf("done.\n")
  end

  AminoAcidStateList(state, logpp, rank)
end

function AminoAcidStateList(data::Matrix{UInt8},
                            partition::Vector{Int},
                            priorR::PriorRowPartition,
                            priorC::AminoAcidPriorCol,
                            settings::KSettings)
  state = Array{AminoAcidState}(settings.popsize)
  logpp = zeros(Float64, settings.popsize)
  rank = Int[i for i in 1:settings.popsize]

  state[1] = AminoAcidState(data, partition, priorR, priorC, settings)
  logpp[1] = state[1].logpp

  if settings.verbose
    @printf("Creating solution set... ")
  end

  R = zeros(Int, size(data, 2))
  for i in 2:settings.popsize
    copy!(R, partition)
    modifypartition!(R, state[1].k)
    state[i] = AminoAcidState(data, R, priorR, priorC, settings)
    logpp[i] = state[i].logpp
  end

  sortperm!(rank, logpp, rev=true, initialized=true)

  if settings.verbose
    @printf("done.\n")
  end

  AminoAcidStateList(state, logpp, rank)
end

function AminoAcidStateList(popsize::Int,
                            s::AminoAcidState)
  state = Array{AminoAcidState}(popsize)
  logpp = zeros(Float64, popsize)
  rank = Int[i for i in 1:settings.popsize]

  for i in 1:popsize
    state[i] = copystate(s)
    logpp[i] = state[i].logpp
  end

  sortperm!(rank, logpp, rev=true, initialized=true)

  AminoAcidStateList(state, logpp, rank)
end

function copystatelist!(dest::AminoAcidStateList,
                        src::AminoAcidStateList,
                        popsize::Int)
  for i in 1:popsize
    copystate!(dest.state[i], src.state[src.rank[i]])
    dest.logpp[i] = src.logpp[src.rank[i]]
    dest.rank[i] = i
  end

  nothing
end
