# This file is part of Kpax3. License is MIT.

__precompile__(true)

module Kpax3
  #################
  # Load packages #
  #################
  import StatsBase: sample, sample!, WeightVec, values
  import Clustering: kmedoids
  import Distances: hamming
  import FileIO: File, @format_str
  import JLD: load, save
  import Gadfly: plot, layer, get_stroke_vector, Coord, Geom, Guide, Scale,
                 Theme, mm, @colorant_str, draw, SVG, SVGJS, PNG, PDF, PS, PGF

  ####################
  # Export functions #
  ####################
  export

  # types
  KSettings,
  AminoAcidData,
  AminoAcidState,

  # distances
  distaamtn84,

  # algorithms
  kpax3mcmc,
  kpax3ga,

  # estimation
  processchain,
  kpax3estimate,
  computelocalmode!,

  # i/o
  initializepartition,
  normalizepartition,
  readfasta,
  categorical2binary,
  save,
  loadaa,
  writeresults,
  readposteriork,
  readposteriorP,
  readposteriorC,

  # diagnostics
  traceR,
  traceC,

  # plot
  plotk,
  plotP,
  plotC

  ##############
  # Load Types #
  ##############
  include("types/types.jl")

  ########################
  # Load basic functions #
  ########################
  include("misc/misc.jl")

  #################
  # Load the rest #
  #################
  include("data_processing/data_processing.jl")
  include("distances/distances.jl")
  include("model/model.jl")
  include("optimizer/optimizer.jl")
  include("mcmc/mcmc.jl")
  include("estimate/estimate.jl")
  include("plots/plots.jl")
end
