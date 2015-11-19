# This file is part of K-Pax3. License is MIT.

function save(outfile::AbstractString,
              x::Kpax3Data,
              varname::AbstractString="kpax3data")
  save(File(format"JLD", outfile), varname , x)
end
