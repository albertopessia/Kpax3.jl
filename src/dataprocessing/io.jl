# This file is part of Kpax3. License is MIT.

function save(outfile::AbstractString,
              x::Kpax3Data)
  save(File(format"JLD", outfile),
       "data", x.data, "id", x.id, "ref", x.ref, "val", x.val, "key", x.key,
       compress=true)
end

function loadnt(infile::AbstractString)
  (data, id, ref, val, key) = load(infile, "data", "id", "ref", "val", "key")
  NucleotideData(data, id, ref, val, key)
end

function loadaa(infile::AbstractString)
  (data, id, ref, val, key) = load(infile, "data", "id", "ref", "val", "key")
  AminoAcidData(data, id, ref, val, key)
end
