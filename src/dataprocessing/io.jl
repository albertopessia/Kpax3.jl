# This file is part of Kpax3. License is MIT.

function save(outfile::AbstractString,
              x::KData)
  # create directory if it does not exist
  dirpath = dirname(outfile)
  if !isdir(dirpath)
    mkpath(dirpath)
  end

  save(File(format"JLD", outfile),
       "data", x.data, "id", x.id, "ref", x.ref, "val", x.val, "key", x.key,
       compress=true)
end

function loadnt(infile::AbstractString)
  # open infile for reading and immediately close it. We do this to throw a
  # proper Julia standard exception if something is wrong
  f = open(infile, "r")
  close(f)

  (data, id, ref, val, key) = load(infile, "data", "id", "ref", "val", "key")
  NucleotideData(data, id, ref, val, key)
end

function loadaa(infile::AbstractString)
  # open infile for reading and immediately close it. We do this to throw a
  # proper Julia standard exception if something is wrong
  f = open(infile, "r")
  close(f)

  (data, id, ref, val, key) = load(infile, "data", "id", "ref", "val", "key")
  AminoAcidData(data, id, ref, val, key)
end
