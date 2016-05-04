# This file is part of Kpax3. License is MIT.

function save(ofile::AbstractString,
              x::KData)
  # create directory if it does not exist
  dirpath = dirname(ofile)
  if !isdir(dirpath)
    mkpath(dirpath)
  end

  save(File(format"JLD", ofile),
       "data", x.data, "id", x.id, "ref", x.ref, "val", x.val, "key", x.key,
       compress=true)
end

function loadnt(ifile::AbstractString)
  # open ifile for reading and immediately close it. We do this to throw a
  # proper Julia standard exception if something is wrong
  f = open(ifile, "r")
  close(f)

  (data, id, ref, val, key) = load(ifile, "data", "id", "ref", "val", "key")
  NucleotideData(data, id, ref, val, key)
end

function loadaa(ifile::AbstractString)
  # open ifile for reading and immediately close it. We do this to throw a
  # proper Julia standard exception if something is wrong
  f = open(ifile, "r")
  close(f)

  (data, id, ref, val, key) = load(ifile, "data", "id", "ref", "val", "key")
  AminoAcidData(data, id, ref, val, key)
end
