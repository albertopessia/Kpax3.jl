# This file is part of Kpax3. License is MIT.

"""
# Genetic data

## Description

DNA data and its metadata.

## Fields

* `data` Multiple sequence alignment (MSA) encoded as a binary (UInt8) matrix
* `id` units' ids
* `ref` reference sequence, i.e. a vector of the same length of the original
sequences storing the values of homogeneous sites. SNPs are instead represented
by a value of 29
* `val` vector with unique values per MSA site
* `key` vector with indices of each value

## Details

Let `n` be the total number of units and `ml` be the total number of unique
values observed at SNP `l`. Define m = m1 + ... + mL, where L is the total
number os SNPs.

`data` is a `m`-by-`n` indicator matrix, i.e. `data[j, i]` is `1` if unit `i`
possesses value `j`, `0` otherwise.

The value associated with column `j` can be obtained by `val[j]` while the SNP
position by `find(ref == 29)[key[j]]`.

## References

Pessia A., Grad Y., Cobey S., Puranen J. S. and Corander J. (2015). K-Pax2:
Bayesian identification of cluster-defining amino acid positions in large
sequence datasets. _Microbial Genomics_ **1**(1).
<http://dx.doi.org/10.1099/mgen.0.000025>.
"""
immutable NucleotideData <: Kpax3Data
  data::Array{UInt8, 2}
  id::Array{ASCIIString, 1}
  ref::Array{UInt8, 1}
  val::Array{UInt8, 1}
  key::Array{Int, 1}
end

"""
# Genetic data

## Description

Amino acid data and its metadata.

## Fields

* `data` Multiple sequence alignment (MSA) encoded as a binary (UInt8) matrix
* `id` units' ids
* `ref` reference sequence, i.e. a vector of the same length of the original
sequences storing the values of homogeneous sites. SNPs are instead represented
by a value of 29
* `val` vector with unique values per MSA site
* `key` vector with indices of each value

## Details

Let `n` be the total number of units and `ml` be the total number of unique
values observed at SNP `l`. Define m = m1 + ... + mL, where L is the total
number os SNPs.

`data` is a `m`-by-`n` indicator matrix, i.e. `data[j, i]` is `1` if unit `i`
possesses value `j`, `0` otherwise.

The value associated with column `j` can be obtained by `val[j]` while the SNP
position by `find(ref == 29)[key[j]]`.

## References

Pessia A., Grad Y., Cobey S., Puranen J. S. and Corander J. (2015). K-Pax2:
Bayesian identification of cluster-defining amino acid positions in large
sequence datasets. _Microbial Genomics_ **1**(1).
<http://dx.doi.org/10.1099/mgen.0.000025>.
"""
immutable AminoAcidData <: Kpax3Data
  data::Array{UInt8, 2}
  id::Array{ASCIIString, 1}
  ref::Array{UInt8, 1}
  val::Array{UInt8, 1}
  key::Array{Int, 1}
end

# TODO: Julia v0.4
# we can't document this constructor: the type has already been documented

#=
# Constructor of an object of type NucleotideData

## Description

Create a new NucleotideData object.

## Usage

NucleotideData(f; miss=['?', '\*', '#', 'b', 'd', 'h', 'k', 'm', 'n', 'r', 's',
'v', 'w', 'x', 'y'], l=100000000, t=500)

## Arguments

* `infile` Path to the input data file
* `miss` Values to be considered missing.
* `l` Sequence length. If unknown, it is better to choose a value which is
surely greater than the real sequence length. If `l` is found to be
insufficient, the array size is dynamically increased (not recommended from a
performance point of view). Default value should be sufficient for most
datasets
* `t` If positive, print a status report every `t` read sequences
=#
function NucleotideData(infile::ASCIIString;
                        miss::Array{Char, 1}=['?', '*', '#', '-', 'b', 'd', 'h',
                                              'k', 'm', 'n', 'r', 's', 'v', 'w',
                                              'x', 'y'],
                        l::Int=100000000,
                        t::Int=500)
  missuint = zeros(UInt8, length(miss))
  i = 0

  for c in miss
    missuint[i += 1] = UInt8(c)
  end

  (data, id, ref) = readfasta(infile; dna=true, miss=missuint, l=l, t=t)
  (bindata, val, key) = categorical2binary(data)

  NucleotideData(bindata, id, ref, val, key)
end

#=
# Constructor of an object of type AminoAcidData

## Description

Create a new AminoAcidData object.

## Usage

AminoAcidData(infile; miss=['?', '\*', '#', 'b', 'j', 'x', 'z'],
l=100000000, t=500)

## Arguments

* `infile` Path to the input data file
* `miss` Values to be considered missing.
* `l` Sequence length. If unknown, it is better to choose a value which is
surely greater than the real sequence length. If `l` is found to be
insufficient, the array size is dynamically increased (not recommended from a
performance point of view). Default value should be sufficient for most
datasets
* `t` If positive, print a status report every `t` read sequences
=#
function AminoAcidData(infile::ASCIIString;
                       miss::Array{Char, 1}=['?', '*', '#', '-',
                                             'b', 'j', 'x', 'z'],
                       l::Int=100000000,
                       t::Int=500)

  missuint = zeros(UInt8, length(miss))
  i = 0

  for c in miss
    missuint[i += 1] = UInt8(c)
  end

  (data, id, ref) = readfasta(infile; dna=false, miss=missuint, l=l, t=t)
  (bindata, val, key) = categorical2binary(data)

  AminoAcidData(bindata, id, ref, val, key)
end
