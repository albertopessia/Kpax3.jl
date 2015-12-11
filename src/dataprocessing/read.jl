# This file is part of Kpax3. License is MIT.

"""
# Read a FASTA formatted file

## Description

Read data in FASTA format and convert it to an integer matrix. Sequences are
required to be aligned. Only polymorphic columns are stored.

## Usage

readfasta(infile, dna, miss, l, verbose, verbosestep)

## Arguments

* `infile` Path to the input data file
* `dna` `true` if reading DNA data or `false` if reading protein data
* `miss` Characters (as `UInt8`) to be considered missing. Use
`miss = zeros(UInt8, 1)` if all characters are to be considered valid. Default
characters for `miss` are:

  - DNA data: _?, \*, #, -, b, d, h, k, m, n, r, s, v, w, x, y_
  - Protein data: _?, \*, #, -, b, j, x, z_

* `l` Sequence length. If unknown, it is better to choose a value which is
surely greater than the real sequence length. If `l` is found to be
insufficient, the array size is dynamically increased (not recommended from a
performance point of view). Default value should be sufficient for most
datasets
* `verbose` If `true`, print status reports
* `verbosestep` Print a status report every `verbosestep` read sequences

## Details

When computing evolutionary distances, don't put the gap symbol `-` among the
missing values. Indeed, indels are an important piece of information for genetic
distances.

Standard conversion tables for FASTA data:

+----------------------------------------+
|         Conversion table (DNA)         |
+----------------------------------------+
|          Nucleotide |  Code |  Integer |
+---------------------+-------+----------+
|           Adenosine |   A   |     1    |
|            Cytosine |   C   |     2    |
|             Guanine |   G   |     3    |
|           Thymidine |   T   |     4    |
|              Uracil |   U   |     4    |
|     Purine (A or G) |   R   |     5    |
| Pyrimidine (C or T) |   Y   |    19    |
|                Keto |   K   |    13    |
|               Amino |   M   |    14    |
|  Strong Interaction |   S   |    17    |
|    Weak Interaction |   W   |    18    |
|               Not A |   B   |    23    |
|               Not C |   D   |     7    |
|               Not G |   H   |    10    |
|          Not T or U |   V   |    20    |
|                 Any |   N   |     6    |
|                 Gap |   -   |    26    |
|              Masked |   X   |    28    |
+----------------------------------------+

+------------------------------------------------+
|           Conversion table (PROTEIN)           |
+------------------------------------------------+
|                  Amino Acid |  Code |  Integer |
+-----------------------------+-------+----------+
|                     Alanine |   A   |     1    |
|                    Arginine |   R   |     5    |
|                  Asparagine |   N   |     6    |
|               Aspartic acid |   D   |     7    |
|                    Cysteine |   C   |     2    |
|                   Glutamine |   Q   |     8    |
|               Glutamic acid |   E   |     9    |
|                     Glycine |   G   |     3    |
|                   Histidine |   H   |    10    |
|                  Isoleucine |   I   |    11    |
|                     Leucine |   L   |    12    |
|                      Lysine |   K   |    13    |
|                  Methionine |   M   |    14    |
|               Phenylalanine |   F   |    15    |
|                     Proline |   P   |    16    |
|                 Pyrrolysine |   O   |    21    |
|              Selenocysteine |   U   |    22    |
|                      Serine |   S   |    17    |
|                   Threonine |   T   |     4    |
|                  Tryptophan |   W   |    18    |
|                    Tyrosine |   Y   |    19    |
|                      Valine |   V   |    20    |
| Asparagine or Aspartic acid |   B   |    23    |
|  Glutamine or Glutamic acid |   Z   |    24    |
|       Leucine or Isoleucine |   J   |    25    |
|                         Gap |   -   |    26    |
|            Translation stop |   *   |    27    |
|                         Any |   X   |    28    |
+------------------------------------------------+

## Value

A tuple containing the following variables:

* `data` Multiple Sequence Alignment (MSA) encoded as a UInt8 matrix
* `id` Units' ids
* `ref` Reference sequence, i.e. a vector of the same length of the original
sequences storing the values of homogeneous sites. SNPs are instead represented
by a value of 29
"""
function readfasta(infile::AbstractString,
                   dna::Bool,
                   miss::Vector{UInt8},
                   l::Int,
                   verbose::Bool,
                   verbosestep::Int)
  # we read the file twice
  # the first time we check the total number of sequences, if each sequence
  # has the same length and at what location the SNPs are
  # the second time we store the data, keeping only the SNPs

  # note: this function has been written with a huge dataset in mind

  # disable status reports if verbosestep is not positive
  verbose = verbose && (verbosestep > 0)

  #    ?    *    #    -    b    d    h    k    m    n    r    s    v    w    x
  # 0x3f 0x2a 0x23 0x2d 0x62 0x64 0x68 0x6b 0x6d 0x6e 0x72 0x73 0x76 0x77 0x78
  #    y    j    z
  # 0x79 0x6a 0x7a
  if length(miss) == 0
    if dna
      miss = [0x3f, 0x2a, 0x23, 0x62, 0x64, 0x68, 0x6b, 0x6d, 0x6e, 0x72, 0x73,
              0x76, 0x77, 0x78, 0x79]
    else
      miss = [0x3f, 0x2a, 0x23, 0x62, 0x6a, 0x78, 0x7a]
    end
  else
    if length(miss) == 1
      if miss[1] == 0x00
        # we can't use 0x00 as an index in enc[miss]. Use instead 0x01: it is
        # not an ASCII character allowed in a sequence in a proper FASTA file
        miss[1] = 0x01
      else
        if zero(UInt8) < miss[1] < UInt8(128)
          if UInt8(64) < miss[1] < UInt8(91)
            # convert to lowercase
            miss[1] += UInt8(32)
          end
        else
          throw(KDomainError(string("Value 'miss[1]' is not in the range ",
                                    "[1, ..., 127]:", Int(miss[1]), ".")))
        end
      end
    else
      for i in 1:length(miss)
        if zero(UInt8) < miss[i] < UInt8(128)
          if UInt8(64) < miss[i] < UInt8(91)
            # convert to lowercase
            miss[i] += UInt8(32)
          end
        else
          throw(KDomainError(string("Value 'miss[", i, "]' is not in the ",
                                    "range [1, ..., 127]:", Int(miss[i]), ".")))
        end
      end
    end
  end

  if l < 1
    throw(KDomainError("Argument 'l' is not positive."))
  end

  f = open(infile, "r")

  s = strip(readuntil(f, '>'))::ASCIIString
  if length(s) == 0
    close(f)
    throw(KFASTAError("No sequence has been found."))
  elseif length(s) > 1
    close(f)
    throw(KFASTAError("First non empty row is not a sequence id."))
  end

  # we now know that there is a first sequence to read... but is the ID empty?
  sid = strip(readuntil(f, '\n'))::ASCIIString

  if length(sid) == 0
    close(f)
    throw(KFASTAError("Missing sequence identifier. Sequence: 1."))
  end

  seqlen = 0
  seqref = zeros(UInt8, l)
  missseqref = falses(length(seqref))

  # support variables
  c = '\0'
  w = 0
  u = 0x00

  # start reading the first sequence
  keepgoing = true
  while keepgoing
    s = strip(readline(f))::ASCIIString

    if length(s) > 0
      if s[1] != '>'
        w = seqlen + length(s)

        # do we have enough space to store the first sequence?
        if w > length(seqref)
          tmp = zeros(UInt8, w + l)
          seqref = copy!(tmp, seqref)

          tmp = falses(length(seqref))
          missseqref = copy!(tmp, missseqref)
        end

        for c in lowercase(s)
          u = UInt8(c)

          if !((u == UInt8(32)) || (UInt8(9) <= u <= UInt8(13))) # skip blanks
            seqlen += 1

            if !in(u, miss)
              seqref[seqlen] = u
              missseqref[seqlen] = false
            else
              seqref[seqlen] = UInt8('?')
              missseqref[seqlen] = true
            end
          end
        end
      else
        keepgoing = false
      end
    elseif eof(f)
      keepgoing = false
    end
  end

  if seqlen == 0
    close(f)
    throw(KFASTAError("Missing sequence. Sequence: 1."))
  end

  # at least a sequence has been found
  n = 1

  seqref = seqref[1:seqlen]
  missseqref = missseqref[1:seqlen]

  if s[1] == '>'
    if length(s) > 1
      sid = lstrip(s, '>')::ASCIIString
    else
      # there is only the '>' character
      close(f)
      throw(KFASTAError("Missing sequence identifier. Sequence: 2."))
    end
  end

  curlen = 0

  snp = falses(seqlen)
  seq = zeros(UInt8, seqlen)
  missseq = falses(seqlen)

  for line in eachline(f)
    s = strip(line)::ASCIIString

    if length(s) > 0
      if s[1] != '>'
        for c in lowercase(s)
          u = UInt8(c)

          if !((u == UInt8(32)) || (UInt8(9) <= u <= UInt8(13))) # skip blanks
            curlen += 1

            if curlen > seqlen
              close(f)
              throw(KFASTAError(string("Different sequence length: sequence ",
                                       n + 1, " (", sid, ") ", "is longer ",
                                       "than expected.")))
            end

            if !in(u, miss)
              seq[curlen] = u
              missseq[curlen] = false
            else
              seq[curlen] = UInt8('?')
              missseq[curlen] = true
            end
          end
        end
      else
        # we just finished scanning the previous sequence
        if curlen != seqlen
          close(f)
          throw(KFASTAError(string("Different sequence length: sequence ",
                                   n + 1, " (", sid, ") ", "is shorter than ",
                                   "expected.")))
        end

        n += 1

        for b in 1:seqlen
          if !missseq[b]
            if missseqref[b]
              # this sequence has a non-missing value where it is missing in the
              # reference sequence
              seqref[b] = seq[b]
              missseqref[b] = false
            else
              # to be compared, both values must be non-missing
              snp[b] = snp[b] || (seq[b] != seqref[b])
            end
          end
        end

        if verbose && (n % verbosestep == 0)
          println(n, " sequences have been pre-processed.")
        end

        if length(s) > 1
          sid = lstrip(s, '>')::ASCIIString
          curlen = 0
        else
          # there is only the '>' character
          close(f)
          throw(KFASTAError(string("Missing identifier at sequence ", n + 1,
                                   ".")))
        end
      end
    end
  end

  # by construction, the last sequence has not been pre-processed
  if curlen != seqlen
    close(f)
    throw(KFASTAError(string("Different sequence length: sequence ", n + 1,
                             "(", sid, ") ", "is shorter than expected.")))
  end

  n += 1

  for b in 1:seqlen
    if !missseq[b]
      if missseqref[b]
        # this sequence has a non-missing value where it is missing in the
        # reference sequence
        seqref[b] = seq[b]
        missseqref[b] = false
      else
        # to be compared, both values must be non-missing
        snp[b] = snp[b] || (seq[b] != seqref[b])
      end
    end
  end

  m = sum(snp)

  if verbose
    println("Found ", n, " sequences: ", m, " SNPs out of ", seqlen,
            " total sites.\nProcessing data...")
  end

  # sequences are now guaranteed to be encoded in ASCII. We can safely use the
  # ASCII table to convert characters into integers. The null character '\0' is
  # not allowed and we won't have problems with BoundsErrors

  #   ?  a  c   g   t   r   n   d   q   e   h   i   l   k   m   f   p   s   w
  #  63 97 99 103 116 114 110 100 113 101 104 105 108 107 109 102 112 115 119
  #
  #   y   v   o   u  b   z   j  -  *   x
  # 121 118 111 117 98 122 106 45 42 120
  enc = zeros(UInt8, 127)
  enc[[ 63,  97,  99, 103, 116, 114, 110, 100, 113, 101, 104, 105, 108, 107,
       109, 102, 112, 115, 119, 121, 118, 111, 117,  98, 122, 106,  45,  42,
       120]] = round(UInt8, 0:28)
  enc[miss] = zero(UInt8)

  if dna
    enc[117] = UInt8(4) # 'u'
  end

  data = zeros(UInt8, m, n)
  id = Array(ASCIIString, n)

  # go back at the beginning of the file and start again
  seekstart(f)

  # we already checked that the first non-empty element is the id
  i = 1

  s = strip(readuntil(f, '>'))::ASCIIString
  id[i] = strip(readuntil(f, '\n'))::ASCIIString

  i1 = 0
  i2 = 0
  for line in eachline(f)
    s = strip(line)::ASCIIString

    if length(s) > 0
      if s[1] != '>'
        for c in lowercase(s)
          u = UInt8(c)

          if !((u == UInt8(32)) || (UInt8(9) <= u <= UInt8(13))) # skip blanks
            i1 += 1

            if snp[i1]
              i2 += 1
              data[i2, i] = enc[u]
            end
          end
        end
      else
        i += 1

        if verbose && (i % verbosestep == 0)
          println(i, " sequences have been processed.")
        end

        # move on with the next sequence
        id[i] = lstrip(s, '>')::ASCIIString

        i1 = 0
        i2 = 0
      end
    end
  end

  if verbose
    println("All ", i, " sequences have been processed.")
  end

  close(f)

  ref = fill(UInt8(29), seqlen)
  ref[!snp] = enc[seqref[!snp]]

  (data, id, ref)
end
