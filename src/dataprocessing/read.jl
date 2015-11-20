# This file is part of K-Pax3. License is MIT.

"""
# Read a FASTA formatted file

## Description

Read data in FASTA format and convert it to an integer matrix. Sequences are
required to be aligned. Only polymorphic columns are stored.

## Usage

readfasta(infile; dna=true, miss=zeros(UInt8, 0), l=100000000, t=500)

## Arguments

* `infile` Path to the input data file
* `dna` `true` if reading DNA data or `false` if reading protein data
* `miss` Values to be considered missing. Default values are (UInt8 version):
 - DNA data

        miss = ['?', '\*', #', '-', 'b', 'd', 'h', 'k', 'm', 'n', 'r', 's', 'v',
                'w', 'x', 'y']

 - Protein data

        miss = ['?', '\*', '#', '-', 'b', 'j', 'x', 'z']

* `l` Sequence length. If unknown, it is better to choose a value which is
surely greater than the real sequence length. If `l` is found to be
insufficient, the array size is dynamically increased (not recommended from a
performance point of view). Default value should be sufficient for most
datasets
* `t` If positive, print a status report every `t` read sequences

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
function readfasta(infile::AbstractString;
                   dna::Bool=true,
                   miss::Array{UInt8, 1}=zeros(UInt8, 0),
                   l::Int=100000000,
                   t::Int=500)
  # we read the file twice
  # the first time we check the total number of sequences, if each sequence
  # has the same length and at what location the SNPs are
  # the second time we store the data, keeping only the SNPs

  # note: this function has been written with a huge dataset in mind

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
    # convert to lowercase
    for i in 1:length(miss)
      miss[i] = (0x40 < miss[i] < 0x5b) ? (miss[i] + 0x20) : miss[i]
    end
  end

  f = open(infile, "r")

  s = strip(readuntil(f, '>'))::ASCIIString
  if length(s) == 0
    close(f)
    throw(Kpax3FASTAError("No sequence has been found."))
  elseif length(s) > 1
    close(f)
    throw(Kpax3FASTAError("First non empty row is not a sequence id."))
  end

  # we now know that there is a first sequence to read... but is the ID empty?
  sid = strip(readuntil(f, '\n'))::ASCIIString

  if length(sid) == 0
    close(f)
    throw(Kpax3FASTAError("Missing sequence identifier. Sequence: 1."))
  end

  seqlen = zero(Int)
  seqref = zeros(UInt8, l)
  missseqref = falses(length(seqref))

  # support variables
  c = '\0'
  L = zero(Int)
  w = zero(Int)
  u = zero(UInt8)

  # start reading the first sequence
  keepgoing = true
  while keepgoing
    s = strip(readline(f))::ASCIIString

    if length(s) > 0
      if s[1] != '>'
        w = seqlen + length(s)

        # do we have enough space to store the first sequence?
        if w > length(seqref)
          L = length(seqref) + l * fld(w, length(seqref))
          tmp = zeros(UInt8, L)
          tmp[1:length(seqref)] = seqref
          seqref = tmp

          tmp = falses(length(seqref))
          tmp[1:length(missseqref)] = missseqref
          missseqref = tmp
        end

        for c in lowercase(s)
          u = UInt8(c)

          if !((u == 0x20) || (0x09 <= u <= 0x0d)) # skip blanks
            seqlen += 1

            if !in(u, miss)
              seqref[seqlen] = u
              missseqref[seqlen] = false
            else
              seqref[seqlen] = 0x3f # '?'
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
    throw(Kpax3FASTAError("Missing sequence. Sequence: 1."))
  end

  # at least a sequence has been found
  n = one(Int)

  seqref = seqref[1:seqlen]
  missseqref = missseqref[1:seqlen]

  if s[1] == '>'
    if length(s) > 1
      sid = lstrip(s, '>')::ASCIIString
    else
      # there is only the '>' character
      close(f)
      throw(Kpax3FASTAError("Missing sequence identifier. Sequence: 2."))
    end
  end

  curlen = zero(Int)

  snp = falses(seqlen)
  seq = zeros(UInt8, seqlen)
  missseq = falses(seqlen)
  j = falses(seqlen)

  for line in eachline(f)
    s = strip(line)::ASCIIString

    if length(s) > 0
      if s[1] != '>'
        for c in lowercase(s)
          u = UInt8(c)

          if !((u == 0x20) || (0x09 <= u <= 0x0d)) # skip blanks
            curlen += 1

            if curlen > seqlen
              close(f)
              throw(Kpax3FASTAError(string("Different sequence length: ",
                                           "sequence ", n + 1, " (", sid, ") ",
                                           "is longer than expected.")))
            end

            if !in(u, miss)
              seq[curlen] = u
              missseq[curlen] = false
            else
              seq[curlen] = 0x3f # '?'
              missseq[curlen] = true
            end
          end
        end
      else
        # we just finished scanning the previous sequence
        if curlen != seqlen
          close(f)
          throw(Kpax3FASTAError(string("Different sequence length: ",
                                       "sequence ", n + 1, " (", sid, ") ",
                                       "is shorter than expected.")))
        end

        n += 1

        # maybe this sequence has a non-missing value where it is missing
        # in the reference sequence
        j[:] = missseqref & (!missseq)
        seqref[j] = seq[j]
        missseqref[j] = false

        # to be compared, both values must be non-missing
        j[:] = !(missseqref | missseq)
        snp[j] = snp[j] | (seq[j] .!= seqref[j])

        if (t > 0) && (n % t == 0)
          println(n, " sequences have been pre-processed.")
        end

        if length(s) > 1
          sid = lstrip(s, '>')::ASCIIString
          curlen = 0
        else
          # there is only the '>' character
          close(f)
          throw(Kpax3FASTAError(string("Missing identifier at sequence ",
                                       n + 1, ".")))
        end
      end
    end
  end

  # by construction, the last sequence has not been pre-processed
  if curlen != seqlen
    close(f)
    throw(Kpax3FASTAError(string("Different sequence length: sequence ", n + 1,
                                 "(", sid, ") ", "is shorter than expected.")))
  end

  n += 1

  j[:] = missseqref & (!missseq)
  seqref[j] = seq[j]
  missseqref[j] = false

  j[:] = !(missseqref | missseq)
  snp[j] = snp[j] | (seq[j] .!= seqref[j])

  m = sum(snp)::Int

  if t > 0
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
       120]] = [0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09,
                0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13,
                0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c] # 0:28
  enc[miss] = 0x00

  if dna
    enc[117] = 0x04 # 'u'
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

          if !((u == 0x20) || (0x09 <= u <= 0x0d)) # skip blanks
            i1 += 1

            if snp[i1]
              i2 += 1
              data[i2, i] = enc[u]
            end
          end
        end
      else
        i += 1

        if (t > 0) && (i % t == 0)
          println(i, " sequences have been processed.")
        end

        # move on with the next sequence
        id[i] = lstrip(s, '>')::ASCIIString

        i1 = 0
        i2 = 0
      end
    end
  end

  if t > 0
    println("All ", i, " sequences have been processed.")
  end

  close(f)

  ref = zeros(UInt8, seqlen)
  ref[:] = 29
  ref[!snp] = enc[seqref[!snp]]

  (data, id, ref)
end
