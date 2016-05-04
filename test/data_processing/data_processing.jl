# This file is part of Kpax3. License is MIT.

missdna = UInt8['?', '*', '#', '-', 'b', 'd', 'h', 'k', 'm', 'n', 'r', 's', 'v',
                'w', 'x', 'y', 'j', 'z']
misspro = UInt8['?', '*', '#', '-', 'b', 'j', 'x', 'z']

##############
# Exceptions #
##############
@test_throws KFASTAError readfasta("data/empty_file.fasta", false, missdna,
                                   100000000, false, 0)
@test_throws KFASTAError readfasta("data/no_1st_seq.fasta", false, missdna,
                                   100000000, false, 0)
@test_throws KFASTAError readfasta("data/no_id_char.fasta", false, missdna,
                                   100000000, false, 0)
@test_throws KFASTAError readfasta("data/no_nth_seq.fasta", false, missdna,
                                   100000000, false, 0)

@test_throws TypeError readfasta("data/utf8_id.fasta", false, missdna,
                                 100000000, false, 0)
@test_throws TypeError readfasta("data/utf8_seq.fasta", false, missdna,
                                 100000000, false, 0)

#################################
# FASTA file filled with blanks #
#################################
(data, id, refseq) = readfasta("data/blanks.fasta", false, missdna, 100000000,
                               false, 0)
@test data == UInt8['g' 'a';
                    'a' 'g';
                    'g' 'c';
                    'c' 'g';
                    'a' 'g';
                    'g' 'c']
@test id == ASCIIString["ID1", "ID5"]
@test refseq == UInt8['a', 't', 'g', '.', '.', '.', 'g', '.', '.', 'a', '.',
                      'a']

@test_throws KDomainError categorical2binary(data, UInt8(1), UInt8('?'))

(bindata, val, key) = categorical2binary(data, UInt8(127), UInt8('?'))
@test bindata == UInt8[0 1;
                       1 0;
                       1 0;
                       0 1;
                       0 1;
                       1 0;
                       1 0;
                       0 1;
                       1 0;
                       0 1;
                       0 1;
                       1 0]
@test val == UInt8['a', 'g', 'a', 'g', 'c', 'g', 'c', 'g', 'a', 'g', 'c', 'g']
@test key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]

(data, id, refseq) = readfasta("data/blanks.fasta", false, missdna, 1, false, 0)
@test data == UInt8['g' 'a';
                    'a' 'g';
                    'g' 'c';
                    'c' 'g';
                    'a' 'g';
                    'g' 'c']
@test id == ASCIIString["ID1", "ID5"]
@test refseq == UInt8['a', 't', 'g', '.', '.', '.', 'g', '.', '.', 'a', '.',
                      'a']

#########################
# Proper DNA FASTA file #
#########################
(data, id, refseq) = readfasta("data/proper_nt.fasta", false, missdna,
                               100000000, false, 0)

@test data == UInt8['t' 't' 't' 't' 't' 'c';
                    'g' 'g' 'a' 'g' 'a' 'a';
                    'a' 'a' '?' 'g' 'g' '?';
                    'g' 'g' 'g' 'c' 'c' 'c';
                    'c' 'g' 'a' 'c' 'g' 'c';
                    'a' 'a' 'a' 'g' 'g' 't';
                    'a' 'a' 't' 't' 'a' 'a';
                    'g' 'c' 'g' 'g' 'c' 'g']
@test id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test refseq == UInt8['a', '.', 'g', '.', '.', '.', 'g', '.', '.', '.', '.',
                      'a']

@test_throws KDomainError categorical2binary(data, UInt8(1), UInt8('?'))

(bindata, val, key) = categorical2binary(data, UInt8(127), UInt8('?'))
@test bindata == UInt8[0 0 0 0 0 1;
                       1 1 1 1 1 0;
                       0 0 1 0 1 1;
                       1 1 0 1 0 0;
                       1 1 0 0 0 0;
                       0 0 0 1 1 0;
                       0 0 0 1 1 1;
                       1 1 1 0 0 0;
                       0 0 1 0 0 0;
                       1 0 0 1 0 1;
                       0 1 0 0 1 0;
                       1 1 1 0 0 0;
                       0 0 0 1 1 0;
                       0 0 0 0 0 1;
                       1 1 0 0 1 1;
                       0 0 1 1 0 0;
                       0 1 0 0 1 0;
                       1 0 1 1 0 1]
@test val == UInt8['c', 't', 'a', 'g', 'a', 'g', 'c', 'g', 'a', 'c', 'g', 'a',
                   'g', 't', 'a', 't', 'c', 'g']
@test key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8]

(data, id, refseq) = readfasta("data/proper_nt.fasta", false, missdna, 1, false,
                               0)
@test data == UInt8['t' 't' 't' 't' 't' 'c';
                    'g' 'g' 'a' 'g' 'a' 'a';
                    'a' 'a' '?' 'g' 'g' '?';
                    'g' 'g' 'g' 'c' 'c' 'c';
                    'c' 'g' 'a' 'c' 'g' 'c';
                    'a' 'a' 'a' 'g' 'g' 't';
                    'a' 'a' 't' 't' 'a' 'a';
                    'g' 'c' 'g' 'g' 'c' 'g']
@test id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test refseq == UInt8['a', '.', 'g', '.', '.', '.', 'g', '.', '.', '.', '.',
                      'a']

# consider all characters
(data, id, refseq) = readfasta("data/proper_nt.fasta", false, zeros(UInt8, 1),
                               100000000, false, 0)
@test data == UInt8['t' 't' 't' 't' 't' 'c';
                    'g' 'g' 'a' 'g' 'a' 'a';
                    'a' 'a' 'x' 'g' 'g' 'x';
                    'g' 'g' 'g' 'c' 'c' 'c';
                    'c' 'g' 'a' 'c' 'g' 'c';
                    'a' 'a' 'a' 'g' 'g' 't';
                    'a' 'a' 't' 't' 'a' 'a';
                    'g' 'c' 'g' 'g' 'c' 'g';
                    'a' 'x' 'a' 'a' 'a' 'a']
@test id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test refseq == UInt8['a', '.', 'g', '.', '.', '.', 'g', '.', '.', '.', '.',
                      '.']

@test_throws KDomainError categorical2binary(data, UInt8(1), UInt8('\0'))

(bindata, val, key) = categorical2binary(data, UInt8(127), UInt8('\0'))
@test bindata == UInt8[0 0 0 0 0 1;
                       1 1 1 1 1 0;
                       0 0 1 0 1 1;
                       1 1 0 1 0 0;
                       1 1 0 0 0 0;
                       0 0 0 1 1 0;
                       0 0 1 0 0 1;
                       0 0 0 1 1 1;
                       1 1 1 0 0 0;
                       0 0 1 0 0 0;
                       1 0 0 1 0 1;
                       0 1 0 0 1 0;
                       1 1 1 0 0 0;
                       0 0 0 1 1 0;
                       0 0 0 0 0 1;
                       1 1 0 0 1 1;
                       0 0 1 1 0 0;
                       0 1 0 0 1 0;
                       1 0 1 1 0 1;
                       1 0 1 1 1 1;
                       0 1 0 0 0 0]
@test val == UInt8['c', 't', 'a', 'g', 'a', 'g', 'x', 'c', 'g', 'a', 'c', 'g',
                   'a', 'g', 't', 'a', 't', 'c', 'g', 'a', 'x']
@test key == [1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9]

#############################
# Proper Protein FASTA file #
#############################
(data, id, refseq) = readfasta("data/proper_aa.fasta", true, misspro, 100000000,
                               false, 0)

@test data == UInt8['k' 'k' 'k' 'k' 'k' 'e';
                    'k' 'k' 'i' 'k' 'i' 'i';
                    'l' 'l' '?' 'v' 'v' '?';
                    'l' 'l' 'l' 'v' 'v' 'v';
                    't' 'l' 'c' 't' 'l' 't';
                    'l' 'l' 'l' 'm' 'm' 'e';
                    'c' 'c' 'y' 'y' 'c' 'c';
                    'a' 'a' 't' 'a' 't' 't']
@test id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test refseq == UInt8['m', '.', 'a', '.', '.', '.', 'v', '.', '.', '.', '.',
                      'f']

(bindata, val, key) = categorical2binary(data, UInt8(127), UInt8('?'))
@test bindata == UInt8[0 0 0 0 0 1;
                       1 1 1 1 1 0;
                       0 0 1 0 1 1;
                       1 1 0 1 0 0;
                       1 1 0 0 0 0;
                       0 0 0 1 1 0;
                       1 1 1 0 0 0;
                       0 0 0 1 1 1;
                       0 0 1 0 0 0;
                       0 1 0 0 1 0;
                       1 0 0 1 0 1;
                       0 0 0 0 0 1;
                       1 1 1 0 0 0;
                       0 0 0 1 1 0;
                       1 1 0 0 1 1;
                       0 0 1 1 0 0;
                       1 1 0 1 0 0;
                       0 0 1 0 1 1]
@test val == UInt8['e', 'k', 'i', 'k', 'l', 'v', 'l', 'v', 'c', 'l', 't', 'e',
                   'l', 'm', 'c', 'y', 'a', 't']
@test key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8]

(data, id, refseq) = readfasta("data/proper_aa.fasta", true, misspro, 1,
                             false, 0)
@test data == UInt8['k' 'k' 'k' 'k' 'k' 'e';
                    'k' 'k' 'i' 'k' 'i' 'i';
                    'l' 'l' '?' 'v' 'v' '?';
                    'l' 'l' 'l' 'v' 'v' 'v';
                    't' 'l' 'c' 't' 'l' 't';
                    'l' 'l' 'l' 'm' 'm' 'e';
                    'c' 'c' 'y' 'y' 'c' 'c';
                    'a' 'a' 't' 'a' 't' 't']
@test id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test refseq == UInt8['m', '.', 'a', '.', '.', '.', 'v', '.', '.', '.', '.',
                      'f']

# consider all characters
(data, id, refseq) = readfasta("data/proper_aa.fasta", true, zeros(UInt8, 1),
                               100000000, false, 0)
@test data == UInt8['k' 'k' 'k' 'k' 'k' 'e';
                    'k' 'k' 'i' 'k' 'i' 'i';
                    'l' 'l' 'x' 'v' 'v' 'x';
                    'l' 'l' 'l' 'v' 'v' 'v';
                    't' 'l' 'c' 't' 'l' 't';
                    'l' 'l' 'l' 'm' 'm' 'e';
                    'c' 'c' 'y' 'y' 'c' 'c';
                    'a' 'a' 't' 'a' 't' 't';
                    'f' 'x' 'f' 'f' 'f' 'f']
@test id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test refseq == UInt8['m', '.', 'a', '.', '.', '.', 'v', '.', '.', '.', '.',
                      '.']

(bindata, val, key) = categorical2binary(data, UInt8(127), UInt8('\0'))
@test bindata == UInt8[0 0 0 0 0 1;
                       1 1 1 1 1 0;
                       0 0 1 0 1 1;
                       1 1 0 1 0 0;
                       1 1 0 0 0 0;
                       0 0 0 1 1 0;
                       0 0 1 0 0 1;
                       1 1 1 0 0 0;
                       0 0 0 1 1 1;
                       0 0 1 0 0 0;
                       0 1 0 0 1 0;
                       1 0 0 1 0 1;
                       0 0 0 0 0 1;
                       1 1 1 0 0 0;
                       0 0 0 1 1 0;
                       1 1 0 0 1 1;
                       0 0 1 1 0 0;
                       1 1 0 1 0 0;
                       0 0 1 0 1 1;
                       1 0 1 1 1 1;
                       0 1 0 0 0 0]
@test val == UInt8['e', 'k', 'i', 'k', 'l', 'v', 'x', 'l', 'v', 'c', 'l', 't',
                   'e', 'l', 'm', 'c', 'y', 'a', 't', 'f', 'x']
@test key == [1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9]
