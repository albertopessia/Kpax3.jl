# This file is part of Kpax3. License is MIT.

##############
# Exceptions #
##############
@test_throws KDomainError readfasta("data/proper_nt.fasta", true, UInt8[0, 63],
                                    10, false, 0)
@test_throws KDomainError readfasta("data/proper_nt.fasta", true,
                                    zeros(UInt8, 0), -1, false, 0)

@test_throws KFASTAError readfasta("data/empty_file.fasta", true,
                                   zeros(UInt8, 0), 100000000, false, 0)
@test_throws KFASTAError readfasta("data/no_1st_seq.fasta", true,
                                   zeros(UInt8, 0), 100000000, false, 0)
@test_throws KFASTAError readfasta("data/no_id_char.fasta", true,
                                   zeros(UInt8, 0), 100000000, false, 0)
@test_throws KFASTAError readfasta("data/no_nth_seq.fasta", true,
                                   zeros(UInt8, 0), 100000000, false, 0)

@test_throws TypeError readfasta("data/utf8_id.fasta", true, zeros(UInt8, 0),
                                 100000000, false, 0)
@test_throws TypeError readfasta("data/utf8_seq.fasta", true, zeros(UInt8, 0),
                                 100000000, false, 0)

#################################
# FASTA file filled with blanks #
#################################
data, id, refseq = readfasta("data/blanks.fasta", true, zeros(UInt8, 0),
                             100000000, false, 0)
@test data == UInt8[3 1;
                    1 3;
                    3 2;
                    2 3;
                    1 3;
                    3 2]
@test id == ASCIIString["ID1", "ID5"]
@test refseq == UInt8[1, 4, 3, 29, 29, 29, 3, 29, 29, 1, 29, 1]

@test_throws KDomainError categorical2binary(data, 0x01)

bindata, val, key = categorical2binary(data, 0x1c)
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
@test val == UInt8[1, 3, 1, 3, 2, 3, 2, 3, 1, 3, 2, 3]
@test key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]

data, id, refseq = readfasta("data/blanks.fasta", true, zeros(UInt8, 0), 1,
                             false, 0)
@test data == UInt8[3 1;
                    1 3;
                    3 2;
                    2 3;
                    1 3;
                    3 2]
@test id == ASCIIString["ID1", "ID5"]
@test refseq == UInt8[1, 4, 3, 29, 29, 29, 3, 29, 29, 1, 29, 1]

#########################
# Proper DNA FASTA file #
#########################
data, id, refseq = readfasta("data/proper_nt.fasta", true, zeros(UInt8, 0),
                             100000000, false, 0)
@test data == UInt8[4 4 4 4 4 2;
                    3 3 1 3 1 1;
                    1 1 0 3 3 0;
                    3 3 3 2 2 2;
                    2 3 1 2 3 2;
                    1 1 1 3 3 4;
                    1 1 4 4 1 1;
                    3 2 3 3 2 3]
@test id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test refseq == UInt8[1, 29, 3, 29, 29, 29, 3, 29, 29, 29, 29, 1]

@test_throws KDomainError categorical2binary(data, 0x01)

bindata, val, key = categorical2binary(data, 0x1c)
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
@test val == UInt8[2, 4, 1, 3, 1, 3, 2, 3, 1, 2, 3, 1, 3, 4, 1, 4, 2, 3]
@test key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8]

data, id, refseq = readfasta("data/proper_nt.fasta", true, zeros(UInt8, 0), 1,
                             false, 0)
@test data == UInt8[4 4 4 4 4 2;
                    3 3 1 3 1 1;
                    1 1 0 3 3 0;
                    3 3 3 2 2 2;
                    2 3 1 2 3 2;
                    1 1 1 3 3 4;
                    1 1 4 4 1 1;
                    3 2 3 3 2 3]
@test id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test refseq == UInt8[1, 29, 3, 29, 29, 29, 3, 29, 29, 29, 29, 1]

# consider all characters
data, id, refseq = readfasta("data/proper_nt.fasta", true, zeros(UInt8, 1),
                             100000000, false, 0)
@test data == UInt8[ 4  4  4  4  4  2;
                     3  3  1  3  1  1;
                     1  1 28  3  3 28;
                     3  3  3  2  2  2;
                     2  3  1  2  3  2;
                     1  1  1  3  3  4;
                     1  1  4  4  1  1;
                     3  2  3  3  2  3;
                     1 28  1  1  1  1]
@test id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test refseq == UInt8[1, 29, 3, 29, 29, 29, 3, 29, 29, 29, 29, 29]

@test_throws KDomainError categorical2binary(data, 0x01)

bindata, val, key = categorical2binary(data, 0x1c)
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
@test val == UInt8[2, 4, 1, 3, 1, 3, 28, 2, 3, 1, 2, 3, 1, 3, 4, 1, 4, 2, 3, 1,
                   28]
@test key == [1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9]

#############################
# Proper Protein FASTA file #
#############################
data, id, refseq = readfasta("data/proper_aa.fasta", false, zeros(UInt8, 0),
                             100000000, false, 0)
@test data == UInt8[13 13 13 13 13  9;
                    13 13 11 13 11 11;
                    12 12  0 20 20  0;
                    12 12 12 20 20 20;
                     4 12  2  4 12  4;
                    12 12 12 14 14  9;
                     2  2 19 19  2  2;
                     1  1  4  1  4  4]
@test id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test refseq == UInt8[14, 29, 1, 29, 29, 29, 20, 29, 29, 29, 29, 15]

bindata, val, key = categorical2binary(data, 0x1c)
@test bindata == UInt8[0 0 0 0 0 1;
                       1 1 1 1 1 0;
                       0 0 1 0 1 1;
                       1 1 0 1 0 0;
                       1 1 0 0 0 0;
                       0 0 0 1 1 0;
                       1 1 1 0 0 0;
                       0 0 0 1 1 1;
                       0 0 1 0 0 0;
                       1 0 0 1 0 1;
                       0 1 0 0 1 0;
                       0 0 0 0 0 1;
                       1 1 1 0 0 0;
                       0 0 0 1 1 0;
                       1 1 0 0 1 1;
                       0 0 1 1 0 0;
                       1 1 0 1 0 0;
                       0 0 1 0 1 1]
@test val == UInt8[9, 13, 11, 13, 12, 20, 12, 20, 2, 4, 12, 9, 12, 14, 2, 19, 1,
                   4]
@test key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8]

data, id, refseq = readfasta("data/proper_aa.fasta", false, zeros(UInt8, 0), 1,
                             false, 0)
@test data == UInt8[13 13 13 13 13  9;
                    13 13 11 13 11 11;
                    12 12  0 20 20  0;
                    12 12 12 20 20 20;
                     4 12  2  4 12  4;
                    12 12 12 14 14  9;
                     2  2 19 19  2  2;
                     1  1  4  1  4  4]
@test id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test refseq == UInt8[14, 29, 1, 29, 29, 29, 20, 29, 29, 29, 29, 15]

# consider all characters
data, id, refseq = readfasta("data/proper_aa.fasta", true, zeros(UInt8, 1),
                             100000000, false, 0)
@test data == UInt8[13 13 13 13 13  9;
                    13 13 11 13 11 11;
                    12 12 28 20 20 28;
                    12 12 12 20 20 20;
                     4 12  2  4 12  4;
                    12 12 12 14 14  9;
                     2  2 19 19  2  2;
                     1  1  4  1  4  4;
                    15 28 15 15 15 15]
@test id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test refseq == UInt8[14, 29, 1, 29, 29, 29, 20, 29, 29, 29, 29, 29]

bindata, val, key = categorical2binary(data, 0x1c)
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
                       1 0 0 1 0 1;
                       0 1 0 0 1 0;
                       0 0 0 0 0 1;
                       1 1 1 0 0 0;
                       0 0 0 1 1 0;
                       1 1 0 0 1 1;
                       0 0 1 1 0 0;
                       1 1 0 1 0 0;
                       0 0 1 0 1 1;
                       1 0 1 1 1 1;
                       0 1 0 0 0 0]
@test val == UInt8[9, 13, 11, 13, 12, 20, 28, 12, 20, 2, 4, 12, 9, 12, 14, 2,
                   19, 1, 4, 15, 28]
@test key == [1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9]
