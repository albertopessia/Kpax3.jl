# This file is part of Kpax3. License is MIT.

f = "../build/test.bin"

##############
# Exceptions #
##############

@test_throws KFASTAError NucleotideData(KSettings("data/empty_file.fasta", f))
@test_throws KFASTAError NucleotideData(KSettings("data/no_1st_seq.fasta", f))
@test_throws KFASTAError NucleotideData(KSettings("data/no_id_char.fasta", f))
@test_throws KFASTAError NucleotideData(KSettings("data/no_nth_seq.fasta", f))

@test_throws TypeError NucleotideData(KSettings("data/utf8_id.fasta", f))
@test_throws TypeError NucleotideData(KSettings("data/utf8_seq.fasta", f))

@test_throws KFASTAError AminoAcidData(KSettings("data/empty_file.fasta", f))
@test_throws KFASTAError AminoAcidData(KSettings("data/no_1st_seq.fasta", f))
@test_throws KFASTAError AminoAcidData(KSettings("data/no_id_char.fasta", f))
@test_throws KFASTAError AminoAcidData(KSettings("data/no_nth_seq.fasta", f))

@test_throws TypeError AminoAcidData(KSettings("data/utf8_id.fasta", f))
@test_throws TypeError AminoAcidData(KSettings("data/utf8_seq.fasta", f))

#################################
# FASTA file filled with blanks #
#################################

nt = NucleotideData(KSettings("data/blanks.fasta", f))
@test nt.data == UInt8[0 1;
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
@test nt.id == ASCIIString["ID1", "ID5"]
@test nt.ref == UInt8['a', 't', 'g', '.', '.', '.', 'g', '.', '.', 'a', '.',
                      'a']
@test nt.val == UInt8['a', 'g', 'a', 'g', 'c', 'g', 'c', 'g', 'a', 'g', 'c',
                      'g']
@test nt.key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]

nt = NucleotideData(KSettings("data/blanks.fasta", f, l=1))
@test nt.data == UInt8[0 1;
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
@test nt.id == ASCIIString["ID1", "ID5"]
@test nt.ref == UInt8['a', 't', 'g', '.', '.', '.', 'g', '.', '.', 'a', '.',
                      'a']
@test nt.val == UInt8['a', 'g', 'a', 'g', 'c', 'g', 'c', 'g', 'a', 'g', 'c',
                      'g']
@test nt.key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]

aa = AminoAcidData(KSettings("data/blanks.fasta", f))
@test aa.data == UInt8[0 1;
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
@test aa.id == ASCIIString["ID1", "ID5"]
@test aa.ref == UInt8['a', 't', 'g', '.', '.', '.', 'g', '.', '.', 'a', '.',
                      'a']
@test aa.val == UInt8['a', 'g', 'a', 'g', 'c', 'g', 'c', 'g', 'a', 'g', 'c',
                      'g']
@test aa.key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]

aa = AminoAcidData(KSettings("data/blanks.fasta", f, l=1))
@test aa.data == UInt8[0 1;
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
@test aa.id == ASCIIString["ID1", "ID5"]
@test aa.ref == UInt8['a', 't', 'g', '.', '.', '.', 'g', '.', '.', 'a', '.',
                      'a']
@test aa.val == UInt8['a', 'g', 'a', 'g', 'c', 'g', 'c', 'g', 'a', 'g', 'c',
                      'g']
@test aa.key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]

#########################
# Proper DNA FASTA file #
#########################

nt = NucleotideData(KSettings("data/proper_nt.fasta", f))
@test nt.data == UInt8[0 0 0 0 0 1;
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
@test nt.id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test nt.ref == UInt8['a', '.', 'g', '.', '.', '.', 'g', '.', '.', '.', '.',
                      'a']
@test nt.val == UInt8['c', 't', 'a', 'g', 'a', 'g', 'c', 'g', 'a', 'c', 'g',
                      'a', 'g', 't', 'a', 't', 'c', 'g']
@test nt.key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8]

nt = NucleotideData(KSettings("data/proper_nt.fasta", f, l=1))
@test nt.data == UInt8[0 0 0 0 0 1;
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
@test nt.id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test nt.ref == UInt8['a', '.', 'g', '.', '.', '.', 'g', '.', '.', '.', '.',
                      'a']
@test nt.val == UInt8['c', 't', 'a', 'g', 'a', 'g', 'c', 'g', 'a', 'c', 'g',
                      'a', 'g', 't', 'a', 't', 'c', 'g']
@test nt.key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8]

# consider all characters
nt = NucleotideData(KSettings("data/proper_nt.fasta", f, miss=zeros(UInt8, 1)))
@test nt.data == UInt8[0 0 0 0 0 1;
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
@test nt.id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test nt.ref == UInt8['a', '.', 'g', '.', '.', '.', 'g', '.', '.', '.', '.',
                      '.']
@test nt.val == UInt8['c', 't', 'a', 'g', 'a', 'g', 'x', 'c', 'g', 'a', 'c',
                      'g', 'a', 'g', 't', 'a', 't', 'c', 'g', 'a', 'x']
@test nt.key == [1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9]

#############################
# Proper Protein FASTA file #
#############################

aa = AminoAcidData(KSettings("data/proper_aa.fasta", f))
@test aa.data == UInt8[0 0 0 0 0 1;
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
@test aa.id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test aa.ref == UInt8['m', '.', 'a', '.', '.', '.', 'v', '.', '.', '.', '.',
                      'f']
@test aa.val == UInt8['e', 'k', 'i', 'k', 'l', 'v', 'l', 'v', 'c', 'l', 't',
                      'e', 'l', 'm', 'c', 'y', 'a', 't']
@test aa.key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8]

aa = AminoAcidData(KSettings("data/proper_aa.fasta", f, l=1))
@test aa.data == UInt8[0 0 0 0 0 1;
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
@test aa.id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test aa.ref == UInt8['m', '.', 'a', '.', '.', '.', 'v', '.', '.', '.', '.',
                      'f']
@test aa.val == UInt8['e', 'k', 'i', 'k', 'l', 'v', 'l', 'v', 'c', 'l', 't',
                      'e', 'l', 'm', 'c', 'y', 'a', 't']
@test aa.key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8]

# consider all characters
aa = AminoAcidData(KSettings("data/proper_aa.fasta", f, miss=zeros(UInt8, 1)))
@test aa.data == UInt8[0 0 0 0 0 1;
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
@test aa.id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test aa.ref == UInt8['m', '.', 'a', '.', '.', '.', 'v', '.', '.', '.', '.',
                      '.']
@test aa.val == UInt8['e', 'k', 'i', 'k', 'l', 'v', 'x', 'l', 'v', 'c', 'l',
                      't', 'e', 'l', 'm', 'c', 'y', 'a', 't', 'f', 'x']
@test aa.key == [1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9]

####################
# Input and Output #
####################
nt = NucleotideData(KSettings("data/proper_nt.fasta", f))

# TODO: Test exception when saving to a location without writing permissions
save("../build/nt.jld", nt)
@test isfile("../build/nt.jld")

@test_throws SystemError loadnt("../build/non_existent.file")

nt = loadnt("../build/nt.jld")

@test isa(nt, NucleotideData)
@test nt.data == UInt8[0 0 0 0 0 1;
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
@test nt.id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test nt.ref == UInt8['a', '.', 'g', '.', '.', '.', 'g', '.', '.', '.', '.',
                      'a']
@test nt.val == UInt8['c', 't', 'a', 'g', 'a', 'g', 'c', 'g', 'a', 'c', 'g',
                      'a', 'g', 't', 'a', 't', 'c', 'g']
@test nt.key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8]

aa = AminoAcidData(KSettings("data/proper_aa.fasta", f))

# TODO: Test exception when saving to a location without writing permissions
save("../build/aa.jld", aa)
@test isfile("../build/aa.jld")

@test_throws SystemError loadaa("../build/non_existent.file")

aa = loadaa("../build/aa.jld")
@test isa(aa, AminoAcidData)
@test aa.data == UInt8[0 0 0 0 0 1;
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
@test aa.id == ASCIIString["ID1", "ID2", "ID3", "ID4", "ID5", "ID6"]
@test aa.ref == UInt8['m', '.', 'a', '.', '.', '.', 'v', '.', '.', '.', '.',
                      'f']
@test aa.val == UInt8['e', 'k', 'i', 'k', 'l', 'v', 'l', 'v', 'c', 'l', 't',
                      'e', 'l', 'm', 'c', 'y', 'a', 't']
@test aa.key == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 8, 8]