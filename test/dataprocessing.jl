# This file is part of K-Pax3. License is MIT.

@test_throws Kpax3FASTAError readfasta("data/empty_file.fasta")
@test_throws Kpax3FASTAError readfasta("data/no_1st_seq.fasta")
@test_throws Kpax3FASTAError readfasta("data/no_id_char.fasta")
@test_throws Kpax3FASTAError readfasta("data/no_nth_seq.fasta")
@test_throws TypeError readfasta("data/utf8_id.fasta")
@test_throws TypeError readfasta("data/utf8_seq.fasta")

@test readfasta("data/spaces.fasta")
