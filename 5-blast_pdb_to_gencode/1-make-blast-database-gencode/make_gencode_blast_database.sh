# copy the gencode fasta in the local directory so the blast indices reside there
cp ../b_gencode_fasta_same_prot_clustered.v35.fa ./a_gencode_fasta_same_prot_clustered.v35.fa
makeblastdb -in a_gencode_fasta_same_prot_clustered.v35.fa -dbtype prot
