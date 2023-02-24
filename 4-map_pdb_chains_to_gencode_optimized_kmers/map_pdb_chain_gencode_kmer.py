# a more optimized way of aligning the sequences to gencode

#%%


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from difflib import SequenceMatcher
from collections import defaultdict

# Load the GENCODE protein sequences into a dictionary
gencode = {}
for record in SeqIO.parse("../3-map_pdb_chains_to_gencode/gencode.v43.pc_translations.fa", "fasta"):
    gencode[str(record.seq)] = record.id

# Define a function to check if a query sequence has a match in the GENCODE sequences
def check_match(query_seq):
    # Check all GENCODE sequences for exact matches or near matches with up to 3 mismatches
    matches = []
    for gencode_seq in gencode:
        start = gencode_seq.find(query_seq)
        if start >= 0:
            end = start + len(query_seq)
            matches.append((gencode[gencode_seq], start, end))
        if len(gencode_seq) != len(query_seq):
            continue
        mismatches = 0
        for i in range(len(query_seq)):
            if query_seq[i] != gencode_seq[i]:
                mismatches += 1
            if mismatches > 3:
                break
        if mismatches <= 3:
            matches.append((gencode[gencode_seq], None, None))
    return matches

# Load the PDB chain sequence and check for matches in the GENCODE sequences
matches = defaultdict(list)
for record in SeqIO.parse("../1-get_unique_sequences/pdb_human_seqs_unique.fasta", "fasta"):
    for match in check_match(str(record.seq)):
        if match[0] not in matches[str(record.id)]:
            matches[str(record.id)].append(match)

# Output the matches for each PDB chain to a file
with open("pdb_chain_to_gencode_matches.tsv", "w") as f:
    for chain_id in matches:
        f.write("Matches for PDB chain " + chain_id + "\n")
        for match in matches[chain_id]:
            f.write("Transcript name: " + match[0] + "\n")
            if match[1] is not None and match[2] is not None:
                f.write("Start position: " + str(match[1]) + "\n")
                f.write("End position: " + str(match[2]) + "\n")
            f.write("\n")

#%%

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from difflib import SequenceMatcher
from collections import defaultdict

# Load the GENCODE protein sequences into a dictionary
gencode = {}
for record in SeqIO.parse("../3-map_pdb_chains_to_gencode/gencode.v43.pc_translations.fa", "fasta"):
    gencode[str(record.seq)] = record.id

# Define a function to create a k-mer index for the protein sequences
def create_kmer_index(k, seqs):
    index = defaultdict(list)
    for seq in seqs:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            index[kmer].append(seq)
    return index

# Create a k-mer index for the protein sequences
k = 8
index = create_kmer_index(k, gencode.keys())

# Define a function to check if a query sequence has a match in the GENCODE sequences
def check_match(query_seq):
    # Find all protein sequences that contain at least one k-mer that matches a k-mer in the query sequence
    potential_matches = set()
    for i in range(len(query_seq) - k + 1):
        kmer = query_seq[i:i+k]
        for match in index[kmer]:
            potential_matches.add(match)

    # Check potential matches for exact matches or near matches with up to 3 mismatches
    for gencode_seq in potential_matches:
        start = gencode_seq.find(query_seq)
        if start >= 0:
            end = start + len(query_seq)
            return gencode[gencode_seq], start, end
        if len(gencode_seq) != len(query_seq):
            continue
        mismatches = 0
        for i in range(len(query_seq)):
            if query_seq[i] != gencode_seq[i]:
                mismatches += 1
            if mismatches > 3:
                break
        if mismatches <= 3:
            return gencode[gencode_seq], None, None
    return None, None, None

# Load the query sequences and check for matches in the GENCODE sequences
matches = {}
for record in SeqIO.parse("../1-get_unique_sequences/pdb_human_seqs_unique.fasta", "fasta"):
    match, start, end = check_match(str(record.seq))
    if match is not None:
        if str(record.id) in matches:
            matches[str(record.id)].append((match, start, end))
        else:
            matches[str(record.id)] = [(match, start, end)]

# # Output the matches for each query sequence
# for query_id in matches:
#     print("Matches for query sequence", query_id)
#     for match in matches[query_id]:
#         print("Transcript name:", match[0])
#         if match[1] is not None and match[2] is not None:
#             print("Start position:", match[1])
#             print("End position:", match[2])
#         print()

#%%
# Output the matches for each query sequence to a file
with open("pdb_chain_to_gencode_matches.tsv", "w") as f:
    for query_id in matches:
        f.write("Matches for query sequence " + query_id + "\n")
        for match in matches[query_id]:
            f.write("Transcript name: " + match[0] + "\n")
            if match[1] is not None and match[2] is not None:
                f.write("Start position: " + str(match[1]) + "\n")
                f.write("End position: " + str(match[2]) + "\n")
            f.write("\n")

#%%