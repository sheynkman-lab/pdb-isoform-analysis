from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from difflib import SequenceMatcher

# Load the GENCODE protein sequences into a dictionary
gencode = {}
for record in SeqIO.parse("gencode_proteins.fa", "fasta"):
    gencode[str(record.seq)] = record.id

# Define a function to check if a query sequence has a match in the GENCODE sequences
def check_match(query_seq):
    for gencode_seq in gencode:
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
for record in SeqIO.parse("query.fasta", "fasta"):
    match, start, end = check_match(str(record.seq))
    if match is not None:
        if str(record.id) in matches:
            matches[str(record.id)].append((match, start, end))
        else:
            matches[str(record.id)] = [(match, start, end)]

# Output the matches for each query sequence
for query_id in matches:
    print("Matches for query sequence", query_id)
    for match in matches[query_id]:
        print("Transcript name:", match[0])
        if match[1] is not None and match[2] is not None:
            print("Start position:", match[1])
            print("End position:", match[2])
        print()
