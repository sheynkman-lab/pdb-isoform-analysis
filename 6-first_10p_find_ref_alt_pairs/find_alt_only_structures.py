import csv

# Read in the appris transcripts
appris_transcripts = set()
with open('../data/appris_transcript_names.tsv', 'r') as appris_file:
    appris_reader = csv.reader(appris_file, delimiter='\t')
    for row in appris_reader:
        appris_transcripts.add(row[1])

# Read in the nonappris transcripts
nonappris_transcripts = set()
with open('../data/not_appris_transcript_names.tsv', 'r') as nonappris_file:
    nonappris_reader = csv.reader(nonappris_file, delimiter='\t')
    for row in nonappris_reader:
        nonappris_transcripts.add(row[1])

# Read in the PDB to GenCode mapping
pdb_to_gencode = []
with open('../data/pdb_gencode_mapping.tsv', 'r') as pdb_file:
    pdb_reader = csv.reader(pdb_file, delimiter='\t')
    for row in pdb_reader:
        pdb_to_gencode.append(row)

# Group the PDB to GenCode mappings by qseqid
pdb_grouped_by_qseqid = {}
for row in pdb_to_gencode:
    qseqid = row[0]
    if qseqid not in pdb_grouped_by_qseqid:
        pdb_grouped_by_qseqid[qseqid] = []
    pdb_grouped_by_qseqid[qseqid].append(row)

# Filter the PDB to GenCode mapping based on appris and nonappris transcripts
filtered_pdb_to_gencode = []
for qseqid, pdb_group in pdb_grouped_by_qseqid.items():
    # Get the top hit based on bitscore
    top_hit = max(pdb_group, key=lambda row: float(row[9]))
    top_bitscore = float(top_hit[9])
    
    # Get all other hits with the same bitscore as the top hit
    other_hits = [row for row in pdb_group if float(row[9]) == top_bitscore and row != top_hit]
    
    # Determine if all sseqids are nonappris
    all_nonappris = all(row[4] in nonappris_transcripts for row in pdb_group)
    
    # If all sseqids are nonappris, add the gene name to the top hit and output to the filtered list
    if all_nonappris:
        top_hit_with_gene = top_hit + [top_hit[4].split('-')[0]]
        filtered_pdb_to_gencode.append(top_hit_with_gene)
        for other_hit in other_hits:
            other_hit_with_gene = other_hit + [other_hit[4].split('-')[0]]
            filtered_pdb_to_gencode.append(other_hit_with_gene)

# Write out the filtered mapping to a TSV file
with open('../data/filtered_pdb_gencode_mapping.tsv', 'w') as output_file:
    output_writer = csv.writer(output_file, delimiter='\t')
    for row in filtered_pdb_to_gencode:
        output_writer.writerow(row)
