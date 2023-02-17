# cluster sequences by 95 similarity

#%%

import os

# Set the input and output file paths
input_file = '../1-get_unique_sequences/pdb_human_seqs_unique.fasta'
output_file = 'protein_clusters.fasta'

# Set the CD-HIT parameters
c = 0.95  # Cluster sequences at 95% identity
cmd = f'cd-hit -i {input_file} -o {output_file} -c {c}'

# Run the CD-HIT command
os.system(cmd)

# Print the number of clusters and the clustering results
with open(output_file + '.clstr', 'r') as f:
    cluster_dict = {}
    for line in f:
        line = line.strip()
        if line.startswith('>Cluster'):
            cluster_id = line.split()[1]
            cluster_dict[cluster_id] = []
        else:
            seq_id = line.split('>')[1].split('...')[0]
            cluster_dict[cluster_id].append(seq_id)

    num_clusters = len(cluster_dict)

print(f'Number of clusters: {num_clusters}')

with open('pdb_seq_clusters_95.tsv', 'w') as ofile:
	for cluster_id, seq_ids in cluster_dict.items():
		ofile.write(f'Cluster {cluster_id}: {seq_ids}\n')
# 	print(f'Cluster {cluster_id}: {seq_ids}')


#%%


