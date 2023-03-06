# get the list of the genes and isoforms that represent the appris reference
# as well as a list of isoforms that represent the non-appris reference
#%% 


# import packages
from collections import defaultdict 
import os, sys


#input direcotory
gtf_file = "../data/gencode.v43.annotation.gtf"

#gene _>then get the appris id and the isoname 
appris = defaultdict(list)
for line in open(gtf_file): 
    if line.startswith('#'): continue
    wds  = line.split('\t')
    feat = wds[2]
    if feat == 'transcript':
        acc = wds[8]
        gene = acc.split('gene_name "')[1].split('"')[0]
        isoname = acc.split('transcript_name "')[1].split('"')[0]
        if 'appris' in acc:
            appris_tag = acc.split('appris_')[1].split('";')[0]
            # later, will determine appris principle based on alphanum sort
            # therefore, change alternative to x_alt.
            appris_tag = appris_tag.replace("alternative","2_alternative")
            appris_tag = appris_tag.replace("principal","1_principal")
            appris[gene].append([appris_tag, isoname])
        else:
            appris_tag = '3_no_appris_annotation'
            # for some reason, all ENSPs appear to be lncRNAs?
            # additional 40K entries above the 20K PC?
            if gene.startswith('ENSG'): continue
            appris[gene].append([appris_tag, isoname])

def get_the_top_appris_isoname(appris_info) : 
    best_appris_tag = sorted(appris_info)[0][0]
    isonames_w_best_appris_tag = []
    isonames_not_best_appris_tag = []
    #there can be 2 reference names sometimes so we want to get them both 
    for appris_tag, isoname in appris_info: 
        if appris_tag == best_appris_tag: 
            isonames_w_best_appris_tag.append(isoname)
        else:
            isonames_not_best_appris_tag.append(isoname)
    return best_appris_tag, isonames_w_best_appris_tag, isonames_not_best_appris_tag

#%%
#write the transcripts to an output file 
with open('../output/appris_ref_transcripts.tsv', 'w') as ofile: 
    ofile.write('gene\tbest_appris_tag_for_gene\tappris_reference_transcript_names\tnot_appris_reference_transcript_names\n')
    for gene, appris_info in appris.items(): 
        best_appris_tag, appris_isonames, nonappris_isonames = get_the_top_appris_isoname(appris_info)
        ofile.write('{}\t{}\t{}\t{}\n'.format(gene, best_appris_tag, ','.join(appris_isonames), ','.join(nonappris_isonames))) 

         
#%%

# write out list of appris transcript names
with open('../output/appris_transcript_names.tsv', 'w') as ofile:
    for gene, appris_info in appris.items(): 
        best_appris_tag, appris_isonames, nonappris_isonames = get_the_top_appris_isoname(appris_info)
        for appris_isoname in appris_isonames:
            ofile.write(gene + '\t' + appris_isoname + '\n')

# write out list of non-appris transcript names
with open('../output/not_appris_transcript_names.tsv', 'w') as ofile:
    for gene, appris_info in appris.items(): 
        best_appris_tag, appris_isonames, nonappris_isonames = get_the_top_appris_isoname(appris_info)
        for appris_isoname in nonappris_isonames:
            ofile.write(gene + '\t' + appris_isoname + '\n')

