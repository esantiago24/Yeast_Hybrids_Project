import os
from Bio import SeqIO

#Define function to specify dictionary keys to be the first element in the FASTA identifier
def get_acc(identifier):
    parts = identifier.split("|")
    assert len(parts) == 3
    return parts[0]

#Define directories
cds_dir="/mnt/Timina/lmorales/Public/ymez/esantiago/FindClosestRef/data/GenSeqs/CDs/"
genids_dir="/mnt/Timina/lmorales/Public/ymez/esantiago/FindClosestRef/data/GenSeqs/GenIDs/"

#Iterate over the CDs files directory in order to fetch the nt sequences using the gene ID lists from the GenIDs directory
for file in os.listdir(cds_dir):
    strain=file[:-7]
    newdir="/mnt/Timina/lmorales/Public/ymez/esantiago/FindClosestRef/data/GenSeqs/nt_seqs_50/" + strain +"/"
    os.mkdir(newdir) #Create a directory with the strain name to store all of the individual gene fasta files
    seqs_index=SeqIO.index(cds_dir + file,"fasta",key_function=get_acc) #Index the cds fasta file to create a rapidly parseable dictionary of the nt sequences
    gene_ids_file=open(genids_dir + strain + "_50genes.txt")
    genes_of_interest=gene_ids_file.readlines()
    for gi in genes_of_interest: #Iterate over the list of genes of interest, fetch the nt sequence from the dictionary and write each out in a file
        gi=gi.strip()
        gene=seqs_index[gi]
        SeqIO.write(gene,newdir + gi + ".fa","fasta")
