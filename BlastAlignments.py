import os

#Define directories
nt_seqs_dir="/mnt/Timina/lmorales/Public/ymez/esantiago/FindClosestRef/data/GenSeqs/nt_seqs_50/"
assemblies_dir="/mnt/Timina/lmorales/Public/ymez/esantiago/FindClosestRef/data/assemblies/"
alignments_dir="/mnt/Timina/lmorales/Public/ymez/esantiago/FindClosestRef/data/alignments/"

for assembly in os.listdir(assemblies_dir): #Iterate over the de novo assemblies
    assembly_name=assembly[:-6]
    new_dir=alignments_dir + assembly_name + "/"
    os.system("mkdir " + new_dir)    #Create a new directory per assembly to store their corresponding alignments
    for strain in os.listdir(nt_seqs_dir):  #Iterate over each strain
        strain_dir=nt_seqs_dir + strain + "/"
        for file in os.listdir(strain_dir): #Iterate over the 50 fasta files for each strain.
            blst_command="blastn -query " + strain_dir + file + " -subject " + assemblies_dir + assembly + " -outfmt '6 qaccver saccver qlen pident length mismatch qstart qend sstart send evalue bitscore' -max_target_seqs 1 -out " + new_dir + assembly_name + "_" + file[:-3] + ".txt"
            os.system(blst_command)
