#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tuesday May 2, 2023

@author: Erick Santiago

    This script two samples (hybrid and its closest parental) of a phylip file and counts the number of SNPs taking into account IUPAC ambiguity codes.

"""

import pandas as pd
import sys
import os,csv
import subprocess
import re

inputs=sys.argv
hybrid=inputs[1] + " "
input_file=inputs[2]
parental=inputs[3] + " "
species=inputs[4]

ldf=[]
i=0
with open(input_file) as file:
    for line in file.readlines():
        fmtdrow=[]
        if (i==0):
            i+=1
            continue
        else:
            row=line.split()
            SNPS=list(row[1])
            ID=row[0] + " "
            fmtdrow.append(ID)
            fmtdrow.extend(SNPS)
            ldf.append(fmtdrow)
    df=pd.DataFrame(ldf)

#Subset dataframe to only the hybrid and its closest parental strain
subset=df.loc[df.iloc[:,0] == hybrid]
parentalrow=df.loc[df.iloc[:,0] == parental]
subset= subset.append(parentalrow,ignore_index=True)


#Transpose de dataframe to keep samples as columns and SNPs as rows
df_t=subset.transpose()
df_t.columns=df_t.iloc[0]
df_t=df_t[1:]
df_t.reset_index(inplace=True,drop=True)

#Create dataframe to store counts
init=[0]*16
counts=pd.DataFrame(
    {
        "A":init,
        "C":init,
        "G":init,
        "T":init,
        "R":init,
        "Y":init,
        "S":init,
        "W":init,
        "K":init,
        "M":init,
        "B":init,
        "D":init,
        "H":init,
        "V":init,
        "N":init,
        "-":init,
    },
    index=["A","C","G","T","R","Y","S","W","K","M","B","D","H","V","N","-"]
)
aliases={"a":"A","t":"T","g":"G","c":"C"} #Dictionary of aliases to read lowercased nucleotides
counts.index.name="Hybrid_SNPs"

#Create matrix to weigh the SNP count
weights=pd.DataFrame(
    {
        "A":[1,0,0,0,0.5,0,0,0.5,0,0.5,0,0.33,0.33,0.33,0,0],
        "C":[0,1,0,0,0,0.5,0.5,0,0,0.5,0.33,0,0.33,0.33,0,0],
        "G":[0,0,1,0,0.5,0,0.5,0,0.5,0,0.33,0.33,0,0.33,0,0],
        "T":[0,0,0,1,0,0.5,0,0.5,0.5,0,0.33,0.33,0.33,0,0,0],
        "R":[1,0,1,0,1,0,0.5,0.5,0.5,0.5,0.33,0.67,0.33,0.67,0,0],
        "Y":[0,1,0,1,0,1,0.5,0.5,0.5,0.5,0.67,0.33,0.67,0.33,0,0],
        "S":[0,1,1,0,0.5,0.5,1,0,0.5,0.5,0.66,0.33,0.33,0.66,0,0],
        "W":[1,0,0,1,0.5,0.5,0,1,0.5,0.5,0.33,0.67,0.67,0.33,0,0],
        "K":[0,0,1,1,0.5,0.5,0.5,0.5,1,0,0.67,0.67,0.33,0.33,0,0],
        "M":[1,1,0,0,0.5,0.5,0.5,0.5,0,1,0.33,0.33,0.33,0.67,0,0],
        "B":[0,1,1,1,0.5,1,1,0.5,1,0.5,1,0.67,0.67,0.67,0,0],
        "D":[1,0,1,1,1,0.5,0.5,1,1,0.5,0.67,1,0.67,0.67,0,0],
        "H":[1,1,0,1,0.5,1,0.5,1,0.5,1,0.67,0.67,1,0.67,0,0],
        "V":[1,1,1,0,1,0.5,1,0.5,0.5,1,0.67,0.67,0.67,1,0,0],
        "N":init,
        "-":init,
    },
    index=["A","C","G","T","R","Y","S","W","K","M","B","D","H","V","N","-"]
)

#Read the hybrid sequence, compare it to the parental sequence to store the number of matches in the counts dataframe
i=0
for index,row in df_t.iterrows():
    if(i==0):
        i+=1
        continue
    if(row[0] in aliases.keys()):
        key_hyb=aliases.get(row[0])
    else:
        key_hyb=row[0]

    if(row[1] in aliases.keys()):
        key_parent=aliases.get(row[1])
    else:
        key_parent=row[1]

    counts.at[key_hyb,key_parent]+=1
    i+=1

#Multiply the snp count by the matrix of values and sum all the rows to obtain the final number of SNPs
total_snps=counts.mul(weights)
sum=total_snps.sum()
print(sum.sum())


#Guardar la matriz ponderada también
weighted_counts="/mnt/Timina/lmorales/Public/ymez/data/phylogeny/esantiago_HybTrees/" + species + "/" + hybrid.strip() + "/" + species + "_" + parental.strip() + "_WeightedSNPCount_" + "HYB_" + hybrid.strip() + ".tsv"
counts.to_csv(weighted_counts,header=True,sep='\t',doublequote=False)


raw_counts="/mnt/Timina/lmorales/Public/ymez/data/phylogeny/esantiago_HybTrees/" + species + "/" + hybrid.strip() + "/" + species + "_" + parental.strip() + "_RawSNPCount_" + "HYB_" + hybrid.strip() + ".tsv"
counts.to_csv(raw_counts,header=True,sep='\t',doublequote=False)
