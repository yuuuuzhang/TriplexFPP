import numpy as np
import pandas as pd
import math
import random
from Bio import SeqIO

## generate random DNA sites in promotor

# read promotor sites
promotor = pd.read_csv('.../ensembl_promotor.txt', delimiter='\t')
chrom_p = promotor.Chromosome
start_p = promotor.Start
end_p = promotor.End

# read triplexDNA site
dna = []
records = list(SeqIO.parse(".../triplexDNA_2546.fa", "fasta"))
name = []
for j in range(len(records)):
    name.append(records[j].id)

def findnth(haystack, needle, n):
    parts= haystack.split(needle, n+1)
    if len(parts)<=n+1:
        return -1
    return len(haystack)-len(parts[-1])-len(needle)

chrom_dna = []
start_dna = []
end_dna = []
for i in range(len(name)):
    a = findnth(str(name[i]),":",1)
    b = findnth(str(name[i]),":",2)
    c = findnth(str(name[i]),":",3)
    d = findnth(str(name[i]),":",4)
    chrom_dna.append(name[i][a+1:b])
    start_dna.append(name[i][b+1:c])
    end_dna.append(name[i][c+1:d])
#print (len(chrom_dna))

dna_len = [int(end_dna[i])-int(start_dna[i])+1 for i in range(len(end_dna))]
promotor_len = [int(end_p[i])-int(start_p[i])+1 for i in range(len(end_p))]

l = len(promotor_len)
start_r = []
end_r = []
chrom_r = []
for i in range(len(dna_len)):
    while True:
        id_r = random.randint(0, l) # random select a promoter
        if promotor_len[id_r]>dna_len[i]:
            r_t = int(start_p[id_r])+random.randint(0, promotor_len[id_r]-dna_len[i])
            start_r.append(r_t)
            end_r.append(int(r_t)+int(dna_len[i])-1)
            chrom_r.append(chrom_p[id_r])
            break
#print (len(dna_len), len(end_r), dna_len[0], end_r[0]-start_r[0]+1)

# exclude the not wanted cases, usually this would not happen
for i in range(len(chrom_r)):
    for j in range(len(dna_len)):
        if chrom_r[i]==chrom_dna[j]:
            if abs(start_r[i]-int(start_dna[j]))<100 and abs(end_r[i]-int(end_dna[j]))<100:
                print (i,j)

# write to .bed file
rows = zip(chrom_r,start_r,end_r)
with open('.../DNA_neg2547_5.txt', "w") as f:
    writer = csv.writer(f)
    for row in rows:
        writer.writerow(row)

# obtain sequence
server = "https://rest.ensembl.org"
seq = []
for i in range(len(chrom_r)):
    chrom = chrom_r[i]
    loc1 = start_r[i]
    loc2 = end_r[i]
    
    ext = "/sequence/region/human/"+str(chrom)+":"+str(loc1)+".."+str(loc2)+":-1?"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
    
    if not r.ok:
        print(i)
    #r.raise_for_status()
    #sys.exit()
    seq.append(r.text)

# write negative sequence data to file
with open('.../DNA_neg_seq_5.fa', 'w') as f:
    f.write(''.join(seq))
