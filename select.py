#1/usr/bin/env python3
import re
import fileinput
from collections import defaultdict

gene_dict = defaultdict(list)
tx_dict = defaultdict(list)
tx_pos_dict = defaultdict(list)
f=open("/home/galaxy/Desktop/Oryza_sativa.IRGSP-1.0.41.gff3")
f.new=open("/home/galaxy/Desktop/Oryza_sativa.IRGSP-1.0.41.gff3.new","w")
for line in f:
  if line.startswith("#"):
    continue
  content = line.split("\t")
  if content[2] == 'gene':
    gene_id = re.search(r'ID=(.*?)[;\n]',content[8]).group(1)
    gene_dict[gene_id] = []
  if content[2] == 'transcript' or content[2] == 'mRNA':
    tx_id = re.search(r'ID=(.*?)[;\n]',content[8]).group(1)
    tx_parent = re.search(r'Parent=(.*?)[;\n]',content[8]).group(1)
    gene_dict[tx_parent].append(tx_id)
    tx_pos_dict[tx_id] = [content[0],content[3], content[4], content[6]]
  if content[2] == 'CDS':
    width = int (content[4]) - int(content[3])
    cds_parent = re.search(r'Parent=(.*?)[;\n]',content[8]).group(1)
    tx_dict[cds_parent].append(width)

for gene, txs in gene_dict.items():
    tmp = 0
    for tx in txs:
       tx_len = sum(tx_dict[tx])
       if tx_len > tmp:
           lst_tx = tx
           tmp = tx_len
    tx_chrom = tx_pos_dict[lst_tx][0]
    tx_start = tx_pos_dict[lst_tx][1]
    tx_end  = tx_pos_dict[lst_tx][2]
    tx_strand = tx_pos_dict[lst_tx][3]
    f.new.write("{gene}\t{tx}\t{chrom}\t{start}\t{end}\t{strand}\n".format(gene=gene,tx=lst_tx,chrom=tx_chrom,start=tx_start,end=tx_end,strand=tx_strand))

f.new.close()
