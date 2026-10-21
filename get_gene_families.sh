# Location to use: /mnt/das1/TFM/results/pangenome/get_homologues/
# Usage: bash /mnt/das1/TFM/scr/get_gene_families.sh
# Output: 	One text file per gene family with the IDs of the genes it contains in the directory: 
# 				/mnt/das1/TFM/results/pangenome/get_homologues/families/
# 			A tsv file with the name of the gene family in the first column and the name of the gene in the second column in:
# 				/mnt/das1/TFM/results/pangenome/get_homologues/gene_families.tsv

mkdir families

# Obtain the gene IDs from the FASTA headers of the gene families FASTAs and put them in a text file:
for i in data_gbks_homologues/SalinibacterruberSP3059prokka_f0_0taxa_algOMCL_e0_/*.faa
do 
base=$(basename $i .faa)
grep ">" $i | cut -d':' -f2 | cut -d' ' -f1 > families/${base}.txt
done

# Print the name of the family, a tab separator and the name of the gene in a tsv file:
for i in families/*.txt
do
base=$(basename $i .txt)
cat $i | while read line 
do
echo -e $base$'\t'$line >> gene_families.tsv
done
done
