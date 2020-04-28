## MAIN STEPS FOR RNA SEQUENCING

# analysis directory: /Users/dmarti14/Documents/MRC_Postdoc/Projects/Communities/RNA_seq

# this will donwload the cDNA data from C. elegans from the Ensemblgenomes website 
# this is the target transcriptome
# https://www.ensembl.org/Caenorhabditis_elegans/Info/Index
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-44/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz -o celegans.fa.gz

# generate the index of the target transcriptome
salmon index -t celegans.fa.gz -i celegans_index

# OP50
salmon quant -i celegans_index -l A -1 data/Unaligned/N2GCB17_S59_L005_R1_001.fastq.gz -2 data/Unaligned/N2GCB17_S59_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/OP50_1
salmon quant -i celegans_index -l A -1 data/Unaligned/N2GCB20_S63_L005_R1_001.fastq.gz -2 data/Unaligned/N2GCB20_S63_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/OP50_2
salmon quant -i celegans_index -l A -1 data/Unaligned/N2GCB26_S67_L005_R1_001.fastq.gz -2 data/Unaligned/N2GCB26_S67_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/OP50_3
salmon quant -i celegans_index -l A -1 data/Unaligned/N2GCB29_S69_L005_R1_001.fastq.gz -2 data/Unaligned/N2GCB29_S69_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/OP50_4

# OP50_PDXJK
salmon quant -i celegans_index -l A -1 data/Unaligned/N2OP5017_S58_L005_R1_001.fastq.gz -2 data/Unaligned/N2OP5017_S58_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/OP50_PDXJK_1
salmon quant -i celegans_index -l A -1 data/Unaligned/N2OP5020_S60_L005_R1_001.fastq.gz -2 data/Unaligned/N2OP5020_S60_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/OP50_PDXJK_2
salmon quant -i celegans_index -l A -1 data/Unaligned/N2OP5026_S64_L005_R1_001.fastq.gz -2 data/Unaligned/N2OP5026_S64_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/OP50_PDXJK_3
salmon quant -i celegans_index -l A -1 data/Unaligned/N2OP5029_S68_L005_R1_001.fastq.gz -2 data/Unaligned/N2OP5029_S68_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/OP50_PDXJK_4


# MG1655
salmon quant -i celegans_index -l A -1 data/Unaligned/ep2GCB17_S57_L005_R1_001.fastq.gz -2 data/Unaligned/ep2GCB17_S57_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/MG1655_1
salmon quant -i celegans_index -l A -1 data/Unaligned/ep2GCB20_S62_L005_R1_001.fastq.gz -2 data/Unaligned/ep2GCB20_S62_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/MG1655_2
salmon quant -i celegans_index -l A -1 data/Unaligned/ep2GCB26_S65_L005_R1_001.fastq.gz -2 data/Unaligned/ep2GCB26_S65_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/MG1655_3
salmon quant -i celegans_index -l A -1 data/Unaligned/ep2GCB29_S71_L005_R1_001.fastq.gz -2 data/Unaligned/ep2GCB29_S71_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/MG1655_4


# MG1655_PDXJK
salmon quant -i celegans_index -l A -1 data/Unaligned/ep2OP5017_S56_L005_R1_001.fastq.gz -2 data/Unaligned/ep2OP5017_S56_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/MG1655_PDXJK_1
salmon quant -i celegans_index -l A -1 data/Unaligned/ep2OP5020_S61_L005_R1_001.fastq.gz -2 data/Unaligned/ep2OP5020_S61_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/MG1655_PDXJK_2
salmon quant -i celegans_index -l A -1 data/Unaligned/ep2OP5026_S66_L005_R1_001.fastq.gz -2 data/Unaligned/ep2OP5026_S66_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/MG1655_PDXJK_3
salmon quant -i celegans_index -l A -1 data/Unaligned/ep2OP5029_S70_L005_R1_001.fastq.gz -2 data/Unaligned/ep2OP5029_S70_L005_R2_001.fastq.gz -p 8 --validateMappings --gcBias -o quants/MG1655_PDXJK_4
