#genome copy 
#AKR_J
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/AKR_J
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/AKR_J_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/AKR_J/AKR_J.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/AKR_J_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' AKR_J_chromosomes_MT_unplaced.fasta.join.tsv > AKR_J.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/AKR_J_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/AKR_J/AKR_J.bed -fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/AKR_J/AKR_J.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/AKR_J/AKR_J.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/AKR_J/AKR_J.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/AKR_J_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/AKR_J.genes.gtf
./piPipes install -g AKR_J -C & 

#A_J
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/A_J
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/A_J_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/A_J/A_J.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/A_J_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' A_J_chromosomes_MT_unplaced.fasta.join.tsv > A_J.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/A_J_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/A_J/A_J.bed -fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/A_J/A_J.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/A_J/A_J.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/A_J/A_J.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/A_J_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/A_J/A_J.genes.gtf
./piPipes install -g A_J -C


#129S1_SvImJ
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/129S1_SvImJ
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/129S1_SvImJ_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/129S1_SvImJ/129S1_SvImJ.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/129S1_SvImJ_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' 129S1_SvImJ_chromosomes_MT_unplaced.fasta.join.tsv  > 129S1_SvImJ.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/129S1_SvImJ_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/129S1_SvImJ/129S1_SvImJ.bed -fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/129S1_SvImJ/129S1_SvImJ.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/129S1_SvImJ/129S1_SvImJ.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/129S1_SvImJ/129S1_SvImJ.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/129S1_SvImJ_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/129S1_SvImJ.genes.gtf
./piPipes install -g 129S1_SvImJ -C

#"BALB_cJ"
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/BALB_cJ
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/BALB_cJ_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/BALB_cJ/BALB_cJ.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/BALB_cJ_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' BALB_cJ_chromosomes_MT_unplaced.fasta.join.tsv > BALB_cJ.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/BALB_cJ_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/BALB_cJ/BALB_cJ.bed -fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/BALB_cJ/BALB_cJ.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/BALB_cJ/BALB_cJ.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/BALB_cJ/BALB_cJ.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/BALB_cJ_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/BALB_cJ.genes.gtf
./piPipes install -g BALB_cJ -C



#"C3H_HeJ"
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C3H_HeJ
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/C3H_HeJ_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C3H_HeJ/C3H_HeJ.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/C3H_HeJ_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' C3H_HeJ_chromosomes_MT_unplaced.fasta.join.tsv > C3H_HeJ.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/C3H_HeJ_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C3H_HeJ/C3H_HeJ.bed -fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C3H_HeJ/C3H_HeJ.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C3H_HeJ/C3H_HeJ.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C3H_HeJ/C3H_HeJ.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/C3H_HeJ_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C3H_HeJ/C3H_HeJ.genes.gtf
./piPipes install -g C3H_HeJ -C


#"CAST_EiJ"
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CAST_EiJ
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/CAST_EiJ_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CAST_EiJ/CAST_EiJ.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/CAST_EiJ_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' CAST_EiJ_chromosomes_MT_unplaced.fasta.join.tsv  > CAST_EiJ.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/CAST_EiJ_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CAST_EiJ/CAST_EiJ.bed -fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CAST_EiJ/CAST_EiJ.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CAST_EiJ/CAST_EiJ.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CAST_EiJ/CAST_EiJ.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/CAST_EiJ_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CAST_EiJ/CAST_EiJ.genes.gtf
./piPipes install -g CAST_EiJ -C




#"CBA_J"
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CBA_J
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CBA_J
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/CBA_J_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CBA_J/CBA_J.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/CBA_J_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' CBA_J_chromosomes_MT_unplaced.fasta.join.tsv   > CBA_J.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/CBA_J_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CBA_J/CBA_J.bed-fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CBA_J/CBA_J.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CBA_J/CBA_J.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CBA_J/CBA_J.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/CBA_J_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/CBA_J/CBA_J.genes.gtf
./piPipes install -g CBA_J -C


#"DBA_2J"
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/DBA_2J
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/DBA_2J
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/DBA_2J_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/DBA_2J/DBA_2J.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/DBA_2J_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' DBA_2J_chromosomes_MT_unplaced.fasta.join.tsv   > DBA_2J.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/DBA_2J_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/DBA_2J/DBA_2J.bed -fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/DBA_2J/DBA_2J.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/DBA_2J/DBA_2J.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/DBA_2J/DBA_2J.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/DBA_2J_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/DBA_2J/DBA_2J.genes.gtf
./piPipes install -g DBA_2J -C

#"FVB_NJ"
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/FVB_NJ
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/FVB_NJ
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/FVB_NJ_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/FVB_NJ/FVB_NJ.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/FVB_NJ_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' FVB_NJ_chromosomes_MT_unplaced.fasta.join.tsv  > FVB_NJ.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/FVB_NJ_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/FVB_NJ/FVB_NJ.bed -fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/FVB_NJ/FVB_NJ.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/FVB_NJ/FVB_NJ.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/FVB_NJ/FVB_NJ.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/FVB_NJ_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/FVB_NJ/FVB_NJ.genes.gtf
./piPipes install -g FVB_NJ -C



#"LP_J"   #note:  using old genome 
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/LP_J
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/LP_J
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/LP_J_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/LP_J/LP_J.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/LP_J_chromosomes_MT_unplaced.fasta.out 
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/LP_J_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/LP_J/LP_J_chromosomes_MT_unplaced.fasta.join.tsv -fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/LP_J/LP_J.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/LP_J/LP_J.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/LP_J/LP_J.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/LP_J_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/LP_J.genes.gtf
./piPipes install -g LP_J -C



#"NOD_ShiLtJ"
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NOD_ShiLtJ
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NOD_ShiLtJ
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/NOD_ShiLtJ_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NOD_ShiLtJ/NOD_ShiLtJ.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/NOD_ShiLtJ_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' NOD_ShiLtJ_chromosomes_MT_unplaced.fasta.join.tsv  > NOD_ShiLtJ.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/NOD_ShiLtJ_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NOD_ShiLtJ/NOD_ShiLtJ.bed -fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NOD_ShiLtJ/NOD_ShiLtJ.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NOD_ShiLtJ/NOD_ShiLtJ.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NOD_ShiLtJ/NOD_ShiLtJ.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/NOD_ShiLtJ_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NOD_ShiLtJ/NOD_ShiLtJ.genes.gtf
./piPipes install -g NOD_ShiLtJ -C



#"NZO_HlLtJ",   #note: using old genome 
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NZO_HlLtJ
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NZO_HlLtJ
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/NZO_HlLtJ_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NZO_HlLtJ/NZO_HlLtJ.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/NZO_HlLtJ_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' NZO_HlLtJ_chromosomes_MT_unplaced.fasta.join.tsv  > NZO_HlLtJ.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/NZO_HlLtJ_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NZO_HlLtJ/NZO_HlLtJ.bed -fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NZO_HlLtJ/NZO_HlLtJ.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NZO_HlLtJ/NZO_HlLtJ.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NZO_HlLtJ/NZO_HlLtJ.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/NZO_HlLtJ_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/NZO_HlLtJ/NZO_HlLtJ.genes.gtf
./piPipes install -g NZO_HlLtJ -C


#"PWK_PhJ"
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/PWK_PhJ
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/PWK_PhJ
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/PWK_PhJ_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/PWK_PhJ/PWK_PhJ.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/PWK_PhJ_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' PWK_PhJ_chromosomes_MT_unplaced.fasta.join.tsv  > PWK_PhJ.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/PWK_PhJ_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/PWK_PhJ/PWK_PhJ.bed-fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/PWK_PhJ/PWK_PhJ.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/PWK_PhJ/PWK_PhJ.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/PWK_PhJ/PWK_PhJ.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/PWK_PhJ_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/PWK_PhJ/PWK_PhJ.genes.gtf
./piPipes install -g PWK_PhJ -C




#"SPRET_EiJ"
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/SPRET_EiJ
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/SPRET_EiJ
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/SPRET_EiJ_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/SPRET_EiJ/SPRET_EiJ.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/SPRET_EiJ_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' SPRET_EiJ_chromosomes_MT_unplaced.fasta.join.tsv   > SPRET_EiJ.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/SPRET_EiJ_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/SPRET_EiJ/SPRET_EiJ.bed -fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/SPRET_EiJ/SPRET_EiJ.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/SPRET_EiJ/SPRET_EiJ.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/SPRET_EiJ/SPRET_EiJ.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/SPRET_EiJ_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/SPRET_EiJ/SPRET_EiJ.genes.gtf
./piPipes install -g SPRET_EiJ -C




#"WSB_EiJ"
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/WSB_EiJ
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/WSB_EiJ
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/WSB_EiJ_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/WSB_EiJ/WSB_EiJ.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/WSB_EiJ_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' WSB_EiJ_chromosomes_MT_unplaced.fasta.join.tsv  > WSB_EiJ.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/WSB_EiJ_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/WSB_EiJ/WSB_EiJ.bed -fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/WSB_EiJ/WSB_EiJ.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/WSB_EiJ/WSB_EiJ.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/WSB_EiJ/WSB_EiJ.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/WSB_EiJ_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/WSB_EiJ/WSB_EiJ.genes.gtf
./piPipes install -g WSB_EiJ -C




#"C57BL_6NJ"
mkdir /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C57BL_6NJ
cd /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C57BL_6NJ
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/C57BL_6NJ_chromosomes_MT.fasta  /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C57BL_6NJ/C57BL_6NJ.fa
rmToUCSCTables.pl -out /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/repeatMasker/C57BL_6NJ_chromosomes_MT_unplaced.fasta.out 
awk '{if ($2 > $3) print $1"\t"$3"\t"$2; else print $1"\t"$2"\t"$3}' C57BL_6NJ_chromosomes_MT_unplaced.fasta.join.tsv  > C57BL_6NJ.bed
bedtools getfasta -fi /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/REL-2205-Assembly/C57BL_6NJ_chromosomes_MT.fasta  -bed /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C57BL_6NJ/C57BL_6NJ.bed -fo /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C57BL_6NJ/C57BL_6NJ.repBase.fa
cp /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C57BL_6NJ/C57BL_6NJ.repBase.fa /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C57BL_6NJ/C57BL_6NJ.transposon.fa
awk -F '\t' -v OFS='\t' '{ if ($3 == "gene") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'   /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/annotation/C57BL_6NJ_v3.2.gff3 > /mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/workflow/scripts/piPipes/common/C57BL_6NJ/C57BL_6NJ.genes.gtf
./piPipes install -g C57BL_6NJ -C


#repeatmasker file 
#conda install -c bioconda repeatmasker
rmToUCSCTables.pl -out ../../../resources/repeatMasker/AKR_J_chromosomes_MT_unplaced.fasta.out 
awk 'BEGIN{OFS="\t"} !/^#/ {chrom=$5; start=$6; end=$7+1; repeat=$11; print chrom, start, end, repeat}' repeatmasker.out > AKR_J_repeat.bed 

#output save as AKR_J_chromosomes_MT_unplaced.fasta.join.tsv
bedtools getfasta -fi ../../../resources/REL-2205-Assembly/AKR_J_chromosomes_MT.fasta  -bed AKR_J_repeat.bed -fo AKR_J.repBase.fa


#piRNA clusters 


#gff3 to gtf 

awk -F '\t' -v OFS='\t' '{ if ($3 == "exon") { gsub(/;/, "; ", $9); gsub(/=/, " \"", $9); $9 = $9 "\";"; print $0; } }'  ../../../../../resources/annotation/A_J_v3.2.gff3 > A_J_12.5dpp.genes.gtf


#pipies creting genome 

fq1 =  config['resultDir']  + "/cutadapt_se/{sid}-{tp}.{rep}.fastq"


piPipes small \
	-i config['resultDir']  + "/cutadapt_se/{sid}-{tp}.{rep}.fastq" \
	-g A_J \
	-o A_J_12.dpp.piPipes_out \
	-N uniqueXmiRNA \
	-c 12 \
	1> Zamore.SRA.ago3_het.ox.ovary.piPipes.stdout \
	2> Zamore.SRA.ago3_het.ox.ovary.piPipes.stderr