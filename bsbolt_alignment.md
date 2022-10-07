#BSBolt WGBS sequencing alignments

The following code needs to be run for each of the paired read files.

This code will:
1. Prepare the reference genome/database for WGBS alignment using bsbolt
2. Align WGBS files to the reference genome using bsbolt
3. Fixmates to prepare for duplicate removal
4. Sort bam by coordinates for duplicate calling
5. Remove duplicate reads
6. Index bam file for methylation calling
7. Call methylation values for CpG, CHG, and CHH loci

```bash
#bsbolt alignment
python3 -m bsbolt Align -t 15 -OT 20 -F1 path_to_R1_Paired.fastq.gz -F2 path_to_R2_Paired.fastq.gz -DB path_to_bsbolt_db -O output_file_name[bsbolt.bam]

# fixmates to prepare for duplicate removal, use -p to disable proper pair check
samtools fixmate -@ 20 -p -m bsbolt.bam bsbolt.fixmates.bam

# sort bam by coordinates for duplicate calling
samtools sort -@ 20 -o bsbolt.sorted.bam bsbolt.fixmates.bam

# remove duplicate reads
samtools markdup -@ 20 bsbolt.sorted.bam bsbolt.dup.bam

# index bam file for methylation calling
samtools index -@ 20 bsbolt.dup.bam

#methylation calling
python3 -m bsbolt CallMethylation -t 5 -min 1 -I bsbolt.dup.bam -DB path_to_bsbolt_db -O output_file_name[bsbolt_meth]
```
