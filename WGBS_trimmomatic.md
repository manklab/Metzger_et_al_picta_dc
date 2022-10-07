# Trimming WGBS data

Place all .fastq.gz files in the same directory so that they are the only .fastq.gz files in that directory.

Navigate to that directory and run the following code. This code is designed to identify unique "R1" files and then auto populate the corresponding R1, R1 and output file names for each sample.

To test that the code is working as designed, replace everything from java - MINLEN:50 with echo $f1 $f2. This should print every properly paired sequencing file.

```bash
for file in *fastq.gz ; do if [[ "$file" ==  *"R1"* ]] ; then name=`echo $file | grep -o "[A-Za-z].*_R"` ; f1=$name\1.fastq.gz ; f2=$name\2.fastq.gz; java -jar /Linux/Trimmomatic/trimmomatic.jar PE -threads 20 -phred33 -trimlog $name\trim.log $f1 $f2 $name\1_Paired.fastq.gz $name\1_Unpaired.fastq.gz $name\2_Paired.fastq.gz $name\2_Unpaired.fastq.gz ILLUMINACLIP:NEBNext-PE-methylated.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50; fi ; done;
```
