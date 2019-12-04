
# PIPELINE FOR PROCESSING SINGLE END RNA-seq DATA USING HISAT2(MAPPING) 
# AND CUFFLINKS(ASSEMBLY AND QUANTIFICATION)
	
	# PROJECT FOLDER SHOULD CONTAIN TWO FOLDERS:
		# fastq_files - containing SE .fastq.gz data
		# slurm-logs - for stdout/stderr output
	# GENOMES and other reference info is in ~/data/ for easy access; all indices and 
	# other reference file manipulations are saved here
		
	# Each of the following was run as a separate script on the Peloton Cluster (CSE, UC Davis)
	# using the following slurm queue submission:
		# sbatch -p (med/high) -n (cores) script

  #################################
  # RUN HISAT TO MAP SE FRAGMENTS #
  #################################


#!/bin/bash -l

#SBATCH -D ~/ # project directory
#SBATCH --mail-type=ALL # type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ncsierra@ucdavis.edu 
#SBATCH -o ~/projects/misc/slurm-logs/cufflinks-o%j.txt # %j saves job number to filename
#SBATCH -e ~/projects/misc/slurm-logs/cufflinks-e%j.txt
#SBATCH -J hs2map # job name
#SBATCH -t 60:00 # wall time

PATH=$PATH:~/tools/hisat2-2.1.0/hisat2

# Index the genome
cd ~/data
mkdir hisat2_genomeindex
hisat2-build Aurelia.Genome_v1.2_11-27-18.fasta hisat2_genomeindex/Aurelia.Genome_v1.2_11-27-18.index

# Run hisat2
cd ~/projects/misc
mkdir hs2_sam

cd fastq_files
for f in *.fastq.gz
do
 hisat2 --no-softclip \
   -x ~/data/hisat2_genomeindex/Aurelia.Genome_v1.2_11-27-18.index \
   -U ${f} \
   --summary-file summaryfile-${f/.fastq.gz}.txt \
   -S ../hs2_sam/${f/.fastq.gz}.sam
done

# Save individual summary files to a single file
cat summaryfile-*.txt >> ../hs2_sam/SummaryFile_hs2map.txt




  ########################################################
  # RUN CUFFLINKS TO ASSEMBLE FRAGMENTS INTO TRANSCRIPTS #
  ########################################################


#!/bin/bash -l

#SBATCH -D ~/projects/misc/
#SBATCH --mail-type=ALL
#SBATCH -o ~/projects/misc/slurm-logs/cufflinks-o%j.txt
#SBATCH -e ~/projects/misc/slurm-logs/cufflinks-e%j.txt
#SBATCH --mail-user=ncsierra@ucdavis.edu 
#SBATCH -J cufflinks
#SBATCH -t 50:00:00

PATH=$PATH:~/tools/cufflinks-2.2.1.Linux_x86_64

cd ~/projects/misc/
mkdir assembled_transcripts

#loop to sort
cd hs2_sam/
for i in *.sam
do 
 sort -k 3,3 -k 4,4n ${i} > ${i}.sorted
done

#now loop to run cufflinks on sorted data
for i in *.sam
do
 cufflinks -b ~/data/Aurelia.Genome_v1.2_11-27-18.fasta \
   -o ../assembled_transcripts/${i/.sam} ${i}.sorted
done




  ###########################################################
  # RUN CUFFMERGE TO ASSEMBLE AN EXPERIMENTAL TRANSCRIPTOME #
  ###########################################################


#!/bin/bash -l

#SBATCH -D ~/projects/misc/
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ncsierra@ucdavis.edu
#SBATCH -o ~/projects/misc/slurm-logs/cuffmerge-o%j.txt
#SBATCH -e ~/projects/misc/slurm-logs/cuffmerge-e%j.txt
#SBATCH -J cuffmerge
#SBATCH -t 50:00:00

PATH=$PATH:~/tools/cufflinks-2.2.1.Linux_x86_64

#Make the text file containing the list of assembly paths that inputs to cuffmerge
cd ~/projects/misc/hs2_sam
for f in *.sam
do 
 echo "~/projects/misc/assembled_transcripts/"${f/.sam}"/transcripts.gtf" >> ../assembled_transcripts/TEST.txt
done

cd ../assembled_transcripts

cuffmerge -g ~/data/Aurelia.Genome_Annotation_v1.2_11-28-18.gff3 --ref-sequence ~/data/Aurelia.Genome_v1.2_11-27-18.fasta assemblies_path.txt




  ################################################################
  # RUN CUFFQUANT TO QUANTIFY TRANSCRIPT ABUNDANCES PER DATA SET #
  ################################################################


#!/bin/bash -l

#SBATCH -D ~/projects/misc/
#SBATCH --mail-type=ALL
#SBATCH -o ~/projects/misc/slurm-logs/cuffquant-o%j.txt
#SBATCH -e ~/projects/misc/slurm-logs/cuffquant-e%j.txt
#SBATCH --mail-user=ncsierra@ucdavis.edu
#SBATCH -J cuffquant
#SBATCH -t 50:00:00

PATH=$PATH:~/tools/cufflinks-2.2.1.Linux_x86_64

cd ~/projects/misc/hs2_sam/

# Resort sam files
mkdir samtools_sort

for f in *.sam
do 
 samtools sort $f -T ${f}_temp -o samtools_sort/${f/.sam}_sorted.sam
done

# Run cuffquant on newly sorted sam files
mkdir ../cuffq_abundances

for f in *.sam #from still in hs2_sam
do
 cuffquant -u -o ../cuffq_abundances/ ~/data/*.gff3 samtools_sorted/*_sorted.sam
done

