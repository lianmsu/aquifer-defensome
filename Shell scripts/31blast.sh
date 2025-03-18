
mkdir -p /home/lipengwei2024phd/03GW/12minced/blastdb

makeblastdb -in /home/lipengwei2024phd/03GW/12minced/2virus_seq/all_virus_seq.fna -dbtype nucl -out /home/lipengwei2024phd/03GW/12minced/blastdb/prophage.db -hash_index -max_file_sz '4GB'

#!/bin/bash
#SBATCH --job-name=blastn
#SBATCH --partition=cpu 
#SBATCH --nodes=1 
#SBATCH --get-user-env
#SBATCH --ntasks-per-node=12
#SBATCH --mem=125GB 
#SBATCH --output=job-prophage.%j.%N.out

module load BLAST+/2.11.0-gompi-2021a

blastn -db /home/lipengwei2024phd/03GW/12minced/blastdb/all_virus_seq.db -query bin_spacers.fa -out /home/lipengwei2024phd/03GW/12minced/2virus_seq/spacers_prophage.blastn.out -outfmt '6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send' -num_threads 12 -evalue 1e-5

