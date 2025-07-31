#!/bin/bash
#
#SBATCH --job-name=varientcall_scotpearl	### Replace with the name of your job
#SBATCH --mail-user=v.gillman.21@abdn.ac.uk	### Put your email here
#SBATCH --partition=uoa-compute
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/uoa/home/r01vg21/sharedscratch/snpsaurus_margmarg_RADseq/scotpearl_popgen/output/slurm_out/investigateoutliersnps_%A_scotpearl.out		# %A is the master job ID
#SBATCH --mem=15G			# If over 25 than goes into high mem node

date

#file location
dir="/uoa/home/r01vg21/sharedscratch/snpsaurus_margmarg_RADseq/scotpearl_popgen/output"
combiref="/uoa/scratch/users/r01vg21/snpsaurus_margmarg_RADseq/pearlmusselarchive/CombinedReference.fna"

module load bedtools2/2.30.0
module load blast-plus/2.14.1

# bedtools getfasta \
  # -fi $combiref \
  # -bed ${dir}/loc_updated_snps/rda_qvalue_outliers_snp_flanks_100bp.bed \
  # -name -fo ${dir}/loc_updated_snps/rda_qvalue_outliers_100bp.fa
  
#for size in 100 500 1000; do
for size in  1000; do
 bedfile="${dir}/loc_updated_snps/rda_qvalue_outliers_snp_flanks_${size}bp.bed"
  fasta_output="${dir}/loc_updated_snps/rda_qvalue_outliers_${size}bp.fa"
  blast_full="${dir}/loc_updated_snps/blast_results_${size}bp_full.xml"
  blast_mollusc="${dir}/loc_updated_snps/blast_results_${size}bp_mollusc.xml"

#  echo "Extracting FASTA for ${size}bp flanks..."
#  bedtools getfasta -fi $combiref -bed $bedfile -name -fo $fasta_output

#  echo "Running full BLAST for ${size}bp..."
#  blastn -query "$fasta_output" -db nt -remote -out "$blast_full" -outfmt 5

  echo "Running molluscan BLAST for ${size}bp..."
  blastn -query "$fasta_output" -db nt -remote -entrez_query "txid6447[Organism]" -out "$blast_mollusc" -outfmt 5
done

