### NEU - no dedup, broad

sbatch code7_cpg_enrichment/neu/no_dedup/individual_samples/broad/1_across_cpgs_bw_signal_1_rep/1_cpg_enrichment.sh 
sbatch code7_cpg_enrichment/neu/no_dedup/individual_samples/broad/2_across_cpgs_bw_signal_2_rep/1_cpg_enrichment.sh
sbatch code7_cpg_enrichment/neu/no_dedup/individual_samples/broad/3_across_peaks_bw_signal_1_rep/1_cpg_enrichment.sh 
sbatch code7_cpg_enrichment/neu/no_dedup/individual_samples/broad/4_across_peaks_bw_signal_2_rep/1_cpg_enrichment.sh  

sbatch code7_cpg_enrichment/neu/no_dedup/individual_samples/broad/1_across_cpgs_bw_signal_1_rep/2_combine_chunks.sh 
sbatch code7_cpg_enrichment/neu/no_dedup/individual_samples/broad/2_across_cpgs_bw_signal_2_rep/2_combine_chunks.sh
sbatch code7_cpg_enrichment/neu/no_dedup/individual_samples/broad/3_across_peaks_bw_signal_1_rep/2_combine_chunks.sh 
sbatch code7_cpg_enrichment/neu/no_dedup/individual_samples/broad/4_across_peaks_bw_signal_2_rep/2_combine_chunks.sh 


echo "Neu broad peaks enrichment:"
stat -c %y /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/cpg_enrichment/Neu/broad/cpg_enrichment_{1,2}_rep_in_peaks/chunk_0.csv
stat -c %y /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/cpg_enrichment/Neu/broad/cpg_enrichment_{1,2}_rep_in_peaks/cpg_enrichment_parallel.csv

echo -e "\nNeu broad CpG enrichment:"
stat -c %y /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/cpg_enrichment/Neu/broad/cpg_enrichment_{1,2}_rep_in_cpg/chunk_0.csv
stat -c %y /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/cpg_enrichment/Neu/broad/cpg_enrichment_{1,2}_rep_in_cpg/cpg_enrichment_parallel.csv

### NSC

sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/1_across_cpgs_bw_signal_1_rep/1_cpg_enrichment.sh 
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/2_across_cpgs_bw_signal_2_rep/1_cpg_enrichment.sh 
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/3_across_peaks_bw_signal_1_rep/1_cpg_enrichment.sh 
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/4_across_peaks_bw_signal_2_rep/1_cpg_enrichment.sh 

sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/1_across_cpgs_bw_signal_1_rep/2_combine_chunks.sh 
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/2_across_cpgs_bw_signal_2_rep/2_combine_chunks.sh 
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/3_across_peaks_bw_signal_1_rep/2_combine_chunks.sh 
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/4_across_peaks_bw_signal_2_rep/2_combine_chunks.sh

echo "NSC broad peaks enrichment:"
stat -c %y /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/cpg_enrichment/NSC/broad/cpg_enrichment_{1,2}_rep_in_peaks/chunk_0.csv
stat -c %y /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/cpg_enrichment/NSC/broad/cpg_enrichment_{1,2}_rep_in_peaks/cpg_enrichment_parallel.csv

echo -e "\nNSC broad CpG enrichment:"
stat -c %y /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/cpg_enrichment/NSC/broad/cpg_enrichment_{1,2}_rep_in_cpg/chunk_0.csv
stat -c %y /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/cpg_enrichment/NSC/broad/cpg_enrichment_{1,2}_rep_in_cpg/cpg_enrichment_parallel.csv