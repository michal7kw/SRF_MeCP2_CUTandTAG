sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/flanking/1_cpg_enrichment_consistent_peaks.sh
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/flanking/2_combine_chunks.sh 
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/flanking/3_split_the_list.sh 
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/flanking/4_peaks_localization.sh
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/flanking/5_analyze_cpg_genes.sh
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/flanking/6_analyze_coding_sequences.sh

sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/narrow/flanking/1_cpg_enrichment_consistent_peaks.sh
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/narrow/flanking/2_combine_chunks.sh 
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/narrow/flanking/3_split_the_list.sh 
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/narrow/flanking/4_peaks_localization.sh
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/narrow/flanking/5_analyze_cpg_genes.sh
sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/narrow/flanking/6_analyze_coding_sequences.sh

#########################################################################################################################

sbatch code7_cpg_enrichment/nsc/dedup/individual_samples/broad/flanking/1_cpg_enrichment_consistent_peaks.sh
sbatch code7_cpg_enrichment/nsc/dedup/individual_samples/broad/flanking/2_combine_chunks.sh 
sbatch code7_cpg_enrichment/nsc/dedup/individual_samples/broad/flanking/3_split_the_list.sh 
sbatch code7_cpg_enrichment/nsc/dedup/individual_samples/broad/flanking/4_peaks_localization.sh
sbatch code7_cpg_enrichment/nsc/dedup/individual_samples/broad/flanking/5_analyze_cpg_genes.sh
sbatch code7_cpg_enrichment/nsc/dedup/individual_samples/broad/flanking/6_analyze_coding_sequences.sh

sbatch code7_cpg_enrichment/nsc/dedup/individual_samples/narrow/flanking/1_cpg_enrichment_consistent_peaks.sh
sbatch code7_cpg_enrichment/nsc/dedup/individual_samples/narrow/flanking/2_combine_chunks.sh 
sbatch code7_cpg_enrichment/nsc/dedup/individual_samples/narrow/flanking/3_split_the_list.sh 
sbatch code7_cpg_enrichment/nsc/dedup/individual_samples/narrow/flanking/4_peaks_localization.sh
sbatch code7_cpg_enrichment/nsc/dedup/individual_samples/narrow/flanking/5_analyze_cpg_genes.sh
sbatch code7_cpg_enrichment/nsc/dedup/individual_samples/narrow/flanking/6_analyze_coding_sequences.sh

#########################################################################################################################

job1=$(sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/flanking/1_cpg_enrichment_consistent_peaks.sh | cut -d ' ' -f 4)
job2=$(sbatch --dependency=afterok:$job1 code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/flanking/2_combine_chunks.sh | cut -d ' ' -f 4)
job3=$(sbatch --dependency=afterok:$job2 code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/flanking/3_split_the_list.sh | cut -d ' ' -f 4)
job4=$(sbatch --dependency=afterok:$job3 code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/flanking/4_peaks_localization.sh | cut -d ' ' -f 4)
job5=$(sbatch --dependency=afterok:$job4 code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/flanking/5_analyze_cpg_genes.sh | cut -d ' ' -f 4)
sbatch --dependency=afterok:$job5 code7_cpg_enrichment/nsc/no_dedup/individual_samples/broad/flanking/6_analyze_coding_sequences.sh

job1=$(sbatch code7_cpg_enrichment/nsc/no_dedup/individual_samples/narrow/flanking/1_cpg_enrichment_consistent_peaks.sh | cut -d ' ' -f 4)
job2=$(sbatch --dependency=afterok:$job1 code7_cpg_enrichment/nsc/no_dedup/individual_samples/narrow/flanking/2_combine_chunks.sh | cut -d ' ' -f 4)
job3=$(sbatch --dependency=afterok:$job2 code7_cpg_enrichment/nsc/no_dedup/individual_samples/narrow/flanking/3_split_the_list.sh | cut -d ' ' -f 4)
job4=$(sbatch --dependency=afterok:$job3 code7_cpg_enrichment/nsc/no_dedup/individual_samples/narrow/flanking/4_peaks_localization.sh | cut -d ' ' -f 4)
job5=$(sbatch --dependency=afterok:$job4 code7_cpg_enrichment/nsc/no_dedup/individual_samples/narrow/flanking/5_analyze_cpg_genes.sh | cut -d ' ' -f 4)
sbatch --dependency=afterok:$job5 code7_cpg_enrichment/nsc/no_dedup/individual_samples/narrow/flanking/6_analyze_coding_sequences.sh

#########################################################################################################################

job1=$(sbatch code7_cpg_enrichment/nsc/dedup/individual_samples/broad/flanking/1_cpg_enrichment_consistent_peaks.sh | cut -d ' ' -f 4)
job2=$(sbatch --dependency=afterok:$job1 code7_cpg_enrichment/nsc/dedup/individual_samples/broad/flanking/2_combine_chunks.sh | cut -d ' ' -f 4)
job3=$(sbatch --dependency=afterok:$job2 code7_cpg_enrichment/nsc/dedup/individual_samples/broad/flanking/3_split_the_list.sh | cut -d ' ' -f 4)
job4=$(sbatch --dependency=afterok:$job3 code7_cpg_enrichment/nsc/dedup/individual_samples/broad/flanking/4_peaks_localization.sh | cut -d ' ' -f 4)
job5=$(sbatch --dependency=afterok:$job4 code7_cpg_enrichment/nsc/dedup/individual_samples/broad/flanking/5_analyze_cpg_genes.sh | cut -d ' ' -f 4)
sbatch --dependency=afterok:$job5 code7_cpg_enrichment/nsc/dedup/individual_samples/broad/flanking/6_analyze_coding_sequences.sh

job1=$(sbatch code7_cpg_enrichment/nsc/dedup/individual_samples/narrow/flanking/1_cpg_enrichment_consistent_peaks.sh | cut -d ' ' -f 4)
job2=$(sbatch --dependency=afterok:$job1 code7_cpg_enrichment/nsc/dedup/individual_samples/narrow/flanking/2_combine_chunks.sh | cut -d ' ' -f 4)
job3=$(sbatch --dependency=afterok:$job2 code7_cpg_enrichment/nsc/dedup/individual_samples/narrow/flanking/3_split_the_list.sh | cut -d ' ' -f 4)
job4=$(sbatch --dependency=afterok:$job3 code7_cpg_enrichment/nsc/dedup/individual_samples/narrow/flanking/4_peaks_localization.sh | cut -d ' ' -f 4)
job5=$(sbatch --dependency=afterok:$job4 code7_cpg_enrichment/nsc/dedup/individual_samples/narrow/flanking/5_analyze_cpg_genes.sh | cut -d ' ' -f 4)
sbatch --dependency=afterok:$job5 code7_cpg_enrichment/nsc/dedup/individual_samples/narrow/flanking/6_analyze_coding_sequences.sh