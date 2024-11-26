rule extended_analysis:
    input:
        peaks = expand(join(config["peaks_outdir"], "{condition}_{replicate}_peaks.narrowPeak"), 
                      condition=CONDITIONS, replicate=["r1", "r2", "r3"]),
        bams = expand(join(config["aligned_outdir"], "{sample}_sorted.bam"), 
                     sample=SAMPLES)
    output:
        peak_dist = join(config["qc_outdir"], "peak_distribution.pdf"),
        go_enrich = expand(join(config["func_annot_outdir"], "{condition}_GO_enrichment.pdf"),
                          condition=CONDITIONS),
        diff_bind = join(config["diff_meth_outdir"], "PCA_plot.pdf")
    conda:
        "envs/r_analysis.yaml"
    script:
        "scripts/extended_analysis.R" 