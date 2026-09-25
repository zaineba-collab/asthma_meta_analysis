# Appendix and supplementary-material plan

Do not duplicate raw GWAS or FinnGen files. Link or cite the compact outputs below.

## Recommended appendix tables

- Full 405-locus table: `results/loci/asthma_fixed_sentinel_merge_loci.txt`.
- Complete 5,254 independent-lead-SNP table: `results/ld_clumping/independent_lead_snps.tsv`.
- Full MR five-method comparison: `results/mr/finngen/summary/finngen_mr_five_method_comparison.tsv`.
- Full MR sensitivity, single-SNP, and leave-one-out tables under `results/mr/finngen/sensitivity/`.
- Complete LDSC h2 QC table: `results/ldsc/qc/ldsc_h2_summary.tsv` and readiness table.
- LDSC technical/environment audits only if examiner reproducibility detail is required.

## Recommended appendix figures

- LDSC h2 QC plot copied here as `Appendix_figure_LDSC_h2_QC.*`.
- Outcome-specific MR scatter, funnel, leave-one-out, and single-SNP forest plots under `results/mr/finngen/plots/`.
- MR sensitivity overview from `results/mr/finngen/summary/figures/`.

## Annotation hold

- Do not include the active full VEP output, gene table, consequence table, or annotation plot as scientific results.
- The active VEP input used PLINK coordinates later validated as GRCh37/hg19, while VEP used GRCh38.
- A corrected VEP run is required in a separately authorised analysis stage before annotation can enter the dissertation.

## Technical records normally excluded from the main chapter

- Coordinate-repair debugging details, P=0 underflow audit, HDL/LDL estimator audit, installation logs, and checksums.
- Retain these in the project archive; include only concise methodological caveats when needed.
