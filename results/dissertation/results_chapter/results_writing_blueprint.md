# Results chapter writing blueprint

Use the exact values below; keep interpretation factual and reserve biological explanation for Discussion.

## 3.1 Overview of the analytical workflow

**Operation carried out**
- Twelve asthma GWAS were harmonised and meta-analysed, followed by signal refinement, FinnGen MR, and LDSC.

**Design features that need mentioning**
- Fixed-effect GWAMA; parallel distance-based locus definition and European-reference LD clumping.

**Results/numbers that MUST appear**
- 12 input GWAS; 34,316,972 meta-analysed variants; 5,254 LD-clumped independent lead SNPs; 12 FinnGen outcomes.

**Statistical analysis**
- State only the project P ≤ 5 × 10⁻⁹ selection threshold here.

**Recommended table/figure**
- Short prose; Figure 3.3 provides the refinement overview.

**Immediate conclusion**
- The workflow produced validated association, MR, and LDSC result sets.

**Logical link to next analysis**
- Begin with the fixed-effect asthma meta-analysis.

## 3.2 Asthma GWAS meta-analysis

**Operation carried out**
- Fixed-effect association estimates were combined across 12 asthma GWAS.

**Design features that need mentioning**
- Distinguish conventional P < 5 × 10⁻⁸ from project signal selection at P ≤ 5 × 10⁻⁹.

**Results/numbers that MUST appear**
- 34,316,972 variants; 50,764 at P < 5 × 10⁻⁸; 39,933 at P ≤ 5 × 10⁻⁹.
- Top association: rs9273374, chromosome 6, corrected GRCh38 position 32,658,837, beta = -0.153052, SE = 0.004504, P = 5.22 × 10⁻²⁵³. The PLINK clumping coordinate was 32,626,614 in GRCh37.
- lambda_GC = 1.3171.

**Statistical analysis**
- Fixed-effect beta, SE, z and P; QQ genomic inflation statistic.

**Recommended table/figure**
- Table 3.1; Figures 3.1 and 3.2.

**Immediate conclusion**
- The strongest association was the chromosome 6 rs9273374 signal.

**Logical link to next analysis**
- Quantify between-study heterogeneity and refine significant signals.

## 3.3 Between-study heterogeneity and association signals

**Operation carried out**
- I2 and Cochran Q summaries were extracted from GWAMA output.

**Design features that need mentioning**
- High heterogeneity was defined in the existing pipeline as I2 ≥ 0.75.

**Results/numbers that MUST appear**
- 229,395 variants (0.668%) had I2 ≥ 0.75.
- 966,173 variants (2.815%) had Cochran Q P < 0.05.

**Statistical analysis**
- Report I2 threshold and Q P-value criterion without adding new tests.

**Recommended table/figure**
- Table 3.1; no additional heterogeneity figure.

**Immediate conclusion**
- High I2 was observed for 0.668% of analysed variants.

**Logical link to next analysis**
- Project-significant variants were reduced to sentinel-defined loci.

## 3.4 Identification of sentinel variants and associated loci

**Operation carried out**
- Distance-based sentinel selection and overlapping ±1 Mb window merging were applied.

**Design features that need mentioning**
- P < 5 × 10⁻⁹; ±1 Mb sentinel windows.

**Results/numbers that MUST appear**
- 1,563 sentinel SNPs and 405 merged loci.
- Median width 4,674,582 bp; mean width 6,047,095 bp; largest 24,489,380 bp; maximum 16 sentinels in one locus.

**Statistical analysis**
- Association P ranking; descriptive locus-size summaries.

**Recommended table/figure**
- Table 3.2 and optional Table 3.2b; Figure 3.3.

**Immediate conclusion**
- Distance-based refinement grouped 1,563 sentinels into 405 loci.

**Logical link to next analysis**
- A separate rsID-mapping and LD-clumping path defined independent lead SNPs.

## 3.5 LD-based refinement of independent association signals

**Operation carried out**
- Project-significant variants were mapped to rsIDs and clumped using a European reference panel.

**Design features that need mentioning**
- P1=P2=5 × 10⁻⁹, r2=0.1, window=99,999 kb.

**Results/numbers that MUST appear**
- 39,933 evaluated; 37,669 rsID matched; 2,235 unmatched; 94.33% matched; 29 duplicate rows collapsed.
- 5,254 independent lead SNPs.

**Statistical analysis**
- PLINK LD clumping parameters and result counts.

**Recommended table/figure**
- Table 3.3; Figure 3.3.

**Immediate conclusion**
- LD clumping yielded 5,254 independent lead SNPs; do not call these 5,254 loci.

**Logical link to next analysis**
- Lead SNPs were intended for functional annotation and MR instrument preparation.

## 3.6 Functional annotation of independent lead variants

**Operation carried out**
- Ensembl VEP/cache release 116 was run offline on the independent lead SNPs with validated GRCh38 coordinates and resolved genomic REF/ALT alleles.

**Design features that need mentioning**
- The authoritative set remained 5,254 SNPs. GRCh37 PLINK positions were not submitted; unresolved coordinates or alleles were retained in QC and excluded from coordinate-based VEP.

**Results/numbers that MUST appear**
- 5,238 of 5,254 SNPs had validated GRCh38 coordinates; 5,050 were VEP-ready; 5,050 were successfully annotated (96.12% of all lead SNPs).
- The five most common raw most-severe consequence categories were intron_variant: 1892 (37.5%); intron_variant,non_coding_transcript_variant: 1259 (24.9%); intergenic_variant: 921 (18.2%); upstream_gene_variant: 242 (4.8%); downstream_gene_variant: 199 (3.9%).
- Report the top variants and gene assignments exactly as supported in Table 3.4; do not infer biological function.

**Statistical analysis**
- Descriptive counts and percentages; one most-severe VEP consequence per successfully annotated SNP.

**Recommended table/figure**
- Table 3.4 and Figure 3.4.

**Immediate conclusion**
- Corrected annotation describes the genomic consequence context of the coordinate- and allele-resolved independent association signals without functional interpretation.

**Logical link to next analysis**
- The annotated lead variants were subsequently used to characterise the functional context of the independent association signals, while the lead SNP set also formed the basis of downstream MR instrument preparation.

## 3.7.1 Instrument preparation and harmonisation

**Operation carried out**
- GRCh38 coordinates were repaired, FinnGen associations extracted, and alleles harmonised.

**Design features that need mentioning**
- Exposure EAF was unavailable; unresolved palindromic variants were removed.

**Results/numbers that MUST appear**
- 5,238/5,254 coordinates recovered (99.70%); 16 unresolved; 3,047 positions changed.
- R13: 810 matched SNPs; R12: 773–798; final MR-ready: 623–652; retention 80.43–80.60%.

**Statistical analysis**
- Descriptive extraction and harmonisation QC.

**Recommended table/figure**
- Table 3.5.

**Immediate conclusion**
- Coordinate repair enabled consistent local FinnGen matching and yielded outcome-specific MR sets.

**Logical link to next analysis**
- The retained instruments entered primary IVW MR.

## 3.7.2 Primary MR results

**Operation carried out**
- IVW estimates were evaluated for 12 FinnGen outcomes.

**Design features that need mentioning**
- Binary outcomes use OR; continuous outcomes use beta; 12-outcome Bonferroni family.

**Results/numbers that MUST appear**
- 0/12 nominal IVW P < 0.05; 0/12 Bonferroni significant. Report all outcome estimates from Table 3.6.

**Statistical analysis**
- IVW beta, SE, 95% CI, P, Bonferroni P and FDR P.

**Recommended table/figure**
- Table 3.6; Figure 3.5.

**Immediate conclusion**
- No primary IVW estimate reached nominal significance; do not state that this proves absence of an effect.

**Logical link to next analysis**
- Sensitivity diagnostics were assessed next.

## 3.7.3 Sensitivity analyses

**Operation carried out**
- IVW heterogeneity, MR-Egger intercept and leave-one-out diagnostics were summarised.

**Design features that need mentioning**
- Diagnostics were retained without SNP exclusion.

**Results/numbers that MUST appear**
- Heterogeneity: 10/12.
- Egger P < 0.05 and FDR significant: 5/12 (Allergic rhinitis; Bronchiectasis; Eosinophilic disease; AST; Triglycerides (fasting)).
- Leave-one-out sign reversal: 3/12 (Eosinophilic disease; NAFLD; Triglycerides (fasting)).

**Statistical analysis**
- Q-test P, Egger intercept P and BH-FDR P, and sign reversal.

**Recommended table/figure**
- Table 3.7; detailed plots in appendix.

**Immediate conclusion**
- Sensitivity analyses identified widespread heterogeneity and directional horizontal pleiotropy flags in several outcomes.

**Logical link to next analysis**
- Genome-wide genetic correlation was assessed to quantify shared genetic architecture without assigning causality.

## 3.8.1 Univariate LDSC quality control

**Operation carried out**
- Observed-scale SNP h2 was estimated for the LDSC asthma dataset and 12 outcomes.

**Design features that need mentioning**
- Asthma was GCST90302886: 48,623 cases, 290,722 controls, N=339,345; it was not the 12-study GWAMA.

**Results/numbers that MUST appear**
- Asthma h2=0.0423, SE=0.0040, Z=10.575; 12/13 traits had positive h2.
- Eosinophilic disease h2=−0.0004, SE=0.0013, Z=−0.308 and was not eligible for stable rg.

**Statistical analysis**
- Univariate LDSC observed-scale h2 and SE; no liability-scale interpretation.

**Recommended table/figure**
- Brief prose; h2 QC table/figure in appendix.

**Immediate conclusion**
- Eleven outcome correlations were estimable; eosinophilic disease was retained as not estimable.

**Logical link to next analysis**
- Run-status-eligible traits were summarised using the completed rg results.

## 3.8.2 Genetic correlation with asthma

**Operation carried out**
- Existing asthma–outcome rg estimates were consolidated.

**Design features that need mentioning**
- 12 planned hypotheses; 11 estimable; primary Bonferroni threshold=0.05/12.
- The Bonferroni threshold was retained at 0.05/12, reflecting the pre-specified set of 12 planned outcomes, rather than recalculated as 0.05/11 after eosinophilic disease was excluded because of non-positive heritability, to avoid post hoc adjustment of the significance criterion.

**Results/numbers that MUST appear**
- Seven Bonferroni- and FDR-significant outcomes: Allergic rhinitis; BMI IRN; NAFLD; ALT; GGT; HDL cholesterol; Fasting triglycerides.
- Strongest positive: allergic rhinitis rg=0.7114, SE=0.0547, 95% CI 0.6042–0.8186, P=1.3027×10⁻³⁸.
- Strongest negative: HDL cholesterol rg=−0.1329, SE=0.0324, 95% CI −0.1964 to −0.0694, P=4.1846×10⁻⁵.

**Statistical analysis**
- rg, SE, 95% CI, P, 12-test Bonferroni P and 11-estimate BH-FDR P.

**Recommended table/figure**
- Table 3.8; Figure 3.6.

**Immediate conclusion**
- Seven outcomes showed significant genetic correlation with asthma.

**Logical link to next analysis**
- Compare the primary MR and LDSC evidence categories.

## 3.9 Comparison of MR and genetic-correlation findings

**Operation carried out**
- Primary Bonferroni significance patterns were joined across MR and LDSC.

**Design features that need mentioning**
- MR and rg answer distinct questions and effect estimates are not directly compared.

**Results/numbers that MUST appear**
- MR significant 0/12; LDSC significant 7/12; both significant 0; LDSC significant/MR not significant 7; neither 4; LDSC not estimable 1.
- LDSC significant/MR not significant: Allergic rhinitis; BMI IRN; NAFLD; ALT; GGT; HDL cholesterol; Fasting triglycerides.
- Neither: Bronchiectasis; CRP; AST; LDL cholesterol; not estimable: Eosinophilic disease.

**Statistical analysis**
- Descriptive cross-classification using prespecified Bonferroni flags.

**Recommended table/figure**
- Table 3.9; Figure 3.7 as main or supplementary.

**Immediate conclusion**
- Seven outcomes showed significant genetic correlation despite no significant primary MR association.

**Logical link to next analysis**
- End Results and move to a separate Discussion chapter.
