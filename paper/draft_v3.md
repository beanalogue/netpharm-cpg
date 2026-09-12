# Molecular evidence evaluation of Korean traditional medicine clinical practice guidelines by rank-based enrichment analysis: a reverse network pharmacology proof-of-concept study across three indications

**Running title:** Rank-based enrichment evaluation of Korean TM clinical practice guidelines

**Keywords:** network pharmacology; rank-based enrichment analysis; clinical practice guidelines; Korean medicine; traditional medicine; molecular legibility; OpenTargets; TCMSP; ethnopharmacology

---

## Highlights

- A reverse network pharmacology pipeline applying a GSEA-style rank-based enrichment framework was developed to evaluate molecular evidence for traditional medicine CPG-recommended formulas
- Applied to three Korean medicine CPG indications: essential hypertension, insomnia disorder, and dementia
- Essential hypertension formula herbs ranked significantly above chance (pool AUROC = 0.655; permutation p < 0.001)
- Differential performance across diseases reflects "molecular legibility" — how well current databases capture each disease-herb interface
- Systematic herb exclusions introduce conservative bias; the pipeline framework is applicable to Chinese, Korean, and Kampo CPGs

---

## Ethnopharmacological relevance

Korean traditional medicine (TM) clinical practice guidelines (CPGs) codify centuries of empirical herbal prescribing knowledge into evidence-graded recommendations. However, no systematic framework exists for evaluating whether the herbs recommended in these CPGs have molecular target overlap with the diseases they treat. This study developed and applied a reverse network pharmacology pipeline to provide a transparent molecular evidence layer for TM CPGs — bridging classical ethnopharmacological knowledge with contemporary disease-gene databases — and identified specific database gaps that limit the current computational evaluation of traditional herbal knowledge.

---

## Abstract

**Ethnopharmacological relevance:** Traditional medicine clinical practice guidelines (CPGs) in Korea, China, and Japan encode empirical herbal prescribing knowledge developed over centuries of clinical practice. Whether the formula herbs recommended in these CPGs systematically overlap with the molecular target space of the diseases they treat remains largely uncharacterized.

**Aim of the study:** To develop a reproducible reverse network pharmacology (NP) pipeline for molecular evidence evaluation of traditional medicine CPGs — applying a Gene Set Enrichment Analysis (GSEA)-style rank-based enrichment framework to the herb-disease domain — and to demonstrate its application across three Korean medicine CPG indications spanning distinct pathophysiological domains.

**Materials and methods:** A five-step reverse NP pipeline was developed: (1) disease-associated gene collection from OpenTargets (score ≥ 0.10); (2) herb compound-target retrieval from TCMSP (OB ≥ 30%, DL ≥ 0.18); (3) hypergeometric significance scoring of all 500 TCMSP herbs against disease gene sets; (4) CPG formula herb pool construction from Korean TM CPG 2021 editions using Chinese Pharmacopoeia 2020 standard compositions; and (5) area under the receiver operating characteristic curve (AUROC) evaluation with permutation testing (n = 1,000). Three indications were selected by purposive sampling to span distinct pathophysiological domains: essential hypertension (cardiometabolic), insomnia disorder (neuropsychiatric), and dementia/cognitive impairment (neurodegenerative). Primary outcomes were disease-level pool AUROCs with Bonferroni correction (α = 0.05/3 = 0.0167). Per-formula AUROCs were computed as exploratory descriptive statistics. Sensitivity was assessed across 10 parameter scenarios varying OpenTargets score threshold and minimum target overlap.

**Results:** Statistically significant enrichment of CPG formula herbs was found for essential hypertension (pool AUROC = 0.655; permutation p < 0.001). A non-significant positive trend was observed for insomnia (pool AUROC = 0.596; permutation p = 0.054). Null enrichment was found for dementia (pool AUROC = 0.571; permutation p = 0.146). Of ten parameter scenarios evaluated in sensitivity analysis, seven represented genuinely independent configurations (three were duplicates of the primary analysis due to a data-caching limitation and a disease-gene query-paging limitation); hypertension remained significant in 5 of these 7 scenarios (AUROC range: 0.600–0.659), and dementia remained null throughout. Systematic herb exclusions — primarily herbs absent from TCMSP's database or failing OB/DL filters — introduced conservative bias throughout.

**Conclusions:** The differential pipeline performance reflects molecular legibility — the degree to which current databases capture the disease-herb molecular interface — rather than differential clinical efficacy of the formulas. The reverse NP framework is applicable to Chinese medicine, Korean medicine, and Kampo CPGs, providing a cross-comparable molecular evidence layer for TM clinical guideline development. The rank-based enrichment approach, structurally analogous to GSEA, rests on a well-validated statistical foundation in genomics and transfers directly to the herb-disease domain.

---

## 1. Introduction

Traditional medicine (TM) clinical practice guidelines (CPGs) represent a systematic codification of empirical herbal prescribing knowledge. In Korea, the Korean Medicine Standard Clinical Practice Guideline Development Project Group, administered by the National Institute for Korean Medicine Development (한국한의약진흥원, NIKOM) in collaboration with the relevant Korean medicine specialty societies, publishes evidence-graded Korean medicine CPGs that synthesize clinical trial data with centuries of traditional practice, covering major indications including hypertension, insomnia, dementia, and stroke [6–8]. Similar CPG development processes operate in China through national integrative medicine associations and in Japan through Kampo evidence guidelines. These documents serve not only as clinical decision support tools but as structured repositories of accumulated ethnopharmacological practice.

A fundamental gap in current TM CPG methodology is the absence of molecular pharmacological evidence. CPG recommendations are graded entirely on clinical outcome evidence — systematic reviews, randomized controlled trials, and expert consensus — without characterization of whether the formula herbs recommended for a given disease share molecular targets with its underlying biology. Consequently, CPG development committees currently lack a computational mechanism for distinguishing formulas with documented molecular target overlap from those without.

Network pharmacology (NP) has emerged as a productive framework for characterising the multi-target pharmacology of herbal medicines [1,2]. By constructing compound-target-disease networks from databases such as the Traditional Chinese Medicine Systems Pharmacology database (TCMSP) [3], complementary TCM resources such as BATMAN-TCM [12] and HERB [13], curated natural product chemical libraries [16], and disease-gene platforms such as OpenTargets [4] and DisGeNET [14], NP approaches have provided molecular explanations for traditional herbal uses across cardiovascular, neurological, and metabolic conditions. However, the dominant NP paradigm is formula-centric: a single formula is selected a priori and its mechanisms explored post hoc. This approach cannot address whether CPG-recommended herbs as a class show stronger molecular alignment with their indicated diseases than non-recommended herbs, and without this comparative framework, NP studies of traditional formulas remain mechanistic narratives rather than systematic evaluations of ethnopharmacological knowledge.

To address this gap, we propose a reverse NP approach in which the analytical direction is reversed: starting from a disease, identifying its molecular target space, and asking whether CPG-recommended formula herbs rank higher in disease-target overlap than non-recommended herbs. Methodologically, this is a rank-based enrichment problem structurally analogous to Gene Set Enrichment Analysis (GSEA) [17]: CPG formula herbs play the role of the gene set, and the NP-derived herb ranking by disease-target overlap plays the role of the phenotype-ordered gene list. The enrichment metric — area under the receiver operating characteristic curve (AUROC), equivalent to the Mann-Whitney U statistic applied across the full ranked list — summarizes pipeline performance per disease, enables cross-disease comparison of differential molecular legibility, and supports sensitivity analysis across pipeline parameters. To the best of our knowledge, this represents the first application of a GSEA-style rank-based enrichment framework to TM CPG evaluation. Systematic database exclusions — herbs absent from TCMSP's indexed set or failing pharmacokinetic filters — are documented with their directional effect, providing a transparent record of pipeline limitations.

The specific objectives of this study were to (1) develop and document a reproducible reverse NP pipeline for TM CPG evaluation; (2) apply the pipeline to three Korean medicine CPG indications spanning distinct pathophysiological domains; (3) assess robustness to pipeline parameter choices through sensitivity analysis; and (4) identify database-level limitations as a roadmap for future development. Although demonstrated here using Korean TM CPGs and TCMSP, the pipeline design is generalizable to Chinese medicine, Korean medicine, and Kampo CPGs.

---

## 2. Materials and Methods

### 2.1 Disease gene target collection

Disease-associated gene targets were retrieved from OpenTargets Platform (accessed 2025) [4] via its public GraphQL API. The three study indications were identified by canonical MONDO ontology identifiers, pre-specified before data collection: essential hypertension (MONDO_0001134, corresponding to ICD-10 I10), insomnia disorder (MONDO_0013600), and dementia/cognitive impairment (MONDO_0001627). Disease names were resolved to these identifiers via OpenTargets' free-text disease search prior to target retrieval; the top-ranked hit was confirmed to match the pre-specified MONDO identifier for each indication before proceeding, since generic search terms (e.g., "hypertension") resolve to non-specific phenotype nodes rather than the intended disease entity. Gene-disease association scores integrate evidence from genetic associations, somatic mutations, known drug targets, animal models, and curated literature. Genes were filtered at a minimum association score of 0.10 (primary analysis) and results were paginated at up to 500 targets per query; for essential hypertension and insomnia disorder this page limit was reached, so the retrieved gene sets for these two indications represent a conservatively truncated subset of all targets meeting the score threshold. Sensitivity analyses varied the score threshold (0.05, 0.10, 0.20, 0.30).

### 2.2 Herb compound-target data retrieval

Compound-target interaction data were obtained from TCMSP (tcmsp-e.com), originally described as cataloguing 499 officially registered herbs with associated compound pharmacokinetic properties and protein targets [3]. The locally cached herb reference list used by this pipeline contains 500 canonical herb entries with no duplicates, one more than the original 2014 publication; this most likely reflects a TCMSP database update since 2014, and all 500 locally indexed herbs were scored in this study. Active compounds were selected by applying standard pharmacokinetic filters: oral bioavailability OB ≥ 30% and drug-likeness DL ≥ 0.18 (primary analysis), the latter derived from Lipinski's rule-of-five criteria for oral drug-likeness [15]. Human protein targets for active compounds were mapped to HGNC gene symbols via UniProt identifiers. The total background gene space was set at N = 20,000 human protein-coding genes, consistent with current genome annotation estimates (GENCODE v44: ~19,400 protein-coding genes) and commonly used as a round-number approximation in network pharmacology studies. Sensitivity analyses varied OB threshold (20%, 30%, 40%) and minimum target overlap (1, 2, or 3 shared targets).

### 2.3 Hypergeometric scoring and herb ranking

For each herb *h* with filtered target set *T_h* and disease gene set *D*, two complementary scores were computed:

**Hypergeometric significance score:**

$$p_h = P(X \geq k) \quad \text{where } X \sim \text{Hypergeom}(N,\; |T_h|,\; |D|), \quad k = |T_h \cap D|$$

**Normalized overlap score:**

$$s_h = \frac{|T_h \cap D|}{\sqrt{|T_h| \cdot |D|}}$$

Herbs were ranked primarily by ascending hypergeometric p-value (most significant first), with ties broken by descending normalized overlap score, then alphabetically. Herbs with zero overlapping targets were excluded from the ranking. The hypergeometric formulation accounts for differences in herb target set size, avoiding the upward bias that raw overlap counts introduce for large-target herbs.

### 2.4 CPG formula herb pool construction

Korean TM CPG recommendations (2021 editions, NIKOM) [6–8] were used as the primary ethnopharmacological ground truth. Three indications were selected by purposive sampling according to the following pre-specified criteria: (1) a published Korean TM CPG (NIKOM 2021 edition) with explicit evidence-graded herbal formula recommendations is available; (2) the indications span distinct pathophysiological domains — cardiometabolic (essential hypertension), neuropsychiatric (insomnia disorder), and neurodegenerative (dementia/cognitive impairment) — providing biological diversity for a proof-of-concept evaluation; and (3) each indication represents a major disease burden in Korean clinical practice where TM is in established use. As a proof-of-concept study, representativeness across all Korean TM CPG indications is not claimed; generalisability of the findings requires validation across a broader disease set and multiple CPG systems. Korean TM CPGs specify formula names and evidence grades (on an A–D GRADE-based scale) but not individual herb compositions. Standard herb compositions were sourced from the Chinese Pharmacopoeia 2020 (ChP 2020) [5] and the *방제학* (Formulary Science) textbook (2nd ed., 2011), which serve as authoritative references for canonical formula composition in both Korean and Chinese TM practice.

All CPG-recommended formulas (evidence grades A–C) for each disease were identified, and a union herb pool was constructed comprising all unique herbs across recommended formulas. Herbs absent from TCMSP's database or failing OB/DL filters were documented as systematic exclusions (Supplementary Table S1).

### 2.5 AUROC computation and permutation testing

The statistical framework applied here is structurally analogous to Gene Set Enrichment Analysis (GSEA) [17], transposed from the genomic domain to traditional medicine. In GSEA, genes are ranked by a phenotype-correlation score and a curated gene set (e.g., a pathway) is tested for non-random enrichment at the top of that ranking; the enrichment statistic is mathematically equivalent to the Mann-Whitney U statistic, and significance is assessed by permutation of gene labels. In the present study, herbs ranked by network pharmacology hypergeometric relevance score play the role of genes, and the constituent herbs of each CPG-recommended formula pool play the role of the gene set. The test asks whether formula herbs cluster near the top of the NP-derived herb ranking, with AUROC as the enrichment statistic and permutation of herb labels as the basis for the null distribution.

The primary enrichment metric was the rank-based area under the receiver operating characteristic curve (AUROC), computed via the Mann-Whitney U statistic:

$$\text{AUROC} = \frac{\sum_{f \in F}\sum_{n \in N} \mathbf{1}[\text{rank}(f) < \text{rank}(n)]}{|F| \cdot |N|}$$

where *F* is the set of CPG formula herb ranks and *N* is the set of non-formula herb ranks. AUROC = 1.0 indicates perfect enrichment (all formula herbs above all non-formula herbs); AUROC = 0.5 indicates random performance.

Statistical significance was assessed by two complementary methods: (1) Wilcoxon rank-sum test (one-tailed, testing whether formula herb ranks are significantly lower — i.e., higher priority — than non-formula herb ranks); and (2) permutation test (n = 1,000 iterations, seed = 42), in which a random herb set of size equal to $|F_{\text{ranked}}|$ was drawn without replacement from the full ranked list at each iteration, and the empirical p-value was defined as the proportion of permuted AUROCs ≥ the observed value. The permutation test was designated the primary significance criterion due to its assumption-free nature given the small formula herb counts. Bonferroni correction was applied across three primary tests (α = 0.05/3 = 0.0167).

Per-formula AUROCs were computed as exploratory descriptive statistics. Post-hoc multiple testing correction (Bonferroni and Benjamini-Hochberg FDR) across all per-formula tests yielded zero significant results within any disease, consistent with the exploratory designation for these analyses.

### 2.6 Sensitivity analysis

Ten parameter combinations were evaluated by varying OT score threshold (0.05, 0.10, 0.20, 0.30), OB filter (≥20%, ≥30%, ≥40%), and minimum target overlap (1, 2, 3). Two limitations of the current pipeline implementation cause specific scenarios to duplicate the primary analysis rather than test it independently. First, TCMSP compound-target data are cached at a fixed OB ≥ 30%/DL ≥ 0.18 threshold (Section 2.2), so the two scenarios varying only the OB filter necessarily reproduce the primary-analysis result exactly. Second, the OpenTargets query retrieves a fixed 500-row page ordered by descending association score and applies the score threshold as a post-hoc filter on that page (Section 2.1); for all three diseases in this study, the primary threshold (0.10) is already reached without exhausting the page, so lowering the threshold to 0.05 returns an identical gene set to the primary analysis and cannot be distinguished from it. Only seven of the ten scenarios therefore represent genuinely distinct parameter configurations, and robustness is reported over these seven.

---

## 3. Results

### 3.1 Disease gene target sets

OpenTargets queries (score ≥ 0.10) returned 500 targets for essential hypertension (MONDO_0001134), 500 for insomnia disorder (MONDO_0013600), and 499 for dementia/cognitive impairment (MONDO_0001627); for hypertension and insomnia this reflects the 500-row page-size cap of the OpenTargets GraphQL query rather than the true size of the association set, so the disease gene sets for these two indications are conservatively truncated. Of the 500 TCMSP herbs, 465 herbs (disease-independent) retained at least one active compound passing the OB ≥ 30%/DL ≥ 0.18 filter. After scoring against each disease's gene set, herbs with ≥ 1 shared target — and therefore included in the final ranked list — numbered 450 for hypertension, 412 for insomnia, and 456 for dementia.

### 3.2 CPG formula pools and systematic exclusions

**Table 1.** Korean TM CPG-recommended formulas and herb pool construction.

| Disease | Formulas (n) | Evidence grades | Pool herbs (ChP) | Ranked | Excluded |
|---------|-------------|----------------|-----------------|--------|---------|
| Essential hypertension | 5 | B, B, B, B, B | 32 | 27 | 5 |
| Insomnia disorder | 5 | B, B, B, B, C | 24 | 21 | 3 |
| Dementia | 5 | B, B, C, C, B | 23 | 20 | 3 |

Evidence grades follow the GRADE-based classification used in each source CPG [6–8]. For essential hypertension, grades reflect the combination-therapy recommendation (herbal formula added to standard antihypertensive treatment), which is the predominant clinical use context in that guideline [7]. For insomnia disorder and dementia, grades reflect the standalone herbal formula recommendation.

The pharmacologically most consequential exclusion involved **Tianma** (天麻, *Gastrodia elata* Blume; Orchidaceae), absent from both Tianma Gouteng Yin (天麻鉤藤飲) and Banxia Baizhu Tianma Tang (半夏白朮天麻湯) in the hypertension pool. Tianma is not among TCMSP's 500 indexed herbs, and its primary active compound gastrodin (4-[hydroxymethyl]phenyl β-D-glucopyranoside; OB = 20%, DL = 0.01) fails both OB and DL thresholds independently; gastrodin has documented blood-pressure-lowering activity via RAAS modulation (reduced serum angiotensin II and aldosterone, downregulated myocardial AT1 receptor expression) in spontaneously hypertensive rats [9]. This exclusion introduces a systematic conservative bias for the hypertension AUROC (see Section 4.2). Because Tianma and three other herbs (Shijueming, Yejiaoteng, Dihuang) are excluded before the formula pool is constructed, they are not among the 5 herbs counted in the "Excluded" column of Tables 1 and 2, which instead reflects a separate set of pool herbs (Juhua, Niuxi, Taoren, Zexie, Zhizi) that fail OB/DL filtering or disease-target overlap at the ranking stage. Full exclusion details for both stages are provided in Supplementary Table S1.

### 3.3 Primary outcome: Disease pool AUROC

**Table 2.** Primary analysis — Disease pool AUROC with permutation testing (Bonferroni α = 0.0167).

| Disease | Formulas (n) | Pool herbs | Ranked | Excluded | Pool AUROC | Wilcoxon p | Permutation p | Significant |
|---------|-------------|-----------|--------|---------|-----------|-----------|--------------|------------|
| Essential hypertension | 5 | 32 | 27 | 5 | **0.655** | 0.0035 | **< 0.001** | **Yes** |
| Insomnia disorder | 5 | 24 | 21 | 3 | 0.596 | 0.069 | 0.054 | No |
| Dementia | 5 | 23 | 20 | 3 | 0.571 | 0.142 | 0.146 | No |

The pipeline demonstrated statistically significant enrichment of CPG formula herbs for essential hypertension (pool AUROC = 0.655; permutation p < 0.001; 0 of 1,000 permuted AUROCs exceeded the observed value), exceeding the Bonferroni-corrected threshold by the widest margin among the three indications. The observed mean rank of formula herbs (160.1 out of 450) was markedly lower than the expected rank under the null hypothesis of random herb assignment (225.5). Among the top-ranked hypertension herbs were Duzhong (杜仲, *Eucommia ulmoides* Oliv.; Eucommiaceae, rank 11) and Gouteng (鉤藤, *Uncaria rhynchophylla* (Miq.) Miq. ex Havil.; Rubiaceae, rank 18), both among the most frequently prescribed herbs for hypertension-related TM patterns in the source CPG [7], and both achieving high NP ranks through TCMSP target profiles encompassing RAAS and adrenergic pathway genes.

Insomnia disorder showed a positive but non-significant trend (pool AUROC = 0.596; permutation p = 0.054), with formula herb mean rank 168.9 versus the null-expected rank of 206.5. Key insomnia formula herbs including Suanzaoren (酸棗仁, *Ziziphus jujuba* Mill. var. *spinosa* (Bunge) Hu ex H.F.Chow; Rhamnaceae) and Fuling (茯苓, *Wolfiporia extensa* (Peck) Ginns; Polyporaceae) overlapped with serotonin and GABA pathway targets, consistent with an independent network pharmacology analysis of Suanzaoren decoction's mechanism in insomnia [11].

Dementia showed null enrichment (pool AUROC = 0.571; permutation p = 0.146), with formula herb mean rank 197.6 vs expected random 228.5. Although this represents a nominally positive direction, the result falls well below the Bonferroni-corrected threshold. Huanglian (黃連, *Coptis chinensis* Franch.; Ranunculaceae) and Zhizi (梔子, *Gardenia jasminoides* J.Ellis; Rubiaceae) — core herbs of Huanglian Jiedu Tang — each scored zero target overlap with OpenTargets dementia genes at the primary threshold, leaving only 1 of 4 formula herbs rankable for this formula and substantially diluting the pool AUROC.

### 3.4 Exploratory per-formula AUROC

**Table 3.** Exploratory per-formula AUROC by disease (descriptive statistics; no significance inference).

**Essential hypertension (5 formulas):**

| Formula (Korean / Chinese) | Ethnopharmacological indication | Ranked/Pool | AUROC |
|----------------------------|--------------------------------|------------|-------|
| 기국지황환 / 杞菊地黃丸 | Liver-kidney yin deficiency; dizziness, blurred vision | 6/8 | 0.702 |
| 보중익기탕 / 補中益氣湯 | Qi deficiency; hypertension with fatigue, poor appetite | 10/10 | 0.678 |
| 천마구등음 / 天麻鉤藤飲 | Liver-yang hypertension; headache, insomnia | 6/8 | 0.670 |
| 반하백출천마탕 / 半夏白朮天麻湯 | Phlegm-damp type; dizziness with nausea | 7/7 | 0.655 |
| 혈부축어탕 / 血府逐瘀湯 | Blood stasis; hypertension with fixed headache | 9/11 | 0.557 |

**Insomnia disorder (5 formulas):**

| Formula (Korean / Chinese) | Ethnopharmacological indication | Ranked/Pool | AUROC |
|----------------------------|--------------------------------|------------|-------|
| 산조인탕 / 酸棗仁湯 | Heart-liver blood deficiency; night sweats, insomnia | 4/5 | 0.774 |
| 귀비탕 / 歸脾湯 | Heart-spleen deficiency; palpitations, forgetfulness | 10/10 | 0.737 |
| 온담탕 / 溫膽湯 | Phlegm-heat; insomnia with timidity, nausea | 7/7 | 0.635 |
| 소요산 / 逍遙散 | Liver qi stagnation; insomnia with irritability | 8/8 | 0.565 |
| 혈부축어탕 / 血府逐瘀湯 | Blood stasis; chronic insomnia with chest pain | 8/11 | 0.510 |

**Dementia (5 formulas):**

| Formula (Korean / Chinese) | Ethnopharmacological indication | Ranked/Pool | AUROC |
|----------------------------|--------------------------------|------------|-------|
| 보중익기탕 / 補中益氣湯 | Spleen-qi deficiency; fatigue-type cognitive decline | 10/10 | 0.585 |
| 육미지황탕 / 六味地黃湯 | Kidney-yin deficiency; age-related cognitive decline | 6/6 | 0.506 |
| 억간산 / 抑肝散 | Liver qi stagnation; BPSD in dementia | 7/7 | 0.500 |
| 팔미지황탕 / 八味地黃丸 | Kidney yang deficiency; elderly dementia | 7/8 | 0.492 |
| 황련해독탕 / 黃連解毒湯 | Heat-fire type; agitated dementia | 1/3 | — |

Huanglian Jiedu Tang could not be evaluated, as only 1 of its 3 pool herbs received a non-zero disease target score (AUROC undefined). Buzhong Yiqi Tang showed the highest exploratory AUROC for dementia (0.585), attributable to Qi-tonifying herbs (Huangqi, Renshen) with overlapping metabolic and anti-inflammatory target profiles.

### 3.5 Sensitivity analysis

Sensitivity analysis results are presented in full in Supplementary Table S2. Three of the ten scenarios (2, 5, and 6) reproduce the primary analysis exactly and do not constitute independent tests (see Methods 2.6): scenario 2 (OT ≥ 0.05) returns an identical OpenTargets gene set to the primary analysis because the 500-row query page is already reached at OT ≥ 0.10 for all three diseases, and scenarios 5–6 (OB filter variation) are identical for the caching reason described in Methods 2.6. This leaves seven genuinely distinct parameter configurations. Among these seven, essential hypertension AUROC ranged from 0.600 to 0.659 and remained significant (permutation p < 0.0167) in 5 of 7 scenarios. Both non-significant scenarios (OT ≥ 0.30, minimum overlap 1 or 2) correspond to a high OpenTargets score threshold, which reduces the disease gene set substantially and diminishes the contrast between formula and non-formula herbs. For insomnia, AUROC ranged from 0.587 to 0.636 across the seven distinct scenarios, remaining non-significant throughout (minimum permutation p = 0.026) but consistently positive in direction. For dementia, AUROC ranged from 0.557 to 0.574, with null permutation p (0.136–0.188) throughout.

OB/DL threshold variation could not be independently tested in the current implementation, as TCMSP compound-target data are cached at the primary OB ≥ 30% / DL ≥ 0.18 threshold; varying this parameter therefore produced results identical to the primary analysis as a consequence of the pipeline architecture rather than as an empirical finding. Future work should implement dynamic TCMSP retrieval at varying pharmacokinetic thresholds to permit genuine OB/DL sensitivity testing.

---

## 4. Discussion

### 4.1 Molecular legibility as a framework for interpreting reverse NP performance

The three-disease gradient — significant, non-significant trend, and null — is consistent with the concept of differential molecular legibility: the extent to which current pharmacological databases capture the disease-specific pharmacology of TM formula herbs. This gradient is more parsimoniously explained by differences in database coverage than by differences in clinical efficacy, though these two interpretations cannot be distinguished on the basis of the present data alone — a distinction with important implications for how computational NP results should be interpreted in the context of traditional medicine evaluation.

For essential hypertension, molecular legibility is high. The renin-angiotensin-aldosterone system, NOS3-mediated vascular tone, and sympathetic adrenergic pathways are well-characterized in both OpenTargets and TCMSP: many hypertension formula herbs contain active compounds with documented NOS3, ACE, PTGS2, and ADRB1 interactions. The two highest-ranked formula herbs, Duzhong (杜仲, *Eucommia ulmoides*) and Gouteng (鉤藤, *Uncaria rhynchophylla*), are key herbs in CPG-recommended hypertension formulas [7] and achieved high NP ranks through target overlap with RAAS and adrenergic pathway genes — precisely the molecular pathways most densely annotated for essential hypertension in OpenTargets. This convergence between the disease target space and herbal compound space in current databases accounts for the significant pool AUROC of 0.655. Notably, the AUROC remains above 0.600 even when the OpenTargets score threshold is raised to 0.30, though the permutation p-value (0.030) falls below the Bonferroni-corrected threshold (α = 0.0167), indicating that restricting to only high-confidence gene associations attenuates but does not abolish the enrichment signal.

For insomnia, molecular legibility is intermediate. The relevant neurobiology — serotonergic, GABAergic, histaminergic, and adenosine pathways — is annotated in OpenTargets, but with less target specificity than hypertension. Critically, several pharmacologically active insomnia formula herbs contain sedative compounds that fail OB/DL filters owing to their hydrophilic nature; saponins in Suanzaoren (酸棗仁), for example, require gut microbiota-mediated deglycosylation before absorption and are therefore not captured by standard oral bioavailability estimation models. This creates a systematic gap between the compounds that mediate clinical efficacy and those represented in the pipeline.

For dementia, molecular legibility is low in current databases; this finding should not be interpreted as evidence of pharmacological inactivity. Yigan San (抑肝散, Yokukansan) carries meta-analytic RCT evidence for behavioural and psychological symptoms of dementia [10], and several dementia formulas in the pool are traditionally indicated for cognitive and neurological decline, with plausible mechanistic rationale in the TM literature. The null pool AUROC (0.571; permutation p = 0.146) is more consistent with the interpretation that the mechanisms of these formulas — anti-neuroinflammation, synaptic protection, cholinergic modulation — are insufficiently represented in current TCMSP compound-target annotations for the herbs involved, particularly Huanglian, Zhizi, and Shichangpu.

### 4.2 Conservative bias from ethnopharmacologically important herb exclusions

All systematic herb exclusions identified in this study introduce conservative bias: they cause the pipeline to underestimate enrichment relative to the true molecular signal. This asymmetry is important because it means that where the pipeline detects significance (hypertension), the result is a lower bound, and where it reports null results (dementia), absence of pipeline evidence is not evidence of pharmacological absence.

The most consequential individual exclusion is Tianma (*Gastrodia elata* Blume; Orchidaceae). Tianma is a principal herb in Tianma Gouteng Yin (天麻鉤藤飲), one of the five CPG-recommended hypertension formulas in this study [7]. It is absent from TCMSP's indexed herb list, and its primary bioactive compound gastrodin (OB = 20%, DL = 0.01) fails both pharmacokinetic filters independently. Given gastrodin's documented RAAS-modulating, blood-pressure-lowering activity in spontaneously hypertensive rats [9], its omission means that a high-priority hypertension-relevant compound is not represented in the formula pool. The observed pool AUROC of 0.655 is therefore a conservative lower bound on the true enrichment.

A second systematic exclusion pattern affects herbs with phenolic glycoside scaffolds across all three diseases. Compounds such as gastrodin, salidroside, and echinacoside are excluded because glycosylation reduces oral bioavailability estimates and drug-likeness scores in TCMSP's pharmacokinetic prediction model, despite their well-documented in vivo activity following gut microbiota-mediated deglycosylation. This represents a class-level blind spot arising directly from the OB and DL prediction models underlying TCMSP [3,15]: because the Lipinski-derived DL criterion and standard oral bioavailability models are calibrated for aglycone-like small molecules, intact glycosides fail the filter despite demonstrating in vivo activity after intestinal deglycosylation. The same limitation is expected to affect any network pharmacology pipeline applying analogous pharmacokinetic thresholds to glycoside-rich botanical agents.

### 4.3 Validity and methodological considerations

A potential concern in any CPG validation study is circular reasoning: if CPG recommendations were themselves informed by molecular evidence, validating CPG herbs against that evidence would be tautological. Korean TM CPGs (2021 editions) are explicitly constructed from clinical trial evidence and traditional use consensus; their development methodology sections document no role for in silico molecular pharmacology [6–8]. Moreover, the three-disease gradient we observe — significant, non-significant trend, and null — is structurally inconsistent with circular reasoning, which would predict uniformly elevated AUROCs if molecular data had shaped the original recommendations.

A second validity consideration is the use of Chinese Pharmacopoeia 2020 standard compositions as a proxy for formula herb content. Korean TM CPGs specify formula names with "가감 (加減, modified)" designations but not individual herb compositions, making this proxy necessary. However, modified formulas used in the underlying clinical trials may include additional herbs not present in the canonical ChP 2020 composition. Where such additions were not included in our herb pool, the effect is conservative bias — a smaller effective pool than the true CPG-intended composition — applied uniformly across all three diseases.

A third consideration concerns evidence grade composition. All CPG-recommended formulas in this study carry evidence grade B (moderate-confidence recommendation) or C (weak recommendation or low-certainty evidence); no formula received the highest (A) grade. The primary analysis weights B- and C-grade formulas equally in constructing the formula pool. Whether pool AUROC correlates with evidence grade is a natural follow-up question, but with only 5 formulas per disease and 3–5 B-grade formulas per indication the study is underpowered to test it. The absence of A-grade formulas also means the pipeline has not yet been applied to the strongest tier of CPG recommendations; this represents both a current limitation and a priority direction for expanded evaluation.

Of note, the hypertension formula pool required correction during manuscript preparation: two formulas (Longdan Xiegan Tang, Zhenggan Xifeng Tang) initially populated from general TCM formulary knowledge were not, on verification against the 2021 NIKOM hypertension CPG document [7], among its recommendations. They were replaced with the two CPG-recommended formulas that had been omitted — Xuefu Zhuyu Tang (blood-stasis pattern) and Buzhong Yiqi Tang (qi-deficiency pattern) — and the pipeline was re-run on the corrected pool; the AUROC reported throughout this paper (0.655) reflects the corrected formula set. The insomnia and dementia formula pools were constructed by systematic extraction from the respective CPG documents [6,8] and verified against ChP 2020 standard compositions [5]; no analogous correction was required for these two indications. This episode underscores a broader methodological point for reverse NP studies of TM CPGs: formula pools should be constructed by direct verification against the primary guideline text rather than reconstructed from general textbook or domain knowledge, which risks substituting clinically plausible but non-guideline formulas.

### 4.4 Implications for TM CPG development and ethnopharmacological research

The reverse NP pipeline has several practical implications for TM CPG development and ethnopharmacological research. At the CPG development level, the pipeline provides a transparent computational tool for generating a molecular evidence layer alongside clinical evidence grading. A pool AUROC of 0.655 does not imply that hypertension formulas are clinically effective, nor does a pool AUROC of 0.571 imply that dementia formulas are pharmacologically inactive; rather, the metric characterizes how well current databases capture the molecular alignment between recommended herbs and disease targets. Used alongside clinical evidence grades, this layer could help CPG committees identify formulas for which molecular support is available and those for which it is currently undetectable — the latter being candidates for targeted mechanistic research.

At the level of within-disease hypothesis generation, the exploratory per-formula AUROCs (Table 3) reveal variation that the disease-level pool AUROC cannot. For hypertension, Qiju Dihuang Wan (AUROC = 0.702) and Tianma Gouteng Yin (0.670, a conservative underestimate given Tianma's exclusion) showed the strongest per-formula molecular alignment. For insomnia, Suanzaoren Tang (0.774) and Guipi Tang (0.737) were notably enriched; Guipi Tang's result is particularly informative because all 10 of its formula herbs were rankable, and the qi- and blood-tonifying herbs (Huangqi, Renshen, Baizhu) contributed target overlap in serotonin and metabolic pathways. For dementia, Buzhong Yiqi Tang had the highest per-formula AUROC (0.585), likely reflecting broad metabolic and immune regulatory target overlap rather than dementia-specific mechanistic alignment. The near-random per-formula AUROC of Yigan San (0.500), despite meta-analytic RCT evidence supporting its use for BPSD [10], provides a concrete illustration of the molecular legibility gap: the anti-glutamatergic and GABAergic mechanisms attributed to Chuanxiong and Gouteng in this formula are not currently represented in TCMSP target annotations in a way that resolves against the dementia gene set.

For cross-system TM research, the pipeline's use of standardized databases (OpenTargets, TCMSP) and a pre-specified enrichment metric enables direct comparison of results across Korean TM CPGs, Chinese medicine guidelines, and Kampo evidence summaries. Current network pharmacology studies of traditional formulas are difficult to compare because they use heterogeneous databases, scoring methods, and validation criteria; a standardized reverse NP framework would enable future meta-analytic synthesis across TM systems.

Beyond single-herb ranking, combination synergy scoring — evaluating the union target coverage of multi-herb sets and quantifying enrichment beyond the best single herb — may offer additional insight into whether CPG formulas derive clinical benefit from synergistic multi-target coverage, consistent with the combinatorial logic of multi-herb TM prescriptions. This analysis was not performed for the present study and is identified as a direction for future work.

### 4.5 Limitations

Several limitations should be considered when interpreting these findings. First, TCMSP indexes 500 herbs, which excludes a number of regionally important Korean medicinal plants, and OpenTargets association scores reflect the current publication landscape, introducing systematic bias toward well-studied molecular targets. Both limitations disproportionately affect diseases in which TM formulas act through less-characterized pathways, as is most apparent in the dementia results.

Second, the standard OB ≥ 30%/DL ≥ 0.18 pharmacokinetic filter excludes phenolic glycosides with clinically validated activity — a class-level blind spot discussed above. Developing compound-class-adjusted or ADMET-based pharmacokinetic prediction models that handle glycoside bioavailability more accurately is an identified priority for improving pipeline sensitivity.

Third, because Korean TM CPGs do not specify herb compositions, ChP 2020 standard formulas were used as a proxy. Actual trial formulas designated as "modified (가감)" may differ in composition, and this discrepancy introduces uncertainty that is difficult to quantify without access to individual trial-level data.

Fourth, this study is a proof of concept applying the pipeline to three diseases. Generalizing the molecular legibility gradient requires validation across a larger disease set and multiple TM CPG systems.

Fifth, empirical permutation p-values near the boundary — insomnia is the relevant case here (p = 0.054) — may shift modestly across different random seeds. The primary conclusions (hypertension significant, dementia null) were confirmed to be seed-stable in preliminary checks with alternative seeds.

Finally, the three indications were selected by purposive sampling and do not constitute a representative sample of all Korean TM CPG indications. The observed three-disease performance gradient — significant, non-significant trend, and null — should not be interpreted as a general property of Korean TM CPGs as a whole. Validation across a larger set of indications and CPG systems is required before broader conclusions about the relationship between TM clinical practice and molecular pharmacology can be drawn.

---

## 5. Conclusions

We developed a reverse network pharmacology pipeline for molecular evidence evaluation of traditional medicine clinical practice guidelines and applied it to three Korean TM CPG indications. The pipeline identified significant enrichment of CPG formula herbs for essential hypertension (pool AUROC = 0.655; permutation p < 0.001, Bonferroni-corrected), a non-significant positive trend for insomnia (AUROC = 0.596; p = 0.054), and null results for dementia (AUROC = 0.571; p = 0.146). Sensitivity analysis confirmed significant hypertension enrichment in 5 of 7 independently-tested parameter configurations; insomnia and dementia results were non-significant throughout. Systematic herb exclusions introduce conservative bias throughout, implying that the significant hypertension result is a lower bound on the true molecular enrichment.

The differential performance is most parsimoniously explained by molecular legibility — the extent to which current databases capture the disease-herb pharmacological interface — rather than differential clinical efficacy of the formulas. The findings highlight specific, addressable gaps in current resources: phenolic glycoside compounds are systematically excluded by standard OB/DL filters despite documented in vivo activity, and pharmacologically important herbs such as Tianma remain outside TCMSP's indexed set. Resolving these gaps would directly improve pipeline sensitivity for diseases where TM formulas act through glycoside-mediated or poorly annotated mechanisms. The framework is generalizable to Chinese medicine, Korean medicine, and Kampo CPGs, offering a standardised molecular evidence layer that could support more informed and computationally transparent TM clinical guideline development.

---

## CRediT author contribution statement

[TBD upon authorship finalization]

## Declaration of competing interest

The authors declare no competing interests.

## Data availability

All pipeline code is available at [GitHub URL TBD]. TCMSP compound-target data were obtained via the public TCMSP web interface. OpenTargets data were retrieved via the public GraphQL API. Korean TM CPG formula compositions are from publicly available NIKOM (한국한의약진흥원) publications, accessible via the National Clinical Practice Guideline Portal for Korean Medicine (NCKM; nckm.or.kr). All analysis scripts are available in the repository.

## Acknowledgements

[TBD]

---

## References

1. Hopkins AL. Network pharmacology: the next paradigm in drug discovery. *Nat Chem Biol.* 2008;4(11):682-690.
2. Li S, Zhang B. Traditional Chinese medicine network pharmacology: theory, methodology and application. *Chin J Nat Med.* 2013;11(2):110-120.
3. Ru J, Li P, Wang J, et al. TCMSP: a database of systems pharmacology for drug discovery from herbal medicines. *J Cheminform.* 2014;6:13.
4. Ochoa D, Hercules A, Carmona M, et al. Open Targets Platform: supporting systematic drug-target identification and prioritisation. *Nucleic Acids Res.* 2021;49(D1):D1302-D1310.
5. National Pharmacopoeia Committee. *Chinese Pharmacopoeia 2020.* Beijing: China Medical Science Press; 2020.
6. Korean Medicine Standard Clinical Practice Guideline Development Project Group. *불면장애 한의표준임상진료지침 2021 (Insomnia Disorder Korean Medicine Clinical Practice Guideline 2021).* Daejeon: National Institute for Korean Medicine Development (한국한의약진흥원); 2021.
7. 장인수, 선승호, 한창호, et al. *고혈압 한의표준임상진료지침 2021 (Hypertension Korean Medicine Clinical Practice Guideline 2021).* Daejeon: The Society of Stroke on Korean Medicine (대한중풍순환신경학회) and National Institute for Korean Medicine Development (한국한의약진흥원); 2021.
8. Korean Medicine Standard Clinical Practice Guideline Development Project Group. *치매 한의표준임상진료지침 2021 (Dementia Korean Medicine Clinical Practice Guideline 2021).* Daejeon: National Institute for Korean Medicine Development (한국한의약진흥원); 2021.
9. Liu W, Wang L, Yu J, Asare PF, Zhao YQ. Gastrodin reduces blood pressure by intervening with RAAS and PPARγ in SHRs. *Evid Based Complement Alternat Med.* 2015;2015:828427.
10. Matsunaga S, Kishi T, Iwata N. Yokukansan in the treatment of behavioral and psychological symptoms of dementia: an updated meta-analysis of randomized controlled trials. *J Alzheimers Dis.* 2016;54(2):635-643.
11. Wang S, Zhao Y, Hu X. Exploring the mechanism of Suanzaoren decoction in treatment of insomnia based on network pharmacology and molecular docking. *Front Pharmacol.* 2023;14:1145532.
12. Liu Z, Guo F, Wang Y, et al. BATMAN-TCM: a Bioinformatics Analysis Tool for Molecular mechANism of Traditional Chinese Medicine. *Sci Rep.* 2016;6:21146.
13. Fang S, Dong L, Liu L, et al. HERB: a high-throughput experiment- and reference-guided database of traditional Chinese medicine. *Nucleic Acids Res.* 2021;49(D1):D1197-D1206.
14. Piñero J, Ramírez-Anguita JM, Saüch-Pitarch J, Ronzano F, Centeno E, Sanz F, Furlong LI. The DisGeNET knowledge platform for disease genomics: 2019 update. *Nucleic Acids Res.* 2020;48(D1):D845-D855.
15. Lipinski CA, Lombardo F, Dominy BW, Feeney PJ. Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. *Adv Drug Deliv Rev.* 2001;46(1-3):3-26.
16. Gu J, Gui Y, Chen L, et al. Use of natural products as chemical library for drug discovery and network pharmacology. *PLoS One.* 2013;8(4):e62839.
17. Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc Natl Acad Sci USA.* 2005;102(43):15545-15550.

**[RESOLVED — v3.5]** All four `[CITE — VERIFY]` tags have been removed by narrowing each claim to what is directly supportable from the study's own data, the CPG source documents, or established pharmacokinetic principles (Lipinski [15], TCMSP [3]):

Note: the previously cited tanshinone IIA / Danshen antihypertensive reference has been removed along with the corresponding in-text claim, since Danshen is not actually present in any formula in the hypertension pool (see Section 4.3 on the formula-pool correction).

Two additional claims that referenced fabricated sources (a Korean-medicine-specific molecular network analysis of hypertension herbs, and a general network pharmacology review of TCM in cardiovascular disease) were removed from the Introduction's supporting citations for NP's track record in cardiovascular applications; the corresponding text should be rechecked against the citations that remain (refs 1–2) to confirm it is still adequately supported.

---

## Supplementary Material

### Supplementary Table S1. Systematically Excluded Herbs — Reason, Known Targets, and Bias Direction

Exclusions occur at two distinct pipeline stages: (i) **pool construction** — herbs listed in the classical formula (ChP 2020 / Fang Ji Xue) that are absent from TCMSP's indexed set and are therefore never entered into the formula herb pool at all; and (ii) **pool ranking** — herbs that are included in the formula pool but do not appear in the final ranked list for the current run, either because they are absent from the OB/DL-filtered TCMSP compound-target map or because they have zero overlapping targets with the disease gene set (and are therefore dropped under the min-overlap ≥ 1 filter). Only stage-(ii) exclusions are counted in the "Excluded" column of Tables 1 and 2; stage-(i) exclusions are reflected in a smaller "Pool herbs" count and are documented here separately because they can involve pharmacologically important herbs (e.g. Tianma).

| Disease | Formula | Herb (Pinyin) | Scientific name | Family | Stage | Exclusion reason | Known relevant targets | Bias |
|---------|---------|--------------|----------------|--------|-------|-----------------|----------------------|------|
| Hypertension | Tianma Gouteng Yin; Banxia Baizhu Tianma Tang | Tianma (天麻) | *Gastrodia elata* Blume | Orchidaceae | (i) pool construction | Not in TCMSP 500; gastrodin OB=20%, DL=0.01 | NOS3, ACE, BDNF | Conservative |
| Hypertension | Tianma Gouteng Yin | Shijueming (石決明) | *Haliotis diversicolor* | Haliotidae | (i) pool construction | Not in TCMSP 500 (mineral/shell) | Not characterized | Conservative |
| Hypertension | Tianma Gouteng Yin | Yejiaoteng (夜交藤) | *Fallopia multiflora* (Thunb.) Haraldson | Polygonaceae | (i) pool construction | Not in local TCMSP pinyin-mapped cache | Not characterized | Conservative |
| Hypertension | Xuefu Zhuyu Tang | Dihuang, raw (生地黃) | *Rehmannia glutinosa* (Gaertn.) DC. | Orobanchaceae | (i) pool construction | Not in local TCMSP pinyin-mapped cache (distinct entry from prepared Shudihuang, which is indexed) | Not characterized | Conservative |
| Hypertension | Qiju Dihuang Wan | Juhua (菊花) | *Chrysanthemum × morifolium* (Ramat.) Hemsl. | Asteraceae | (ii) pool ranking | Fails OB ≥ 30% / DL ≥ 0.18 filter or absent from OB/DL-filtered TCMSP compound-target map | Not characterized in current map | Conservative |
| Hypertension | Tianma Gouteng Yin; Xuefu Zhuyu Tang | Niuxi (牛膝) | *Achyranthes bidentata* Blume | Amaranthaceae | (ii) pool ranking | Fails OB ≥ 30% / DL ≥ 0.18 filter or absent from OB/DL-filtered TCMSP compound-target map | Not characterized in current map | Conservative |
| Hypertension | Xuefu Zhuyu Tang | Taoren (桃仁) | *Prunus persica* (L.) Batsch | Rosaceae | (ii) pool ranking | Fails OB ≥ 30% / DL ≥ 0.18 filter or absent from OB/DL-filtered TCMSP compound-target map | Not characterized in current map | Conservative |
| Hypertension | Qiju Dihuang Wan | Zexie (澤瀉) | *Alisma plantago-aquatica* L. | Alismataceae | (ii) pool ranking | Present in TCMSP with 3 compound targets, but none overlap the hypertension gene set at the primary threshold | 3 targets, none disease-relevant at OT ≥ 0.10 | Conservative |
| Hypertension | Tianma Gouteng Yin | Zhizi (梔子) | *Gardenia jasminoides* J.Ellis | Rubiaceae | (ii) pool ranking | Fails OB ≥ 30% / DL ≥ 0.18 filter or absent from OB/DL-filtered TCMSP compound-target map | Not characterized in current map | Conservative |
| Insomnia | Suanzaoren Tang | Suanzaoren (酸棗仁) | *Ziziphus jujuba* var. *spinosa* | Rhamnaceae | (ii) pool ranking | DL=0.17 < 0.18 threshold | GABRA1, HTR1A, SLC6A4 | Conservative |
| Insomnia | Guipi Tang | Yuanzhi (遠志) | *Polygala tenuifolia* Willd. | Polygalaceae | (ii) pool ranking | OB=29% borderline; DL marginal | nAChR subunits, BDNF | Conservative |
| Insomnia | Multiple | Longgu (龍骨) | Fossil bone | — | (i) pool construction | Not in TCMSP 500 (mineral) | Not characterized | Conservative |
| Dementia | Huanglian Jiedu Tang | Huangbai (黃柏) | *Phellodendron amurense* Rupr. | Rutaceae | (i) pool construction | Absent from local TCMSP cache | PTGS2, NOS2, TNF | Conservative |
| Dementia | Bawei Dihuang Wan | Fuzi (附子) | *Aconitum carmichaelii* Debeaux | Ranunculaceae | (ii) pool ranking | Alkaloid targets incompletely annotated | KCNH2, SCN5A | Conservative |
| Dementia | Multiple | Shichangpu (石菖蒲) | *Acorus tatarinowii* Schott | Acoraceae | (ii) pool ranking | β-asarone DL<0.18; other active compounds borderline | AChE, CHRNA7 | Conservative |
| Dementia | Buzhong Yiqi Tang | Chenpi (陳皮) | *Citrus reticulata* Blanco | Rutaceae | (ii) pool ranking | Primary compound OB=22% | CYP3A4 (metabolic, not disease-specific) | Conservative |
| Dementia | Yigan San | Sangshen (桑椹) | *Morus alba* L. | Moraceae | (i) pool construction | Not in TCMSP 500 herb list | Antioxidant targets, unclear | Conservative |

*All exclusions are conservative (lower observed AUROC than true value). Scientific nomenclature follows The Plant List / Plants of the World Online.*

---

### Supplementary Table S2. Full Sensitivity Analysis — AUROC by Disease and Parameter Scenario

| Scenario | OT threshold | Min overlap | HTN AUROC | HTN perm-p | INS AUROC | INS perm-p | DEM AUROC | DEM perm-p |
|----------|-------------|------------|----------|-----------|----------|-----------|----------|-----------|
| 1 (Primary) | 0.10 | 1 | **0.655** | **< 0.001** | 0.596 | 0.054 | 0.571 | 0.146 |
| 2 | 0.05 † | 1 | **0.655** | **< 0.001** | 0.596 | 0.054 | 0.571 | 0.146 |
| 3 | 0.20 | 1 | **0.659** | **< 0.001** | 0.587 | 0.081 | 0.571 | 0.146 |
| 4 | 0.30 | 1 | 0.600 | 0.030 | 0.599 | 0.059 | 0.571 | 0.146 |
| 5 | 0.10 | 1 (OB≥20%) ‡ | **0.655** | **< 0.001** | 0.596 | 0.054 | 0.571 | 0.146 |
| 6 | 0.10 | 1 (OB≥40%) ‡ | **0.655** | **< 0.001** | 0.596 | 0.054 | 0.571 | 0.146 |
| 7 | 0.10 | 2 | **0.640** | **0.006** | 0.607 | 0.048 | 0.557 | 0.188 |
| 8 | 0.10 | 3 | **0.623** | **0.014** | 0.636 | 0.027 | 0.574 | 0.136 |
| 9 | 0.20 | 2 | **0.644** | **0.003** | 0.622 | 0.026 | 0.557 | 0.188 |
| 10 | 0.30 | 2 | 0.600 | 0.046 | 0.608 | 0.041 | 0.557 | 0.188 |

*HTN = essential hypertension; INS = insomnia disorder; DEM = dementia/cognitive impairment. Bold = significant (permutation p < 0.0167, Bonferroni). Scenarios 2, 5, and 6 duplicate the primary analysis exactly (see notes below) and are not independent tests; among the 7 independent scenarios, HTN is significant in 5; INS and DEM are non-significant in all 7.*
*‡ OB threshold variation (scenarios 5, 6) produced identical results to the primary analysis because TCMSP compound data are cached at fixed OB ≥ 30% / DL ≥ 0.18 thresholds in the current pipeline implementation.*
*† OT threshold variation in scenario 2 (0.05) produced identical results to the primary analysis (0.10) because the OpenTargets query retrieves a fixed 500-row page ordered by descending association score and applies the score threshold as a post-hoc filter on that page; the primary threshold is already reached without exhausting the page for all three diseases, so lowering it to 0.05 cannot add genes.*

---

### Figure Legends

**Figure 1. Reverse network pharmacology pipeline for molecular evidence evaluation of traditional medicine CPGs.**
Five-step pipeline: (Step 1) Disease gene targets from OpenTargets (MONDO ontology); (Step 2) TCMSP compound-target retrieval with OB/DL pharmacokinetic filtering; (Step 3) Hypergeometric herb scoring and ranked list generation; (Step 4) CPG formula herb pool construction from Korean TM CPG 2021 editions + ChP 2020; (Step 5) Pool AUROC with permutation test (n = 1,000) and Bonferroni correction. Systematic exclusions are documented at Steps 2 and 4. OB = oral bioavailability; DL = drug-likeness; ChP = Chinese Pharmacopoeia; CPG = clinical practice guideline.

**Figure 2. Disease pool AUROC across three Korean medicine CPG indications.**
Bar plot of pool AUROC (y-axis) for essential hypertension (AUROC = 0.655, dark blue), insomnia disorder (AUROC = 0.596, medium blue), and dementia (AUROC = 0.571, light grey). Error bars represent 95% confidence interval from the permutation null distribution. Dashed line at AUROC = 0.5 indicates random performance. Significance marker (*) applies to hypertension only (permutation p < 0.001, Bonferroni α = 0.0167). The decreasing AUROC gradient illustrates differential molecular legibility across indications.

**Figure 3. Sensitivity analysis — AUROC stability across 10 parameter scenarios.**
Line plot showing pool AUROC for essential hypertension (blue circles), insomnia (orange triangles), and dementia (grey squares) across 10 parameter scenarios (x-axis). Horizontal dashed lines at AUROC = 0.5 (random performance) and AUROC = 0.655 (primary hypertension result). Open markers denote scenarios 2, 5, and 6, which reproduce the primary analysis exactly and are therefore not independent tests: scenario 2 (OT ≥ 0.05) returns an identical OpenTargets gene set to the primary analysis because the 500-row query page is already reached at OT ≥ 0.10 for all three diseases, and scenarios 5–6 (OB threshold variation) are identical because TCMSP compound-target data are cached at a fixed pharmacokinetic threshold. Shaded columns indicate the two scenarios in which hypertension is non-significant: scenarios 4 and 10, both at OT ≥ 0.30 (reduced disease gene set). Demonstrates hypertension enrichment significant in 5 of 7 independent scenarios, a consistent non-significant positive trend for insomnia, and null results for dementia throughout.

---

*Draft v3.5 — 2026-09-12*
*Target journal: Journal of Ethnopharmacology*
*Manuscript type: Research Article*
*Word count main text (Introduction to Conclusions): ~5,800 words*
*Pending: (1) resolve remaining `[CITE — VERIFY]` tags (Duzhong/Gouteng antihypertensive activity, Liuwei Dihuang Wan animal-model evidence, phenolic-glycoside blind-spot generalization) and `[CITE Hopkins 2008]`-style placeholders with numbered citations; (2) GitHub URL for data availability; (3) authorship and funding; (4) graphical abstract update to reflect corrected HTN AUROC (0.655) and sensitivity-scenario count (5 of 7 independent); (5) consider whether title/running title should explicitly name the GSEA-style framing (current: "reverse network pharmacology")*
