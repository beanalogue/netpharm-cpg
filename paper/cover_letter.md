**Cover Letter**

---

To the Editor-in-Chief,
*Journal of Ethnopharmacology*

I am writing to submit the manuscript entitled **"Rank-based evaluation of herb–disease molecular target overlap in Korean traditional medicine clinical practice guidelines: a reverse network pharmacology pipeline demonstrated across three indications"** for consideration as a Research Article in the *Journal of Ethnopharmacology*.

**Ethnopharmacological relevance and scope fit**

Korean traditional medicine (TM) clinical practice guidelines (CPGs) codify centuries of empirical herbal prescribing knowledge into evidence-graded recommendations. This manuscript develops and applies a reverse network pharmacology pipeline that asks whether CPG-recommended formula herbs show systematically stronger disease-gene target overlap than non-recommended herbs — providing a transparent molecular evidence layer for TM CPG evaluation. The work sits at the intersection of ethnopharmacological knowledge systems and contemporary computational pharmacology, which I believe aligns directly with the aims and scope of the *Journal of Ethnopharmacology*.

**What the study does and why it matters**

Current network pharmacology studies of traditional herbal formulas are formula-centric: a single formula is selected a priori and its mechanisms explored post hoc. This approach cannot address whether CPG-recommended herbs as a class have stronger molecular target alignment with their indicated diseases than non-recommended herbs. The present study reverses this direction: starting from three Korean TM CPG indications — essential hypertension, insomnia disorder, and dementia/cognitive impairment — it scores all 500 TCMSP-indexed herbs by hypergeometric overlap with disease gene sets from OpenTargets and evaluates whether CPG formula herbs are enriched at the top of this ranking, quantified by the area under the receiver operating characteristic curve (AUROC) with permutation-based significance testing.

The key findings are: (1) significant enrichment for essential hypertension (pool AUROC = 0.655; permutation p < 0.001, Bonferroni-corrected), stable across 5 of 7 independent sensitivity configurations; (2) a non-significant positive trend for insomnia (AUROC = 0.596; p = 0.054); and (3) null results for dementia (AUROC = 0.571; p = 0.146). The differential performance is interpreted as reflecting differential database coverage of the herb–disease pharmacological interface — a pattern consistent with the density of approved drug classes targeting each disease — rather than differential clinical efficacy of the formulas. Systematic herb exclusions are documented with their directional effects, providing a transparent record of current pipeline limitations and a roadmap for future improvement.

**Originality and methodological transparency**

The study makes three contributions: (i) a reproducible, pre-specified enrichment framework enabling cross-disease and cross-CPG-system comparison; (ii) explicit documentation of conservative biases introduced by database coverage gaps, with their directionality; and (iii) identification of the specific database-side limitations — phenolic glycoside exclusion, missing herbs, incomplete dementia compound annotation — that currently limit computational evaluation of traditional medicine knowledge.

**Ethical considerations**

This study is entirely computational and uses publicly available databases (OpenTargets, TCMSP) and published CPG documents. No human participants, animal experiments, or patient data were involved. No ethics committee approval was required.

**Data and code availability**

All pipeline code, analysis scripts, CPG formula data, and the manuscript are publicly available at: https://github.com/beanalogue/netpharm-cpg

I confirm that this manuscript has not been published previously and is not under consideration for publication elsewhere. I have no conflicts of interest to disclose.

Yours sincerely,

**Chan-Young Kwon, KMD, PhD**
Department of Oriental Neuropsychiatry
Dong-eui University College of Korean Medicine
Busan, Republic of Korea
