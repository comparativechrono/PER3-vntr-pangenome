====
Association statistics for coding VNTRs in UKB
companion to R. E. Mukamel*, R. E. Handsaker*, M. A. Sherman, A. R. Barton, Y. Zheng, S. A. McCarroll**, P.-R. Loh**, "Protein-coding repeat polymorphisms strongly shape diverse human phenotypes", Science (2021).
*These authors contributed equally to this work
**These authors co-supervised this work
Correspondence to: REM (rmukamel@broadinstitute.org), REH (handsake@broadinstitute.org), SAM (mccaroll@genetics.med.harvard.edu), or P-RL (poruloh@broadinstitute.org)
Date: 2021-11-05
====

This tar archive contains the results from association analyses of coding VNTRs in the UK Biobank detailed in the manuscript Mukamel et al., Science 2021.  The archive contains one file for each of 787 phenotypes tested, with each file listing the statistics of association with 118 protein-coding VNTRs.  Association statistics were computed using linear regression with BOLT-LMM v2.3.6, and analyses included assessment center, genotyping array, sex, age, age squared, and 20 PCs as covariates.  Analysis was restricted to up to N=337,466 unrelated, White British UKB participants.  VNTR genotypes were converted to real-valued dosages by summing the estimated lengths of the two parental alleles and linearly rescaling so the minimum and maximum dosages were 0 and 2 respectively.  Reported betas roughly correspond to the difference in phenotype for carriers of the shortest vs. longest VNTR alleles.  For additional details, please see the corresponding manuscript.
