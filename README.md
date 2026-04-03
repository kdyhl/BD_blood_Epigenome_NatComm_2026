
# Scripts for blood epigenome in bipolar disorder published in Nature communications, 2026

## environment 
most R codes is under R-3.6 unless it is mentioned in the script

---
## data 
Raw ChIP-seq data (in BAM format), processed count matrices, genotype matrix, and metadata including sex, age, race, and ethnics are available at GEO under accession number GSE316497.
Other data including brain eQTL, hQTL mapped could be found at Zenodo  

---
## folders overview
- `./deconv`: deconvolution — cell fraction estimation  
- `./eQTLvsGWAS`: genetic evidence of candidate genes based on eQTL–GWAS colocalization  
- `./haQTL`: hQTL mapping  
- `./haQTLvseQTL`: CRE–gene linking inference based on hQTL–eQTL MR and colocalization  
- `./haQTLvsGWAS`: driver CRE inference based on hQTL–GWAS MR and colocalization  
- `./heritability`: heritability enrichment using LDSC  
- `./hiC.DiffPeakandDEGandTF`: integrate Hi-C data to identify dCREs, infer BD epigenetic signatures, and perform drug repurposing  
- `./peakMergeNormal`: data processing and differential signal mapping  
- `./peakVariationAcrossTiss`: assess CRE activity across tissues  
- `./PRS`: generate polygenic risk scores  
- `./sampleManifold`: analyze patient heterogeneity 

---

## results generation based on figures

## Fig 1 CRE detection and annotation

### CRE counts matrix and peak annotation
./peakMergeNormal
- a.mergePeak2ReadsCounts_V4.2_filteredMergedPeak_H3K27ac.R
- a.mergePeak2ReadsCounts_V4.2_filteredMergedPeak.R
- a.2.normalizePeakCounts.V1.4_H3K27ac.R
- a.2.normalizePeakCounts.V1.4_other.R
- b.annotatePeak.byHomer.R => #Supp Fig.1a

output folder: ./peakMergeNormal/a_2_peakHeightsNormalized_V1.4/ #all CREs detected

### CRE activity across different tissues in EpiMap
./peakVariationAcrossTiss/
- a.pairwiseOvlpPeaks.R
- a.mergePeaksAcross3actMarks.R
- a.2.extractRefSignalForMergedPeaksFromEpimap_V1.5_HM.R
- a.3.identifyAREModulesAcrossEpimap_V1.6.1_merge3ActiveMarks.R
- a.4.getAnnotationForModule_3actHM_onEpiMap.R =>#Fig 1b

- a.4.runTFMotifsEnrichForAREModlue_byHomer.sh
- a.4.getTFMotifsEnrichedForAREModule.byHomer_V1.1_EpiMapMotif.R
- a.4.2.visTFMotifsEnrichedInAREModule.byHomer_EpimapMotif.R => #Supp Fig. 1c&d

- a.4.annotateModuleByrGREAT_3actMarks_Epimap.R => #Supp Fig. 1e

---
## Fig 2 CRE differential signal

### presence between case and control
./peakMergeNormal
- b.cmpPeakPresenceInd.R #(Supp Fig. 2)

### cell fraction estimation
./deconv
step a, signature peaks
- a.selectPeakAndSimulation_V1.7.2_colineariyAuto_rowNorm_6Celltypes.R 
- output folder: ./deconv/a_selectPeaksAndSimulation_V1.7.2_colinearAuto_rowNorm_6cellType

step b, estimate cell fraction
- b.real.estFract_V1.4.4_nnPoisR_noIntercept_6CellType.R #(Supp Fig. 3a)
- output folder: ./deconv/b_real_estFract_V1.4.4.1_nnPoisR_noIntcpt_peakNorm_6cellType_autoSig 

### differential CREs
./peakMergeNormal
- c.getDiff_byLimma_V1.4.1_sva_NoEHR.R 
- d.annotateDiffPeak_byrGREAT_V1.4.1_sva_NoEHR.R (Fig.2a)
- d.getTFMotifsEnrichForDP_byHomer.sh
- d.getTFMotifsEnrichedForDP_V1.1_EpimapMotif.R (random bg)
- d.2.visTFMotifsEnrichedForDP_EpimapMotif.R
- d.cmpBulkDPacrossHMs_byLimma_V1.2_cmpNumber.R (Supp Fig. 5a)

### diagnose batch effect
./peakMergeNormal
step a. identifiy dCRE by permuting group labeles of batch
- c.getDiff_byLimma_V1.4.1.1_sva_NoEHR_batchPerm_paral.R
- c.getDiff_byLimma_V1.4.1.1_sva_NoEHR_batchPerm.R
- c.getDiff_byLimma_V1.4.1.1_sva_NoEHR_batchPerm_Summary.R

step b. compare the differential signal when comparing one sample to samples from the same group or from the different group
- c.getDiff_byLimma_V1.4.1.2_sva_NoEHR_cmpBatch_paral.R
- c.getDiff_byLimma_V1.4.1.2_sva_NoEHR_cmpBatch.R
- c.getDiff_byLimma_V1.4.1.2_sva_NoEHR_cmpBatch_Summary.R

step c. PCA of residuals from differential analysis
- c.getDiff_byLimma_V1.4.1_sva_NoEHR.R (Supp Fig. 4c)


### identify genes with BD epi-signatures
./hiC.DiffPeakandDEGandTF 
- a.geneBodyDP.meta_V1.2.1_allGenes.R
- a.2.cmpDPs.geneRegDomainAndGeneBody_V1.1_rmOvlp.R #(Supp Fig. 5c&d)
- a.3.visDEG.metaDP_V1.2_5conds.R #(Fig.2c, Supp Fig. 5e)

---
## Fig 3 genetics

### QTL mapping 
./haQTL
- a.identifyCovariatesForhaQTL_V1.2.3_filtPeakByMedianCV_geneINT.Peer.R
- a.2.callTisshaQTL_V1.2.3_filtPeaksByMedianCV_gINT.Peer_rmSexAge.R
- a.3.visHaQTLNumVSCovNum.R
- a.4.callTisshaQTL_fixedCov_V1.2.3_dp.bg.R
- a.5.multipleTestOnQTL_forPeer_V1.1_cutoffs.R (Supp Fig. 6c)
- b.gCRE.shared.upSet.R #(Fig. 3a)
- b.Manhattan.plot.R (Supp Fig. 6b)

### Peer vs covariates
./haQTL
- a.2.testSampPeersWithCovariates.CellFract.R (Supp Fig. 6a)

### compare with hQTL calling with top 5 genotype PCs
./haQTL
- a.4.callTisshaQTL_fixedCov_V1.2.3.1_dp.bg_genotypePC.R
- c.comphQTL_amongDifferentVersions.R

### sharing across HMs, across tissues and vs. hQTL-GWAS colocalization
based on Xiong, Nature Genetics, 2021, m6A genetics
./haQTL
- b.2.Comp.QTL.DirectionConsistEst.R (Supp Fig. 6d)
- c.comp2haQTL_fromeGTEx_V1.1_sharedgARE.R (Supp Fig. 6f)
./haQTLvsGWAS
- a.2.BDGWASColocgARE.vs.hQTLSharingWithBrain.R

### hQTL vs GWAS
./haQTLvsGWAS
- a.testBulkHaQTLVSGWAS_byColoc_V1.1_allGWAS.R
- a.testBulkHaQTLVSGWAS_byColoc_V1.1_allGWAS_summary.R (Fig.3c)
- a.testBulkHaQTLVSGWAS_byMR_V1.1_allGWAS.R
- b.colocalizeTest.VSGWAS.Vis_perLocus.acrossHM.R (Fig.3d)

### clump and summarize colocalization
./haQTLvsGWAS
- a.ClumpBDGWAS.R (Supp table 4)
- a.testBulkHaQTLVSGWAS_byColoc_V1.1_allGWAS_summary_byGWAS.loci.R

### heritability
./heritability
- a.prepAnnotForLDSC.gARE.R
- a.prepAnnotForLDSC.CREGrpAndbulkDP.R
- a2.LDSCpartitionHeritab_AREGrp.bulkDP.gCREs_slctGWAS.R
- a3.visLDSC_AREGrp.bulkDP.gARE_slctGWAS.R (Fig.3b, Supp Fig. 6e)

---
## Fig 4 causal gene and tissue sharing

### linking: hiC.reg and links
./hiC.DiffPeakandDEGandTF 
- a.buildHiCRegionLink.bulk.R   

### gLink scores
./haQTLvseQTL
1) haQTL vs eQTL coloc
- a.2.colocalizeTest.haQTLVSeQTL.coloc_V1.3.gARE_para.R
- a.2.colocalizeTest.haQTLVSeQTL.coloc_V1.3.gARE_perGene.R
- a.2.colocalizeTest.haQTLVSeQTL.coloc_V1.3.gARE_summary.R

2) haQTL vs eQTL MR
- a.2.linkGenePeak_MR-Egger_V1.1.1_gARE_para.R
- a.2.linkGenePeak_MR-Egger_V1.1.1_gARE_perGene.R
- a.2.linkGenePeak_MR-Egger_V1.1.1_gARE_summary.R

3) FMeQTL linking
- a.linkGenePeak_byFMeQTL_V1.1_gARE.R

4) summary
- b.summaryLinks.coloc.MR.PRS_V1.3_gARE.R
- b.2.cmp.geneticLinks.vsHiC_V1.1_gARE.R (linkType="gene2Neighb") (Supp Fig. 7c)


### candidate genes
based on HiC links with TSS-asscoiated hiC regions, 
both hQTL-GWAS coloc and differential signal for all histone marks
./hiC.DiffPeakandDEGandTF
gene Neigbhood and hic region annotation
- a.gneighb.VisGWASANDBulkAnnot_v1.1_bulkLink.R

DP enrichment at HiC block level: 
- a.2.compDPs.hiCBlock.bulk.R (Supp Fig. 7a&b)

### candidate genes: based on gLink, both hQTL-GWAS coloc, MR and 
differential signal, for only H3K4me1, H3K27me3, and H3K27ac
./haQTLvseQTL
- c.gneighb.VisGWASANDBulkAnnot_linkFromGenetics_V1.1_3HM.R

visualization of shared candidate genes from both linking scores
./haQTLvseQTL/
- c.2.gneighb.VisGWASANDBulkAnnot_HiCAndgLink.R (Fig. 4a)


### eQTL vs GWAS coloc signal
./eQTLvsGWAS
- a.testBlueprint.eQTLVSGWAS.coloc_byCellType.R
- a.testBlueprint.eQTLVSGWAS.coloc_para.R
- a.testBrainScRNAeQTLVSGWAS.coloc_byCellType.R
- a.testBrainScRNAeQTLVSGWAS.coloc_para.R
- a.testGTExeQTLVSGWAS_byColoc.R

### compare eQTL Coloc + gLinking/diffPeak 
./haQTLvseQTL/
- c.2.cmbPotentialGenes_V1.2_bulkAndcellsortedeQTL.R (Fig. 4b-d)

---
## Fig 5 individual heterogeneity and subtype

### cluster and annotation
./sampleManifold
- b.factorAnalysis.across5HMs.MOFA_HiC.DP_V1.2_pvalCutoff.R (Supp Fig. 8b)
- b2.cmpMOFA.Factor.withEHR.R (Supp Fig. 8c)
- b3.clusterPatientsBasedOnMOFA.factor.V1.1_allSamps.R (Fig.5a, Fig.5b, Supp Fig. 8d)
- c.cmpPatientClusterWithEHR_byMOFA_V1.1.2_addNumericFeature.R (Fig.5c)
- c.cmpAllSampClusterWithEHR_byMOFA_V1.1.2_addNumericFeature.R  #check EHR enrichment of all samples clustering 
- c.findMarkerDPforCluster_V1.2_uniq.R (Supp Fig.8e)
- c.2.annotateMarkDP_byrGREAT.R

### PRS: Clumping+Threshold
./PRS
- a.baseGWASSummary.QC.R
- b.genomicRegions.forDecouplePRS.R
- b.targetMayoBD.V1.1.noSwap.QC.R
- c.PRS.byPlink_decoupleImmuneComp.R (ld0.4.clumpDis100kb.wind2kb)
- d.cmpPRS.MayoEHR_grps.R (LD 0.4, clumbDis 100kb, wind 2kb) (Supp Fig. 9a, c, d, Fig. 5d)

### PRS: PRS-CS and compare with PRS-C+T
./PRS
- c.PRS.byPRS-CS.1.postEff.sh
- c.PRS.byPRS-CS.2.decoupleImmuneComp.R
- d.cmpPRS.MayoEHR_grps_V1.2_PRS_CS.R  
- d.cmpPRS.plink.vs.PRS_CS.R(Supp Fig. 9b)

---

## Fig 6 drugs 

### BD drug response
./peakMergeNormal/
- c.getDiff_byLimma_V1.4_sva.R
- d.annotateDiffPeak_byrGREAT_V1.4_sva.R (output dir: d_bulkDP_limma_V1.4_annotByrGREAT/ q0.2) (Fig. 9a)
- d.cmpBulkDP_BDvsDrugs.R (Supp Fig. 10a&b)

### based on BD DEGs 
./hiC.DiffPeakandDEGandTF/
- b.testEpiDEGvsCMapSignatures_V1.1_paral.R
- b.testEpiDEGvsCMapSignatures_V1.1.py

### network visualization
./hiC.DiffPeakandDEGandTF/
- b2.visualizeEpiDEGvsCMap.R (Fig.6b)

### vis signature genes across different perturbations#
./hiC.DiffPeakandDEGandTF/
- b3.extractSlctExpForSignatures.py
- b4.visualizeEpiDEGvsCMap.heatmap.R (Supp Fig. 11a)

### subtype DP -> DEGs for each group of patients
./peakMergeNormal
- c.getDiff_byLimma_V1.5.1_inflamSubgrps.R
./hiC.DiffPeakandDEGandTF/
- a.geneBodyDP.meta_V1.2.1_allGenes_subGrpDP.R
- a.2.cmpDPs.geneRegDomainAndGeneBody_V1.1_rmOvlp_subGrpDP.R
- a.3.visDEG.metaDP_V1.2_5conds_subGrpDP.R

### based on BD DEGs identified for each group of patients
./hiC.DiffPeakandDEGandTF/
- b.testEpiDEGvsCMapSignatures_subGrpDP_V1.1_paral.R
- b.testEpiDEGvsCMapSignatures_subGrpDP_V1.1.py

### network visualization for 
./hiC.DiffPeakandDEGandTF/
- b2.visualizeEpiDEGvsCMap_subGrpDP_V1.2_cytoscape.R
- output dir: b_comp2CompoundSig_cMap_inflamSubGrps/a3_V1.2_subGrps/b2_vis_V1.2_cytoscape (Fig.6c, Supp Fig. 11 b-e)
