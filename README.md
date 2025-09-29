# Reproducible Code for:

**IL-15 Drives CD8⁺ T Cell–Mediated Neurotoxicity in Tauopathy**

---

## Citation

Saef Izzy¹²⁹*, Patrick Da Silva²*, Taha Yahya²*, Maryam Hazim Al Mansi², Francesca Percopo², Daniel Schmidlin⁵, Harm Tjebbe Imkamp², Ronaldo S. Francisco Jr², Michael Aronchik², Sybren Lein Nikola Maas²³⁴, Liza Marie Morsett⁵, Jonathan R. Christenson², Marilia Garcia De Oliveira², Masen Boucher⁵, Ali Al Jarrah¹, Kuan-Jung Lu², Farid Radmanesh¹, Debanjan Mukherjee², Leah C. Beauchamp², Aikaterini Kalemaki¹, Somen Kumar Mistri², Kohei Yamada², Mahsa Khayat-Khoei², Juan Pablo de Rivero Vaccari, Edy Kim, Stelios M. Smirnakis¹, Changiz Geula⁷, Suzanne E. Hickman⁵, Hyun-Sik Yang, Rafael M. Rezende², Howard L. Weiner²#, Joseph El Khoury⁵⁸⁹#

\* Equally contributed co-first authors  

¹ Department of Neurology, Brigham and Women’s Hospital, Boston, MA  
² Ann Romney Center for Neurological Diseases, Brigham and Women’s Hospital, Harvard Medical School, Boston, MA  
³ Department of Pathology, Leiden University Medical Center, Leiden University, Leiden, The Netherlands  
⁴ Department of Pathology, University Medical Center, Utrecht University, Utrecht, The Netherlands  
⁵ Center for Immunology & Inflammatory Diseases, Massachusetts General Hospital, Boston, MA  
⁷ Northwestern University School of Medicine, Chicago, IL  
⁸ Department of Medicine, Division of Infectious Diseases, Massachusetts General Hospital, Boston, MA  
⁹ Harvard Medical School, Boston, MA  

**Corresponding Authors:**  
- Joseph El Khoury, MD — jelkhoury@mgh.harvard.edu  
- Saef Izzy, MD — sizzy@bwh.harvard.edu  

---

## Repository Description

This repository contains reproducible scripts used in the analyses and figure generation for the above publication. Scripts cover preprocessing, differential expression analysis, and visualization for scRNA-seq, snRNA-seq, and mouse model datasets.

---

## Included Scripts

- `DESeq2_mouse_aIL15.R` — Differential expression analysis in mouse models with anti-IL15 treatment  
- `figs_scRNAseq_AD_paper2.R` — Figure generation for AD scRNA-seq datasets (set 2)  
- `figures_scRNAseq_AD_paper.R` — Figure generation for main AD scRNA-seq datasets  
- `scRNAseq_Izzy_mice.R` — scRNA-seq analysis for Izzy et al. mouse models  
- `scRNAseq_izzy_public_AD2-3b.R` — Integration and analysis of public AD scRNA-seq datasets  
- `scRNAseq_public_CSF_AD1.R` — scRNA-seq analysis of public CSF datasets in AD  
- `scRNAseq_rafa_lecanemab_mouse3.R` — scRNA-seq analysis of Lecanemab mouse models (Rafael dataset)  
- `snRNAseq_ALS_brain_PFC_MCX_syn51105515.R` — snRNA-seq analysis of ALS brain tissue (synapse dataset)  
- `snRNAseq_FTD_brain_PFC_MCX_syn51105515.R` — snRNA-seq analysis of FTD brain tissue (synapse dataset)  
- `snRNAseq_public_brain_AD1-2-5.R` — snRNA-seq analysis of multiple public AD brain datasets  

---

## Notes

- Scripts are designed for reproducibility; raw datasets are not hosted here.  
- Data can be accessed through their respective repositories as referenced in the manuscript.  
- Please cite the paper above when using these scripts or derived analyses.  

---
