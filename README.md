# Genomic evidence of local adaptation in Scottish *Margaritifera margaritifera*
Scripts for bioinformatics processing and analysis. Location data is not provided here due to the continued risk of pearl fishing in Scotland. Please contact us if you require it. 

__Contact:__      victoria.gillman.21@abdn.ac.uk

__Citation:__     Gillman V, Cosgrove P, Morrissey B, Lancaster L, Pritchard V, Layton K. 2025. Genomic evidence of local adaptation in Scottish freshwater pearl mussels. Published online as a preprint: [https://doi.org/10.21203/rs.3.rs-7363533/v1](https://doi.org/10.21203/rs.3.rs-7363533/v1)

__Raw Data:__ Raw genomic RADSeq data are available in the NCBI Short Read Archive under accession **PRJNA1370707**. 


### Usage Guidance

- All code is organised in the `code` folder.
- Unix and bioinformatics scripts are in `code/unix`.
- R analysis scripts are stored in `code/r_scripts`.
- Scripts are grouped by major steps and ordered accordingly:
  - Environmental data extraction: `Env_X`
  - Distance extraction between points (including "As-the-fish-swims" and Euclidean distances): `Distances_X`
  - Population genomics: `Gen_X`
  - Redundancy analysis: `RDA_X`


### Acknowledgements  
This research was funded by NatureScot and UHI Inverness. We thank Mark Coulson, Silvia Ferreira Carvalho and Jenny O’Dell for contributions to project development, Iain Sime (NatureScot) and Chris Daphne for sample collection, and Dasha Svobodova and Lydia McGill (UHI Inverness) for laboratory assistance. Paul Etter and Eric Johnson at SNPsaurus LLC facilitated the generation and quality control of RAD-seq data. VG is supported the NERC Scottish Universities Partnership for Environmental Research (SUPER) Doctoral Training Partnership (Grant reference number NE/S007342/1 to KL and website https://superdtp.st-andrews.ac.uk/) and the University of Aberdeen. 