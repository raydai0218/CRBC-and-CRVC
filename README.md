# Crop root bacterial and viral genomes unveil novel species and microbiome patterns
Workflows and code for crop root bacterial and viral genome collection construction and analysis  

## Summary
Reference genomes of root microbes are essential for metagenomic analyses and mechanistic studies of crop root microbiomes. Combining high-throughput bacterial cultivation with metagenomic sequencing, we constructed comprehensive bacterial and viral genome collections from the roots of wheat, rice, maize, and Medicago. The crop root bacterial genome collection (CRBC) significantly expands the quantity and phylogenetic diversity of publicly available crop root bacterial genomes, with 6,699 bacterial genomes (68.9 % from isolates) and 1,817 potentially novel species, expanding crop root bacterial diversity by 290.6%. The crop root viral genome collection (CRVC) contains 9,736 nonredundant viral genomes, with 1,572 previously unreported genus-level clusters in crop roots. From these data, we identified conserved bacterial functions enriched in root microbiomes across soils and host species and uncovered previously unexplored bacteria–virus connections in crop root ecosystems. Together, the CRBC and CRVC serve as valuable resources (www.cropmicrobiome.com/) for investigating microbial mechanisms and applications, supporting sustainable agriculture.

## Material and data availability
- Details on CRBC bacterial isolate strains can be accessed through the Guangdong Microbial Culture Collection Center (GDMCC; https://gdmcc.net/#/index). All cultivated bacterial strains are available upon request by contacting Dr. Yang Bai (ybai@pku.edu.cn).
- The bacterial and viral genomes, along with annotation files, and other related data are accessible on Zenodo (accessions: 14091751, 14095420, 13918137, and 13939322). Raw sequencing data for bacterial genomes and root metagenomes are available on NCBI (accession PRJNA1183633 and PRJNA1184367), and maize rhizosphere metagenome data are available on ENA (accession PRJEB77048). The genomic data in this work are also accessible through our website (http://www.cropmicrobiome.com/).

## Workflow of CRBC
1. CRBC genome construction
2. Taxonomic classification and genome novelty
3. Bacterial gene prediction and functional annotation
4. BGC diversity and novelty analysis
5. Defense system and CRlSPR array
6. Calculate bacterial abundance, prevalence and growth rate
7. Construction of Kraken2 database
  
## Workflow of CRVC
1. CRVC genome construction
2. Viral genome taxonomic classification
3. Define novel viral species and genera
4. Viral gene prediction and functional annotation
5. Viral lifestyle prediction
6. Host prediction
