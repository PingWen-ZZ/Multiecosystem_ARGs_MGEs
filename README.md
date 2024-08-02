# Multiecosystem_ARGs_MGEs
Analysis scripts
### 1. ARG overview

```R
# code
1.ARG_overview.r:  R script for analysing total ARG abundance, diversity and composition in different habitats

# Total abudnance, diversity and composition
input files:
total.abundance.diversity.txt：Total ARG abudance, diversity and composition in 352 samples
drug.composition.txt: Total ARG composition
output files:
ARG_total_abundance_diversity_composition.pdf

# abundance 
input files:
dat_Farmland_bysite.txt/dat_Forest_bysite.txt/dat_Grass_bysite.txt/dat_Gobi_bysite.txt/dat_Sewage_bysite.txt: Total ARG abudance and diversity in different sites
output files:
ARG_abundance_byecosystem.pdf

# composition
input files:
dat_Farmland.txt/dat_Forest.txt/dat_Grass_bysite.txt/dat_Gobi.txt/dat_Sewage.txt: ARG composition in different sites
output files:
ARG_composition_byecosystem.pdf
```



### 2. MGE overview

```R
# code
2.MGE_overview.r:  R script for analysing total MGE annotation, abundance and composition in different habitats

# MGE annotation
input files:
MGE_diversity.txt: MGE types files (plasmid, phage, ICE, transposon, IS, integron)
output files:
MGE_diversity.pdf/MGEs_diversity_percentage.pdf

# MGE annotation ratio: classified MGEs to the total MGEs
input files:
diamond_ratio.csv
output files:
MGE_database_ratio.pdf

# Total abudnance and composition
input files:
1_TotalMGE_TotalMGEtypes_basedContigsUnique.RData/1_All_MGE_mapping_corrected_by_JXforest_1kb.RData
output files:
MGE_abundance_diversity_composition.pdf

# abundance
input files:
1_TotalMGE_TotalMGEtypes_basedContigsUnique.RData/Farmland_location.txt/Forest_location.txt/Tailings_location.txt/Gobi_location.txt/Grass_location.txt/Sewage_location.txt
output files:
MGE_abundance_diversity_bylocation.pdf

# composition
input files:
2_Total_MGE_composition_bysampleType.RData/All_location_df.txt
output files:
MGE_composition_bylocation.pdf
```



### 3.Correlation between ARGs and MGEs

```R
3.Correlation_between_ARG_MGE.r: R script for analysing correlation between ARGs and MGEs using VPA and ABT

# ABT
input files:
1_ARG_MGE_combined_df.Rdata: A combined data contain ARGs and MGEs
output files:
ABT_abundance.pdf

# VPA
input files:
1_ARG_MGE_combined_df.Rdata:A combined data contain ARGs and MGEs
output files:
ARG_differentplasmids_phage_ICE_unmobilizableMGE_logscale.pdf
```



### 4. Profiles of ARG-carrying MGEs contigs

```R
4.Profiles_of_ARGs_carrying_MGEs.r: R script for analysing the composition and abundance of ARG-carrying MGEs contigs

# ARG-carrying MGEs annotation
input files:
ARG_MGE_percentage_yes_df.csv:  ARGs and MGEs located on the same contig
output files:
1_ARG_with_MGE_Number.pdf/1_ARG_with_MGE_Number_percentage.pdf

# ARG-carrying MGEs composition
input files:
drug_depth.txt
output files:
ARG_MGE_MGEtype_all.pdf

# ARG-carrying MGEs abundance
input files:
ARG_MGE_total_abundance.txt
output files:
ARG_with_MGE_abundance.pdf

# Abundance seperate
input files:
ARG_MGE_abundance_type.txt
output files:
ARG_with_MGE_abundance_differentMGE.pdf
```



### 5. Diversity of ARGs carried by various MGEs in different habitat types

```R
5.Profiles_of_ARGs_diversity_inVariousMGEs.r: R script for showing diversity of ARGs carried by various MGEs in different habitat types

# Total MGE 
input files:
1_Combined_results_diversity_total.RData/1_Ranking_total.RData: ARG diversity of different MGEs
output files:
1_MGE_carried_ARG_diversity_total.pdf

# Plasmid
input files:
1_Combined_results_diversity_plasmid.RData/1_Ranking_plasmid.RData:  ARG diversity of different plasmids
output files:
1_plasmid_carried_ARG_diversity_total.pdf

# Venn
input files:
1_Venn_combined_df.RData: special and unique ARGs in different MGEs
output files:
1_plasmids_phage_ICE_diversity_veen.pdf

#circle  bar
input files:
PlotData.farmland.txt/PlotData.sewage.txt/PlotData.grass.txt/PlotData.forest.txt/PlotData.tailings.txt/PlotData.gobi.txt: abundant ARGs in different habitat types
output files:
plasmids_phage_ICE_ARGsubtype_percentage_bar.pdf
```



### 6. Diversity of ARGs carried by plasmids associated with ISs, transposons or integrons  in different habitat types

```R
6.Profiles_of_ARGs_diversity_in_Plasmid_correlatedwith_otherMGE.r: showing diversity of ARGs carried by various MGEs in different habitat types.

# bar chart
input files:
1_Combined_results_diversity.RData/1_Ranking_total.RData
output files:
plasmid_carried_Tn_IS_Integron_diversity_total.pdf

# Venn
input files:
1_Venn_combined_df.RData: special and unique ARGs in different plasmids associated IS/transposon/integron
output files:
Plasmid_Transposon_IS_onPlasmids_diversity_venn.pdf

# circle bar
input files:
PlotData.farmland.txt/PlotData.sewage.txt/PlotData.grass.txt/PlotData.forest.txt/PlotData.tailings.txt/PlotData.gobi.txt: abundant ARGs in different habitat types
output files:
plasmids_with_TnISIntegron_ARGsubtype_percentage_abundance.pdf
```



### 7. Density

```R
7.ARG_density.r: compute ARG density in diferent plasmid types

# ARG density
input files:
Plasmids_carrying_Tn_IS_Integron_number_Length.RData/sample_name_number_df.RData/sample_name_sampleName_df.RData
output files:
DifferentPlasmids_with_Tn_IS_NumARG_perSamples.pdf

```

