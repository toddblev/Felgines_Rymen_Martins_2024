# Felgines_Rymen_Martins_2024
Code and example datasets used in the publication Felgines, Rymen, Martins et al. (2024) _Nature Communications_
Full datasets are available as supplementary data. [10.1038/s41467-024-54268-0](https://doi.org/10.1038/s41467-024-54268-0)

## Available scripts

### Immunoprecipitation and mass spectrometry (IP-MS) analyses

- Volcano plots
  - volcano_plot/volcano_example.R
    - Author: Calvin Matteoli
    - Description: Base script used to generate volcano plots for Fig. 2C, 2D, 3A and 3C
    - Example dataset: Data from Fig. 2C

- Balloon plots
  - balloon_plot_simple/balloonplot_example.R
    - Author: Calvin Matteoli
    - Description: Base script used to generate balloon plots for Fig. 3B and 3D
    - Example dataset: Data from Fig. 3D
  - balloon_plot_split/balloonplot_split_example.R
    - Author: Calvin Matteoli
    - Description: Script used to generate split p-value balloon plots
    - Example dataset: Fig. 2B
   
### Small interfering RNA (siRNA) analyses

- Violin plots
  - siRNA_scripts/Figure_4A_violin_plot/Figure_4A_violin_plot.R
    - Author: Guanghui Xu
    - Description: Base script used to generate violin plots for Fig. 4A
    - Example dataset: Data from Fig. 4A

- Half-violin plots
  - siRNA_scripts/Figure_4B_half_violin_plot/Figure_4B_half_violin_plot.R
    - Author: Guanghui Xu
    - Description: Base script used to generate half-violin plots for Fig. 4B
    - Example dataset: Data from Fig. 4B

- siRNA mapping and quantification
  - siRNA_scripts/siRNA_mapping_and_quantification/splitTagDirectoryByLength.dev2.pl
  - siRNA_scripts/siRNA_mapping_and_quantification/JSON_findPerfectMatches_and_TerminalMisMatches_v3
    - Author: Guanghui Xu
    - Description: Custom scripts used in combination with standard sRNA-seq alignment tools

## Paper Information
_Nature Communications_(2024) 15:10298 [10.1038/s41467-024-54268-0](https://doi.org/10.1038/s41467-024-54268-0)
### CLSY docking to Pol IV requires a conserved domain critical for small RNA biogenesis and transposon silencing.

##### Luisa Felgines<sup>1,a</sup>, Bart Rymen<sup>1,a</sup>, Laura M. Martins<sup>2,a</sup>, Guanghui Xu<sup>2</sup>, Calvin Matteoli<sup>1</sup>, Christophe Himber<sup>1</sup>, Ming Zhou<sup>2,3,b</sup>, Josh Eis<sup>2</sup>, Ceyda Coruh<sup>2</sup>, Marcel Böhrer<sup>1</sup>, Lauriane Kuhn<sup>4</sup>, Johana Chicher<sup>4</sup>, Vijaya Pandey<sup>5</sup>, Philippe Hammann<sup>4</sup>, James Wohlschlegel<sup>5</sup>, Florent Waltz<sup>6</sup>, Julie A. Law<sup>2,7,c</sup>, Todd Blevins<sup>1,c</sup>

1.	Institut de Biologie Moléculaire des Plantes, CNRS, Université de Strasbourg, Strasbourg, F-67084, France
2.	Plant Molecular and Cellular Biology Laboratory, Salk Institute for Biological Studies, La Jolla, CA, 92037, USA
3.	State Key Laboratory of Plant Environmental Resilience, College of Life Sciences, Zhejiang University, Hangzhou 310058, China
4.	Institut de Biologie Moléculaire et Cellulaire, CNRS, Plateforme Protéomique Strasbourg-Esplanade, Strasbourg, F-67084, France
5.	Department of Biological Chemistry, University of California, Los Angeles, CA, 90095, USA
6.	Biozentrum, University of Basel, CH-4056, Basel, Switzerland
7.	Division of Biological Sciences, University of California, San Diego, La Jolla, CA 92093

a. Contributed equally to the work  
b. Current affiliation  
c. Co-Corresponding authors: jlaw@salk.edu, todd.blevins@ibmp-cnrs.unistra.fr  

#### Abstract
Eukaryotes must balance the need for gene transcription by RNA polymerase II (Pol II) against the danger of mutations caused by transposable element (TE) proliferation. In plants, these gene expression and TE silencing activities are divided between different RNA polymerases. Specifically, RNA polymerase IV (Pol IV), which evolved from Pol II, transcribes TEs to generate small interfering RNAs (siRNAs) that guide DNA methylation and block TE transcription by Pol II. While the Pol IV complex is recruited to TEs via SNF2-like CLASSY (CLSY) proteins, how Pol IV partners with the CLSYs remains unknown. Here we identified a conserved CYC-YPMF motif that is specific to Pol IV and is positioned on the complex exterior. Furthermore, we found that this motif is essential for the co-purification of all four CLSYs with Pol IV, but that only one CLSY is present in any given Pol IV complex. These findings support a “one CLSY per Pol IV” model where the CYC-YPMF motif acts as a CLSY-docking site. Indeed, mutations in and around this motif phenocopy _pol iv_ null and _clsy_ quadruple mutants. Together, these findings provide structural and functional insights into a critical protein feature that distinguishes Pol IV from other RNA polymerases, allowing it to promote genome stability by targeting TEs for silencing.
