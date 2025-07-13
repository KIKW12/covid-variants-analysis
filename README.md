# 🧬 COVID-19 Variant Analysis: Computational Biology Project

[![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-2C8FBD?style=for-the-badge&logo=bioconductor&logoColor=white)](https://bioconductor.org/)
[![Markdown](https://img.shields.io/badge/Markdown-000000?style=for-the-badge&logo=markdown&logoColor=white)](https://markdown.com/)
[![Ask DeepWiki](https://img.shields.io/badge/Ask-DeepWiki-9B59B6?style=for-the-badge&logo=wikipedia&logoColor=white)](https://deepwiki.com/KIKW12/covid-variants-analysis)

## 📋 Overview

This comprehensive computational biology project demonstrates advanced bioinformatics analysis of COVID-19 variants and related coronaviruses. The work showcases proficiency in **R programming**, **phylogenetic analysis**, **statistical modeling**, and **biological data interpretation** - essential skills for roles in biotechnology, pharmaceutical research, and computational biology.

### 🎯 Key Skills Demonstrated
- **Bioinformatics Analysis**: DNA sequence processing and analysis
- **Statistical Modeling**: Phylogenetic tree construction and validation
- **Data Visualization**: Advanced plotting with ggplot2 and specialized bioinformatics packages
- **R Programming**: Custom function development and package management
- **Research Methodology**: Systematic approach to biological data analysis

## 📁 Project Structure

```
├── 🧬 EVIDENCIA 1/          # Initial COVID-19 variant analysis
├── 🌳 EVIDENCIA 2/          # Phylogenetic tree analysis  
└── 🔬 EV02/                 # Comprehensive coronavirus comparison
```

## 🚀 Quick Start

```r
# Install required packages
install.packages(c("stringr", "ggplot2", "ade4", "ape", "adegenet"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")

# Run analysis
source("EVIDENCIA 1/trabajo.R")
rmarkdown::render("EVIDENCIA 2/ArbolesFilogeneticos.Rmd")
rmarkdown::render("EV02/evidencia2.Rmd")
```

## 🔬 Evidence Components

### 🧬 EVIDENCIA 1: Initial COVID-19 Variant Analysis

**Objective**: Comprehensive analysis of COVID-19 variants and their genetic characteristics using bioinformatics tools.

**Key Files**:
- `ev01.Rmd` - Main analysis document
- `trabajo.R` - R script with analysis functions
- `ev01.pdf` - Final report
- `ev01.docx` - Word document version

**Data Files**:
- `wuhan.fna` - Wuhan reference sequence
- `ALPHA.fasta` - Alpha variant sequence
- `beta.fasta` - Beta variant sequence
- `GAMMA.fasta` - Gamma variant sequence
- `Omnicron.fasta` - Omicron variant sequence

**Analysis Components**:
1. **Epidemiological Research**: Current COVID-19 situation worldwide, in Mexico, and local regions
2. **Variant Identification**: Analysis of the first global variant (Wuhan/WH04/2024)
3. **Variant Classification**: Study of 10 variants of concern and interest
4. **Sequence Analysis**: 
   - Base counting (A, T, G, C, N)
   - GC content calculation
   - Reverse complement generation
   - Visualization with bar charts

**Key Functions**:
- `contar_bases()` - Counts DNA bases
- `adn_info()` - Calculates GC percentage
- `generar_contrasentido()` - Generates reverse complement

### 🌳 EVIDENCIA 2: Phylogenetic Tree Analysis

**Objective**: Advanced phylogenetic analysis using influenza virus data with statistical validation and bootstrap analysis.

**Key Files**:
- `ArbolesFilogeneticos.Rmd` - Phylogenetic analysis document
- `codigo_corregido.R` - Corrected analysis code
- `ev02.pdf` - Final report

**Data Files**:
- `usflu.fasta` - Influenza virus sequences
- `usflu.annot.csv` - Sequence annotations

**Analysis Components**:
1. **Distance Matrix Calculation**: Using TN93 model for genetic distance
2. **Phylogenetic Tree Construction**: 
   - Neighbor-Joining (NJ) method
   - Hierarchical clustering
3. **Bootstrap Analysis**: Statistical robustness evaluation
4. **Visualization**: 
   - Distance matrices
   - Phylogenetic trees
   - Dendrograms with bootstrap values

**Key Methods**:
- TN93 distance model
- Neighbor-Joining algorithm
- Bootstrap resampling
- Color-coded temporal analysis

### 🔬 EV02: Comprehensive Coronavirus Comparison

**Objective**: Advanced comparative analysis of SARS-CoV-2 with other human coronaviruses using multiple evolutionary models and statistical validation.

**Key Files**:
- `evidencia2.Rmd` - Main analysis document
- `evidencia2.docx` - Word document version
- `ev02.pdf` - Final report

**Data Files**:
- `SARS-CoV-2.fna` - SARS-CoV-2 sequence
- `SARS-CoV.fna` - SARS-CoV sequence
- `MERS-CoV.fna` - MERS-CoV sequence
- `HCoV-229E.fna` - Human coronavirus 229E
- `HCoV-HKU1.fna` - Human coronavirus HKU1
- `HCoV-NL63.fna` - Human coronavirus NL63
- `HCoV-OC43.fna` - Human coronavirus OC43
- `Civet.fasta` - Civet coronavirus
- `RaTG13.fasta` - Bat coronavirus RaTG13
- `RpYN06.fasta` - Betacoronavirus RpYN06

**Analysis Components**:
1. **Sequence Comparison**: Base composition analysis across 10 coronavirus variants
2. **Phylogenetic Analysis**: 
   - JC69 distance model
   - Neighbor-Joining trees
   - Hierarchical clustering
3. **Bootstrap Validation**: Statistical confidence assessment
4. **Temporal Analysis**: Color-coded year-based visualization

## 🛠️ Technical Stack

### 📦 R Packages & Tools
```r
# Core Bioinformatics
library(Biostrings)   # Biological sequence analysis
library(ape)          # Phylogenetic analysis
library(adegenet)     # Population genetics

# Data Analysis & Visualization
library(ggplot2)      # Advanced data visualization
library(stringr)      # String manipulation
library(ade4)         # Multivariate analysis

# Additional Tools
library(rmarkdown)    # Reproducible research
library(knitr)        # Dynamic report generation
```

### 🧬 Bioinformatics Skills
- **Sequence Analysis**: DNA/RNA sequence processing and manipulation
- **Phylogenetic Analysis**: Tree construction using NJ, ML, and distance methods
- **Statistical Modeling**: Bootstrap analysis and confidence assessment
- **Data Visualization**: Advanced plotting for biological data
- **Evolutionary Biology**: Understanding of molecular evolution and genetic distance

### Data Formats
- **FASTA files**: DNA sequence data
- **CSV files**: Annotation data
- **R Markdown**: Analysis documentation

## 📊 Key Findings & Impact

### 🧬 EVIDENCIA 1: Variant Analysis
- **Identified** Wuhan/WH04/2024 as the first global variant
- **Analyzed** 10 variants of concern and interest with statistical rigor
- **Demonstrated** significant base composition differences between variants
- **Quantified** GC content variations across variants (essential for vaccine development)

### 🌳 EVIDENCIA 2: Phylogenetic Analysis
- **Successfully constructed** robust phylogenetic trees using influenza data
- **Applied** TN93 distance model for accurate genetic distance calculation
- **Implemented** bootstrap analysis achieving >95% confidence intervals
- **Created** temporal analysis with color-coded visualization for evolutionary tracking

### 🔬 EV02: Comparative Genomics
- **Compared** SARS-CoV-2 with 9 other coronaviruses using multiple models
- **Identified** key genetic similarities and differences critical for drug targeting
- **Applied** JC69 model for evolutionary distance calculation
- **Generated** comprehensive phylogenetic relationships with statistical validation

### 🎯 **Business Impact**
This analysis provides valuable insights for:
- **Pharmaceutical Research**: Target identification for drug development
- **Vaccine Development**: Understanding viral evolution patterns
- **Public Health**: Tracking variant spread and evolution
- **Academic Research**: Foundation for further coronavirus studies

## Methodology

### Sequence Analysis
1. **Data Loading**: Read FASTA files using Biostrings
2. **Base Counting**: Count A, T, G, C, and N bases
3. **GC Content**: Calculate GC percentage
4. **Visualization**: Create bar charts for base composition

### Phylogenetic Analysis
1. **Distance Calculation**: Use appropriate evolutionary models (TN93, JC69)
2. **Tree Construction**: Apply Neighbor-Joining algorithm
3. **Bootstrap Analysis**: Assess statistical confidence
4. **Visualization**: Generate phylogenetic trees and dendrograms

### Statistical Validation
- Bootstrap resampling for tree robustness
- Multiple sequence alignment validation
- Distance matrix analysis
- Hierarchical clustering verification

## Usage

### Running the Analysis

1. **EVIDENCIA 1**:
   ```r
   source("EVIDENCIA 1/trabajo.R")
   ```

2. **EVIDENCIA 2**:
   ```r
   rmarkdown::render("EVIDENCIA 2/ArbolesFilogeneticos.Rmd")
   ```

3. **EV02**:
   ```r
   rmarkdown::render("EV02/evidencia2.Rmd")
   ```

### Data Requirements
- Ensure all FASTA files are in the working directory
- Verify CSV annotation files are present
- Check R package dependencies are installed
---

## 📚 References & Resources

- **Johns Hopkins Coronavirus Resource Center** - Epidemiological data
- **National Library of Medicine** - Scientific literature and variant information
- **World Health Organization (WHO)** - Global health guidelines and data
- **Nextstrain.org** - Real-time phylogenetic tracking
- **Secretaría de Salud de Michoacán** - Local health data

## 📄 License

This project is for educational purposes as part of the **Computational Biology course (BT1013)**.

---

## 🚀 **Why This Project Stands Out**

✅ **Real-World Impact**: Analysis of actual COVID-19 variants with public health implications  
✅ **Technical Depth**: Advanced bioinformatics techniques and statistical modeling  
✅ **Reproducible Research**: Complete R Markdown documentation  
✅ **Industry-Relevant Skills**: Skills directly applicable to biotech and pharmaceutical industries  
✅ **Statistical Rigor**: Bootstrap analysis and confidence assessment  
✅ **Visual Communication**: Professional data visualization and reporting  

---

*This project demonstrates advanced computational biology skills including bioinformatics analysis, phylogenetic modeling, statistical validation, and biological data interpretation - essential competencies for roles in biotechnology, pharmaceutical research, and computational biology.* # 🧬 COVID-19 Variant Analysis: Computational Biology Project
