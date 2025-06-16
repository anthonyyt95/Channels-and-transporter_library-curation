# Gene Curation and Annotation Pipeline (MATLAB)

This repository contains MATLAB scripts for processing and curating gene lists, especially solute carrier (SLC) family genes and immune-related genes. The pipeline integrates data from public databases (e.g., SLC BioParadigms, HGNC, UniProt), annotates genes with functional and Gene Ontology (GO) information, and filters for biologically relevant candidates using keyword-based searches and expression filters.

---

## 🧬 Purpose

To facilitate large-scale curation and functional annotation of gene lists with applications in:
- Solute carrier/transporter discovery
- Immune cell function screening
- Ion channel and membrane protein prioritization

---

## 🔧 Features

### 1. Data Curation
- Parses gene lists copied from:
  - **SLC BioParadigms**
  - **HGNC** exports
- Removes pseudogenes (e.g., entries starting with `ENSG`)
- Deduplicates entries while preserving order

### 2. Keyword-Based Filtering
- Filters gene names or descriptions for terms of interest (e.g., `channel`, `ion`, `transporter`)
- Supports flexible search terms defined by the user
- Also searches annotated UniProt descriptions and GO categories

### 3. UniProt Annotation
- Scrapes UniProt entries for:
  - Species-specific validation (e.g., *Homo sapiens*)
  - Short functional summaries
  - GO Molecular Function and Biological Process terms
  - Gene aliases

### 4. Batch Processing
- Splits gene list into manageable chunks (e.g., 500 genes at a time)
- Outputs Excel files with annotated gene information:
  - `a1.xlsx`, `a2.xlsx`, ..., `a15.xlsx`

### 5. Set Operations
- Excludes overlapping genes between multiple lists (e.g., removing already studied genes)
- Identifies genes not present in target exclusion lists

### 6. Cleanup and Normalization
- Strips HTML from scraped strings
- Consolidates aliases to prevent redundancy in gene naming

---

## 💻 Inputs

- `list`, `list1`, `list2`: Cell arrays of gene names or annotations
- `query`: Keyword or term of interest
- `ref`: Reference annotation table (e.g., HGNC export with synonyms)

---

## 📤 Outputs

- Annotated gene lists (Excel or `.mat` format)
- Subsets of genes with/without functional terms
- Description and GO-based filtered lists
- Summary tables of expression or function (where applicable)

---

## 🚀 Example Use Cases

- Identify all human SLC transporters with “ion channel” in their UniProt description
- Exclude known immune cell transporters to find novel candidates
- Annotate gene function across 7,000+ genes using UniProt and GO
- Remove pseudogenes and duplicate entries from downloaded gene panels

---

## 🗂️ File Organization

Some major functions used across this pipeline:
- `getGeneInfo.m`: Queries UniProt for short descriptions
- `extractHTMLText.m`: Cleans HTML content from web-scraped text
- `importdata`: Loads annotation or gene list files
- `xlswrite`: Writes output annotations to `.xlsx` files

---

## 📝 Notes
