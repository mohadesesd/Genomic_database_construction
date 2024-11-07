## Introduction

This documentation outlines the structure and usage of the `genomic_variants.db` SQLite3 database and the accompanying **Genome Browser** and **Tkinter GUI**, which are web-based and desktop-based interfaces designed for easy querying and visualization of genomic data. The database consolidates variant data from VCF files and integrates ClinVar annotations to provide comprehensive genomic insights.

---

## Overview

The `genomic_variants.db` database is a centralized repository that manages genomic variant data, incorporating information from VCF files and ClinVar annotations. This integration facilitates detailed genomic analyses and supports various querying methods through user-friendly interfaces.

---

## Database Structure

The `genomic_variants.db` database is organized into four main tables, each serving a distinct purpose in managing and relating genomic data.

### Tables

1. **variants**
2. **samples**
3. **genotype**
4. **clinvar_annotations**

### 1. variants

Stores general information about each genomic variant, including key quality metrics and details extracted from VCF INFO fields.

**Columns:**

| **Column**          | **Description**                                              |
|---------------------|--------------------------------------------------------------|
| `variant_id`        | Unique identifier for each variant                           |
| `chrom`             | Chromosome number (e.g., '1', '2', ..., 'X', 'Y')            |
| `pos`               | Position on the chromosome                                   |
| `ref`               | Reference allele                                             |
| `alt`               | Alternate allele(s)                                          |
| `qual`              | Quality score of the variant                                 |
| `filter`            | Filter status (e.g., "PASS", "q10")                          |
| `DP`                | Read Depth                                                   |
| `AF`                | Allele Frequency                                             |
| `AC`                | Allele Count                                                 |
| `ExcessHet`         | Excess heterozygosity statistic                              |
| `FS`                | Strand bias statistic                                        |
| `MQ`                | Mapping Quality                                              |
| `QD`                | Quality by Depth                                             |
| `SOR`               | Symmetry Odds Ratio                                          |

### 2. samples

Lists the sample names that correspond to specific genotypes.

**Columns:**

| **Column**     | **Description**                                  |
|----------------|--------------------------------------------------|
| `sample_id`    | Unique identifier for each sample                |
| `sample_name`  | Name of the sample (e.g., 'Sample_1')            |

### 3. genotype

Links genotypes to specific variants and samples.

**Columns:**

| **Column**      | **Description**                                     |
|-----------------|-----------------------------------------------------|
| `genotype_id`   | Unique identifier for each genotype entry           |
| `variant_id`    | Foreign key linking to the `variants` table         |
| `sample_id`     | Foreign key linking to the `samples` table          |
| `genotype`      | Genotype information (e.g., '0/1', '1/1')           |

### 4. clinvar_annotations

Provides ClinVar-specific annotations, adding clinical context for certain variants.

**Columns:**

| **Column**               | **Description**                                              |
|--------------------------|--------------------------------------------------------------|
| `annotation_id`          | Unique identifier for each ClinVar annotation               |
| `variant_id`             | Foreign key linking to the `variants` table                 |
| `clinical_significance`  | Clinical significance (e.g., "Pathogenic", "Benign")        |
| `condition`              | Associated conditions or disease names                      |
| `review_status`          | Review status of ClinVar submission                         |
| `AF_EXAC`                | Allele Frequency in ExAC                                     |
| `CLNDISDB`               | ClinVar-disease database references                         |
| `CLNDN`                  | Disease name associated with the variant                    |
| `CLNHGVS`                | Human Genome Variation Society (HGVS) nomenclature           |
| `CLNVC`                  | Variant type (e.g., SNV, InDel)                              |
| `CLNVCSO`                | Sequence ontology term for variant                           |
| `GENEINFO`               | Gene information                                             |
| `MC`                     | Molecular consequence                                        |
| `ORIGIN`                 | Origin of the variant                                        |
| `RS`                     | Reference SNP ID (mirrors the `RS` field in VCF data)        |

---

## User Interfaces

To facilitate easy access and interaction with the genomic data, two user interfaces are provided:

1. **Genome Browser (Flask Web Interface)**
2. **Tkinter GUI (Desktop Application)**

### 1. Genome Browser (Flask Web Interface)

The **Genome Browser** is a custom-built, web-based interface developed using Flask. It allows users to seamlessly query, filter, and visualize variant data from the `genomic_variants.db` database.

### 2. Tkinter GUI

In addition to the web-based Genome Browser, a desktop graphical user interface (GUI) built using Tkinter is available. This Tkinter GUI offers an alternative way to query and interact with the genomic data in a user-friendly environment directly on your local machine.

### 3. Comparison of Interfaces

| **Feature**              | **Flask Web Interface**                    | **Tkinter GUI**                       |
|--------------------------|--------------------------------------------|---------------------------------------|
| **Platform**             | Browser-based                              | Desktop-based                         |
| **Installation**         | Requires Flask and browser                 | Requires Python and Tkinter           |
| **Usability**            | Accessible from any device with a browser  | Accessible on the local machine       |
| **Search and Filtering** | Multi-filter search with sorting options   | Multi-filter search with export options |
| **Data Export**          |                                            | Supports CSV export                   |
| **Customization**        | Uses HTML/CSS for advanced styling         | Limited styling with Tkinter          |
| **Pagination**           | Built-in pagination for large datasets     | Scrollable table for viewing results  |

Both interfaces offer effective ways to query and interact with the genomic database. The Flask interface is ideal for collaborative use or remote access, while the Tkinter GUI is suitable for quick desktop use without a browser dependency.

---

## Setup Instructions

Follow these steps to set up and run the Genomic Database and Genome Browser on your local machine.

### Prerequisites

- **Python 3.7 or Higher**: Ensure Python is installed on your system. Download it from [python.org](https://www.python.org/downloads/).
- **pip**: Python package installer, typically included with Python installations.

### Installation
1. **Create a Virtual Environment** (Optional but Recommended):

   ```bash
   python -m venv venv
   ```

   Activate the virtual environment:

   - **On macOS/Linux:**

     ```bash
     source venv/bin/activate
     ```

   - **On Windows:**

     ```bash
     venv\Scripts\activate
     ```

2. **Install Required Dependencies**:

   Ensure you are in the project directory and run:

   ```bash
   pip install -r requirements.txt
   ```

   This command installs all necessary Python packages listed in the `requirements.txt` file.

### Running the Applications

#### 1. Initialize the Database (If Not Already Initialized)

Ensure that the `genomic_variants.db` SQLite database is present in the project directory. If not, run any provided scripts to set up the database schema and import data.

#### 2. Start the Flask Application (Genome Browser)

```bash
python app.py
```

This command launches the Flask server hosting the Genome Browser.

#### 3. Run the Tkinter GUI

```bash
python GUI_tkinter
\.py
```

This command launches the Tkinter desktop application.

#### 4. Access the Genome Browser

Open your web browser and navigate to [http://127.0.0.1:5000](http://127.0.0.1:5000) to access the Genome Browser interface.

---

## Query Examples

For users who prefer interacting directly with the SQLite database, here are some example SQL queries:

### 1. Retrieve All Variants on a Specific Chromosome and Position

```sql
SELECT *
FROM variants
WHERE chrom = '17' AND pos = 123456;
```

**Description:** This query fetches all details of variants located on chromosome 17 at position 123,456.

### 2. Find Variants with Specific Clinical Significance from ClinVar

```sql
SELECT v.chrom, v.pos, v.ref, v.alt, c.clinical_significance
FROM variants AS v
JOIN clinvar_annotations AS c ON v.variant_id = c.variant_id
WHERE c.clinical_significance = 'Pathogenic';
```

**Description:** Retrieves chromosome, position, reference and alternate alleles, and clinical significance for all variants marked as "Pathogenic" in ClinVar.

### 3. List All Genotypes for a Specific Sample

```sql
SELECT s.sample_name, g.genotype, v.chrom, v.pos, v.ref, v.alt
FROM genotype AS g
JOIN samples AS s ON g.sample_id = s.sample_id
JOIN variants AS v ON g.variant_id = v.variant_id
WHERE s.sample_name = 'Sample_1';
```

**Description:** Lists all genotype information for "Sample_1", including the associated chromosome, position, reference, and alternate alleles.

---

## Using the User Interfaces

### Genome Browser

The **Genome Browser** provides a user-friendly interface to explore the genomic data without direct interaction with the database. Here's how to utilize its features:

1. **Access the Interface**:

   Open your web browser and navigate to [http://127.0.0.1:5000](http://127.0.0.1:5000).

2. **Search and Filter**:

   - **Chromosome**: Select the chromosome of interest (e.g., 17).
   - **Position**: Enter a specific genomic position or a range.
   - **Clinical Significance**: Filter variants based on their clinical impact (e.g., "Pathogenic", "Benign").
   - **Sample Genotype**: Choose specific genotype information related to samples.
   - **Additional Filters**: Utilize other available fields such as allele frequency, gene information, etc., to refine your search.

3. **View Results**:

   - **Table View**: Examine variant details in a structured table format with sortable columns.
   - **Data Visualization**: If available, view charts or plots that illustrate data distributions and patterns.

4. **Navigate Through Data**:

   - **Pagination**: Use pagination controls to navigate through large datasets efficiently.

### Tkinter GUI

The **Tkinter GUI** serves as a straightforward, standalone application for users who prefer a desktop interface over a web-based one. It provides various features for querying and exporting genomic variant data, leveraging the `genomic_variants.db` SQLite database.

#### Setup and Usage

1. **Prerequisites**:

   - Ensure that Python 3 and Tkinter are installed on your system. Tkinter typically comes pre-installed with Python, but if not, install it using the package manager specific to your operating system.

2. **Running the Application**:

   - Run the Tkinter application by executing the GUI script:

     ```bash
     python GUI_tkinter.py
     ```

   - The application window will open, allowing you to interact with the database.

3. **Using the Tkinter GUI**:

   - **Filter Criteria**: Select filter fields (e.g., chromosome, position) and define criteria (e.g., equals, contains) to narrow down your search.
   - **Apply Filter**: Click the "Apply Filter" button to execute the query and view the results in the table.
   - **Add Additional Criteria**: Use the "Add Criteria" button to apply multiple filters at once.
   - **Export Results**: Click the "Export Results" button to save the query results in a CSV format.

4. **Navigating Results**:

   - The results table allows easy navigation through rows of data. You can scroll vertically and horizontally to view all available fields.
---

## Notes

- **File Paths**: Ensure that file paths within scripts and configurations are correctly specified to avoid any runtime errors. Adjust configurations as needed based on your system environment.

- **Tkinter Script Paths**: Specifically, ensure all file paths in the Tkinter script are updated to match the location of `genomic_variants.db` on your machine. This application provides a convenient desktop-based alternative for querying genomic variant data, making it accessible even without web server deployment.

---
