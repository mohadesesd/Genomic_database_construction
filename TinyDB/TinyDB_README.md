## Overview

This documentation outlines the structure and usage of a genomic database managed with TinyDB and a web-based genome browser developed using Flask. The database focuses on variant data from VCF files and incorporates ClinVar annotations specifically for chromosome 17. The genome browser provides an intuitive interface for querying, filtering, and visualizing this data.

---

## Database Structure

The genomic database is stored as a TinyDB JSON file named `genomic_db.json`. It organizes data into variant records, enriched with ClinVar annotations for chromosome 17. Below is a breakdown of the database's structure.

### Variant Records

Each variant record in the database includes the following fields:

| **Field** | **Description** |
|-----------|------------------|
| `chrom`   | Chromosome number (normalized format without the "chr" prefix) |
| `pos`     | Position on the chromosome |
| `ref`     | Reference allele |
| `alt`     | Alternate allele(s) |
| `qual`    | Quality score for the variant |
| `filter`  | Filter status (e.g., "PASS" or specific filter criteria) |

### INFO Fields from VCF Files

Additional information extracted from VCF files includes:

| **INFO Field** | **Description** |
|----------------|------------------|
| `AC`           | Allele Count in genotypes |
| `AF`           | Allele Frequency |
| `AN`           | Total number of alleles in called genotypes |
| `ExcessHet`    | Excess heterozygosity statistic |
| `FS`           | Strand bias statistic |
| `MLEAC`        | Maximum likelihood estimate of allele count |
| `MLEAF`        | Maximum likelihood estimate of allele frequency |
| `MQ`           | RMS mapping quality |
| `QD`           | Quality by depth |
| `SOR`          | Symmetry odds ratio |
| `RS`           | Reference SNP ID (if available) |
| `source`       | Indicates whether the variant is from the VCF or ClinVar source |

### ClinVar Annotations for Chromosome 17

For variants on chromosome 17, additional ClinVar-specific fields are included:

| **ClinVar Field**     | **Description** |
|-----------------------|------------------|
| `clinvar_id`          | ClinVar Allele ID |
| `clinical_significance` | Clinical significance (e.g., "Pathogenic", "Benign") |
| `conditions`          | Associated conditions or disease names |
| `review_status`       | Review status of ClinVar submission |
| `AF_EXAC`             | Allele Frequency in ExAC |
| `CLNDISDB`            | ClinVar-disease database references |
| `CLNDN`               | Disease name associated with the variant |
| `CLNHGVS`             | Human Genome Variation Society (HGVS) nomenclature |
| `CLNVC`               | Variant type (e.g., SNV, InDel) |
| `CLNVCSO`             | Sequence ontology term for variant |
| `GENEINFO`            | Gene information |
| `MC`                  | Molecular consequence |
| `ORIGIN`              | Origin of the variant |
| `RS`                  | Reference SNP ID (mirrors the `RS` field in VCF data) |

---

## Genome Browser

The **Genome Browser** is a web-based interface developed with Flask that enables users to seamlessly query, filter, and visualize variant data from the genomic database. Key features include:

- **Filtering Options**: Users can filter variants based on chromosome, position, clinical significance, and other relevant fields.
- **Interactive Visualization**: View variants in a clean, organized layout with options to sort and navigate through data.
- **User-Friendly Interface**: Designed for ease of use, eliminating the need for in-depth knowledge of the underlying database structure.

---

## Setup Instructions

Follow these steps to set up and run the Genomic Database and Genome Browser on your local machine.

### Prerequisites

- **Python 3.7 or higher**: Ensure Python is installed on your system. You can download it from [python.org](https://www.python.org/downloads/).
- **pip**: Python package installer, typically included with Python installations.

### Installation

1. **Create a Virtual Environment (Optional)**:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows, use venv\Scripts\activate
   ```

2. **Install Required Dependencies**:
   Install the necessary Python packages listed in the `requirements.txt` file:
   ```bash
   pip install -r requirements.txt
   ```

### Running the Application

1. **Start the Flask Application**:
   ```bash
   python app.py
   ```

2. **Access the Genome Browser**:
   Open your web browser and navigate to [http://127.0.0.1:5000](http://127.0.0.1:5000) to access the genome browser interface.

---

## Example Queries

For users who prefer interacting directly with the TinyDB database, here are some example queries:

### Retrieve All Variants on Chromosome 17

```python
from tinydb import TinyDB, Query

# Initialize the database
db = TinyDB('genomic_db.json')
Variant = Query()

# Search for all variants on chromosome 17
results = db.search(Variant.chrom == '17')

# Example: Print the number of variants found
print(f"Total variants on chromosome 17: {len(results)}")
```

### Find ClinVar Variants Marked as "Pathogenic"

```python
from tinydb import TinyDB, Query

db = TinyDB('genomic_db.json')
Variant = Query()

# Search for pathogenic ClinVar variants
results = db.search((Variant.source == 'clinvar') & (Variant.clinical_significance == 'Pathogenic'))

print(f"Pathogenic ClinVar variants found: {len(results)}")
```

### List Variants with a Specific Allele Frequency Range

```python
from tinydb import TinyDB, Query

db = TinyDB('genomic_db.json')
Variant = Query()

# Define allele frequency range
min_af = 0.05
max_af = 0.10

# Search for variants within the specified allele frequency range
results = db.search((Variant.AF >= min_af) & (Variant.AF <= max_af))

print(f"Variants with AF between {min_af} and {max_af}: {len(results)}")
```

---

## Using the Genome Browser

The Genome Browser offers an intuitive interface to explore the genomic data without direct interaction with the database. Here's how to utilize its features:

1. **Access the Interface**:
   Navigate to [http://127.0.0.1:5000](http://127.0.0.1:5000) in your web browser.

2. **Search and Filter**:
   - **Chromosome**: Select the chromosome of interest (e.g., 17).
   - **Position**: Specify a genomic position or range.
   - **Clinical Significance**: Filter variants based on their clinical impact (e.g., "Pathogenic", "Benign").
   - **Additional Filters**: Utilize other fields such as allele frequency, gene information, etc., to narrow down results.

3. **View Results**:
   - **Table View**: Examine variant details in a structured table.
   - **Visualization**: Use integrated charts or plots (if available) to visualize data distributions and patterns.

4. **Export Data**:
   - **Download Options**: Export query results in formats like CSV for further analysis.
---

**Note**: Ensure that file paths within scripts are correctly specified to avoid any runtime errors. Adjust configurations as needed based on your system environment.

---
