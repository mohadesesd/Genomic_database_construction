#!/usr/bin/env python3
import os
import logging
from tinydb import TinyDB, Query
from cyvcf2 import VCF
from tqdm import tqdm
from logging.handlers import RotatingFileHandler
import sys

# ---------------------------- Configuration ---------------------------- #

# Paths for VCF and ClinVar VCF files
VCF_DIRECTORY = '/home/mohadese/Downloads/Genomic Database Task files/'
CLINVAR_VCF_PATH = '/home/mohadese/Downloads/clinvar.vcf.gz'
DB_PATH = 'genomic_dab.json'
LOG_FILE = 'integration.log'
BATCH_SIZE = 1000

# ---------------------------- Logging Setup ---------------------------- #

# Initialize logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Create rotating file handler
file_handler = RotatingFileHandler(LOG_FILE, maxBytes=5 * 1024 * 1024, backupCount=5)
file_handler.setLevel(logging.DEBUG)
file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(file_formatter)

# Create console handler
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_formatter = logging.Formatter('%(levelname)s - %(message)s')
console_handler.setFormatter(console_formatter)

# Add handlers to logger
logger.addHandler(file_handler)
logger.addHandler(console_handler)

# ---------------------------- Helper Functions ---------------------------- #

def normalize_chrom(chrom):
    """
    Normalize chromosome names by removing 'chr' prefix and uppercasing.

    Args:
        chrom (str): Chromosome name (e.g., 'chr1', '1', 'chrX', 'x').

    Returns:
        str: Normalized chromosome name (e.g., '1', 'X').
    """
    if chrom.lower().startswith('chr'):
        return chrom[3:].upper()
    return chrom.upper()

def insert_records(db, records):
    """
    Insert records into TinyDB in batches.

    Args:
        db (TinyDB): TinyDB database instance.
        records (list): List of dictionaries representing variant records.
    """
    if records:
        try:
            db.insert_multiple(records)
            logger.debug(f"Inserted {len(records)} records.")
            records.clear()
        except Exception as e:
            logger.error(f"Failed to insert records: {e}")

# process_vcf.py

def get_int(info_fields, key, default=0):
    """
    Retrieve an integer value from info_fields.
    
    Args:
        info_fields (dict): Dictionary containing INFO fields from VCF.
        key (str): The key to retrieve.
        default (int): Default value if key is missing or conversion fails.
    
    Returns:
        int: The integer value of the field.
    """
    try:
        return int(info_fields.get(key, default))
    except (ValueError, TypeError):
        return default

def get_float(info_fields, key, default=0.0):
    """
    Retrieve a float value from info_fields.
    
    Args:
        info_fields (dict): Dictionary containing INFO fields from VCF.
        key (str): The key to retrieve.
        default (float): Default value if key is missing or conversion fails.
    
    Returns:
        float: The float value of the field.
    """
    try:
        return float(info_fields.get(key, default))
    except (ValueError, TypeError):
        return default


# ---------------------------- Processing Functions ---------------------------- #

def load_existing_variants(db):
    """
    Load existing variants from TinyDB into a set for quick lookup.

    Args:
        db (TinyDB): TinyDB database instance.

    Returns:
        set: Set of tuples representing existing variants (chrom, pos, ref, alt).
    """
    logger.info("Loading existing variants from the database into memory...")
    Variant = Query()
    existing = set()
    try:
        all_records = db.all()
        for record in all_records:
            key = (record.get('chrom'), record.get('pos'), record.get('ref'), record.get('alt'))
            existing.add(key)
        logger.info(f"Loaded {len(existing)} existing variants.")
    except Exception as e:
        logger.error(f"Error loading existing variants: {e}")
    return existing

def parse_vcf(db, vcf_directory, existing_variants):
    """
    Parse VCF files and insert variant data into TinyDB.

    Args:
        db (TinyDB): TinyDB database instance.
        vcf_directory (str): Path to directory containing VCF files.
        existing_variants (set): Set of existing variant keys to avoid duplication.
    """
    vcf_files = [os.path.join(vcf_directory, f) for f in os.listdir(vcf_directory) if f.endswith('.vcf.gz')]
    Variant = Query()

    for vcf_file in vcf_files:
        logger.info(f"Processing VCF file: {vcf_file}")
        vcf = None
        try:
            vcf = VCF(vcf_file)
            records = []
            for variant in tqdm(vcf, desc=f"Processing {os.path.basename(vcf_file)}", unit="variants"):
                chrom = normalize_chrom(variant.CHROM)
                pos, ref, alt_list = variant.POS, variant.REF, variant.ALT
                info_fields = dict(variant.INFO)
                filter_status = ";".join(variant.FILTER) if variant.FILTER else "PASS"

                for alt in alt_list:
                    key = (chrom, pos, ref, alt)
                    if key in existing_variants:
                        continue  # Skip existing variant

                    record = {
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "qual": variant.QUAL if variant.QUAL not in ('.', None) else None,
                        "filter": filter_status,
                        "AC": get_int(info_fields, "AC"),
                        "AF": get_float(info_fields, "AF"),
                        "AN": get_int(info_fields, "AN"),
                        "ExcessHet": get_float(info_fields, "ExcessHet"),
                        "FS": get_float(info_fields, "FS"),
                        "MLEAC": get_int(info_fields, "MLEAC"),
                        "MLEAF": get_float(info_fields, "MLEAF"),
                        "MQ": get_float(info_fields, "MQ"),
                        "QD": get_float(info_fields, "QD"),
                        "SOR": get_float(info_fields, "SOR"),
                        "RS": get_int(info_fields, "RS"),
                        "source": "vcf"  # Indicate origin
                    }

                    records.append(record)
                    existing_variants.add(key)  # Add to existing to prevent future duplicates

                    if len(records) >= BATCH_SIZE:
                        insert_records(db, records)

            # Insert any remaining records after processing the file
            insert_records(db, records)
            logger.info(f"Completed processing {vcf_file}")
        except Exception as e:
            logger.error(f"Error processing {vcf_file}: {e}", exc_info=True)
        finally:
            if vcf:
                vcf.close()

def parse_clinvar(db, clinvar_vcf, existing_variants):
    """
    Parse ClinVar VCF and update annotations in TinyDB for existing variants.
    Excludes inserting variants that are only present in ClinVar.

    Args:
        db (TinyDB): TinyDB database instance.
        clinvar_vcf (str): Path to ClinVar VCF file.
        existing_variants (set): Set of existing variant keys to identify updates.
    """
    logger.info("Processing ClinVar VCF for chromosome 17...")
    clinvar = None
    unmatched_variants = []
    Variant = Query()

    try:
        clinvar = VCF(clinvar_vcf)
        records_updated = 0
        for variant in tqdm(clinvar, desc="Processing ClinVar VCF", unit="variants"):
            chrom_normalized = normalize_chrom(variant.CHROM)
            if chrom_normalized != '17':
                continue  # Skip chromosomes other than 17

            pos, ref, alt_list = variant.POS, variant.REF, variant.ALT
            info_fields = variant.INFO

            for alt in alt_list:
                key = (chrom_normalized, pos, ref, alt)
                clinvar_data = {
                    "clinvar_id": info_fields.get("ALLELEID"),
                    "clinical_significance": info_fields.get("CLNSIG"),
                    "conditions": info_fields.get("CLNDBN"),
                    "review_status": info_fields.get("CLNREVSTAT"),
                    "AF_EXAC": info_fields.get("AF_EXAC"),
                    "CLNDISDB": info_fields.get("CLNDISDB"),
                    "CLNDN": info_fields.get("CLNDN"),
                    "CLNHGVS": info_fields.get("CLNHGVS"),
                    "CLNVC": info_fields.get("CLNVC"),
                    "CLNVCSO": info_fields.get("CLNVCSO"),
                    "GENEINFO": info_fields.get("GENEINFO"),
                    "MC": info_fields.get("MC"),
                    "ORIGIN": info_fields.get("ORIGIN"),
                    "RS": info_fields.get("RS"),
                }

                if key in existing_variants:
                    # Update existing record with ClinVar data
                    db.update(clinvar_data, (Variant.chrom == chrom_normalized) & 
                             (Variant.pos == pos) & 
                             (Variant.ref == ref) & 
                             (Variant.alt == alt))
                    records_updated += 1
                else:
                    # Do not insert new records for ClinVar-only variants
                    unmatched_variants.append(f"{chrom_normalized}:{pos}:{ref}>{alt}")

        clinvar.close()
        logger.info(f"Successfully processed ClinVar VCF file. Updated {records_updated} records.")

    except Exception as e:
        logger.error(f"Error processing ClinVar VCF file {clinvar_vcf}: {e}", exc_info=True)
    finally:
        # Write unmatched variants to a separate log file
        if unmatched_variants:
            with open('unmatched_variants.log', 'a') as f:
                for variant in unmatched_variants:
                    f.write(f"{variant}\n")
            logger.warning(f"Unmatched variants found: {len(unmatched_variants)}. Details in 'unmatched_variants.log'.")

# ---------------------------- Main Execution ---------------------------- #

def main():
    """
    Main function to orchestrate the integration of VCF files and ClinVar data into TinyDB.
    """
    try:
        logger.info("Starting integration process...")
        with TinyDB(DB_PATH) as db:
            # Load existing variants to minimize database searches
            existing_variants = load_existing_variants(db)
            
            # Parse and insert variants from VCF files
            parse_vcf(db, VCF_DIRECTORY, existing_variants)
            
            # Parse and integrate ClinVar annotations
            parse_clinvar(db, CLINVAR_VCF_PATH, existing_variants)
        
        logger.info("Integration process complete.")
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}", exc_info=True)

# ---------------------------- Entry Point ---------------------------- #

if __name__ == "__main__":
    main()
