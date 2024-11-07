import sqlite3
import json
from cyvcf2 import VCF
import os
import sys
import logging
import numpy as np


# ---------------------------- Configuration ---------------------------- #

# Directory containing your VCF and VCF.GZ files (excluding ClinVar)
VCF_DIRECTORY = '/home/mohadese/Downloads/Genomic Database Task files'  # Update with your actual path

# Path to the ClinVar VCF.GZ file
CLINVAR_VCF_PATH = '/home/mohadese/Downloads/clinvar.vcf.gz'  # Update with your actual path

# SQLite database file
DATABASE_PATH = 'genomic_variants.db'

# Log file path
LOG_FILE = 'insert_vcfs.log'

# ---------------------------- Logging Setup ---------------------------- #

logging.basicConfig(
    filename=LOG_FILE,
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.DEBUG
)

# ---------------------------- Helper Functions ---------------------------- #

def serialize_info(info):
    """
    Convert INFO field to a JSON-serializable dictionary.
    """
    def convert(value):
        if isinstance(value, np.generic):
            return value.item()
        elif isinstance(value, (list, tuple)):
            return [convert(v) for v in value]
        elif isinstance(value, dict):
            return {k: convert(v) for k, v in value.items()}
        else:
            return value

    return {k: convert(v) for k, v in info.items()}

def normalize_chrom(chrom):
    """
    Normalize chromosome names by removing 'chr' prefix if present and converting to uppercase.
    """
    chrom = chrom.strip()
    if chrom.lower().startswith('chr'):
        chrom = chrom[3:]
    if chrom.upper() == 'M':
        return 'MT'
    return chrom.upper()

def connect_db(db_path=DATABASE_PATH):
    """
    Connect to the SQLite database.
    """
    try:
        conn = sqlite3.connect(db_path)
        conn.execute("PRAGMA foreign_keys = ON;")
        logging.info(f"Connected to SQLite database at {db_path}.")
        return conn
    except sqlite3.Error as e:
        logging.error(f"Error connecting to database: {e}")
        sys.exit(1)

def initialize_database(conn):
    """
    Initialize the database by creating necessary tables and indexes.
    """
    cursor = conn.cursor()
    try:
        cursor.executescript("""
        DROP TABLE IF EXISTS variants;
        DROP TABLE IF EXISTS samples;
        DROP TABLE IF EXISTS genotype;
        DROP TABLE IF EXISTS clinvar_annotations;

        CREATE TABLE IF NOT EXISTS variants (
            variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            qual REAL,
            filter TEXT,
            info TEXT,
            DP INTEGER,
            AF REAL,
            AC INTEGER,
            AN INTEGER,
            ExcessHet REAL,
            FS REAL,
            MLEAC INTEGER,
            MLEAF REAL,
            MQ REAL,
            QD REAL,
            SOR REAL,
            ANN TEXT,
            RS INTEGER,
            UNIQUE(chrom, pos, ref, alt)
        );

        CREATE TABLE IF NOT EXISTS samples (
            sample_id INTEGER PRIMARY KEY AUTOINCREMENT,
            sample_name TEXT UNIQUE NOT NULL
        );

        CREATE TABLE IF NOT EXISTS genotype (
            genotype_id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id INTEGER NOT NULL,
            sample_id INTEGER NOT NULL,
            genotype TEXT,
            FOREIGN KEY (variant_id) REFERENCES variants(variant_id),
            FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
            UNIQUE(variant_id, sample_id)
        );

        CREATE TABLE IF NOT EXISTS clinvar_annotations (
            annotation_id INTEGER PRIMARY KEY AUTOINCREMENT,
            variant_id INTEGER NOT NULL,
            clinvar_id TEXT,
            clinical_significance TEXT,
            condition TEXT,
            review_status TEXT,
            CLNREVSTAT TEXT,
            CLNSIG TEXT,
            CLNVC TEXT,
            CLNVCSO TEXT,
            GENEINFO TEXT,
            MC TEXT,
            ORIGIN TEXT,
            ALLELEID INTEGER,
            CLNDISDB TEXT,
            CLNDN TEXT,
            CLNHGVS TEXT,
            AF_EXAC REAL,
            FOREIGN KEY (variant_id) REFERENCES variants(variant_id),
            UNIQUE(variant_id)
        );

        -- Create indexes
        CREATE INDEX IF NOT EXISTS idx_variants_chrom_pos ON variants (chrom, pos);
        CREATE INDEX IF NOT EXISTS idx_variants_chrom_pos_ref_alt ON variants (chrom, pos, ref, alt);
        CREATE INDEX IF NOT EXISTS idx_clinvar_variant_id ON clinvar_annotations (variant_id);
        CREATE INDEX IF NOT EXISTS idx_genotype_variant_id ON genotype (variant_id);
        CREATE INDEX IF NOT EXISTS idx_genotype_sample_id ON genotype (sample_id);
        """)
        conn.commit()
        logging.info("Database initialized successfully with required tables and indexes.")
    except sqlite3.Error as e:
        logging.error(f"Error initializing database: {e}")
        conn.rollback()
        sys.exit(1)
    finally:
        cursor.close()

def parse_ann_field(ann_field):
    """
    Parse the ANN field into a list of dictionaries.
    """
    annotations = []
    if ann_field:
        if isinstance(ann_field, str):
            ann_entries = [ann_field]
        else:
            ann_entries = ann_field

        for ann_entry in ann_entries:
            ann_parts = ann_entry.split('|')
            ann_keys = [
                'Allele', 'Consequence', 'Impact', 'Symbol',
                'Gene', 'Feature_type', 'Feature', 'Biotype',
                'Rank', 'HGVS.c', 'HGVS.p', 'cDNA_position', 'CDS_position',
                'Protein_position', 'Distance', 'ERRORS', 'WARNINGS'
            ]
            ann_dict = dict(zip(ann_keys, ann_parts))
            annotations.append(ann_dict)
    return annotations

def insert_variant(cursor, variant):
    """
    Insert a variant into the variants table with normalized chromosome names and serialized INFO field.
    Also extract specific INFO fields into separate columns.
    """
    chrom = normalize_chrom(variant.CHROM)
    pos = variant.POS
    ref = variant.REF.strip()
    alt_list = [allele.strip() for allele in variant.ALT] if variant.ALT else ['.']
    qual = float(variant.QUAL) if variant.QUAL not in ('.', None) else None
    filter_status = ';'.join(variant.FILTER) if variant.FILTER else 'PASS'

    # Serialize INFO field
    try:
        info_dict = dict(variant.INFO)
        info_serialized = serialize_info(info_dict)
        info_json = json.dumps(info_serialized)
    except Exception as e:
        logging.error(f"Failed to serialize INFO field for variant {chrom}:{pos}:{ref}>{','.join(alt_list)}: {e}", exc_info=True)
        info_serialized = {}
        info_json = json.dumps(info_serialized)

    # Extract and handle specific INFO fields
    def handle_field(value, field_name, expected_type):
        if value is None:
            return None
        if isinstance(value, list):
            logging.debug(f"Field '{field_name}' is a list: {value}. Using first element.")
            value = value[0]
        try:
            return expected_type(value)
        except ValueError:
            logging.warning(f"Could not convert field '{field_name}' value '{value}' to {expected_type}. Setting as None.")
            return None

    AC = handle_field(info_serialized.get('AC'), 'AC', int)
    AF = handle_field(info_serialized.get('AF'), 'AF', float)
    AN = handle_field(info_serialized.get('AN'), 'AN', int)
    DP = handle_field(info_serialized.get('DP'), 'DP', int)
    ExcessHet = handle_field(info_serialized.get('ExcessHet'), 'ExcessHet', float)
    FS = handle_field(info_serialized.get('FS'), 'FS', float)
    MLEAC = handle_field(info_serialized.get('MLEAC'), 'MLEAC', int)
    MLEAF = handle_field(info_serialized.get('MLEAF'), 'MLEAF', float)
    MQ = handle_field(info_serialized.get('MQ'), 'MQ', float)
    QD = handle_field(info_serialized.get('QD'), 'QD', float)
    SOR = handle_field(info_serialized.get('SOR'), 'SOR', float)
    RS = handle_field(info_serialized.get('RS'), 'RS', int)

    ANN_raw = info_serialized.get('ANN')
    ANN_json = None
    if ANN_raw:
        ANN_parsed = parse_ann_field(ANN_raw)
        ANN_json = json.dumps(ANN_parsed)

    variant_id = None

    try:
        for alt in alt_list:
            cursor.execute("""
                INSERT INTO variants (
                    chrom, pos, ref, alt, qual, filter, info, DP, AF, AC, AN,
                    ExcessHet, FS, MLEAC, MLEAF, MQ, QD, SOR, ANN, RS
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                chrom, pos, ref, alt, qual, filter_status, info_json, DP, AF, AC, AN,
                ExcessHet, FS, MLEAC, MLEAF, MQ, QD, SOR, ANN_json, RS
            ))
            variant_id = cursor.lastrowid
    except sqlite3.IntegrityError:
        # Variant already exists
        cursor.execute("""
            SELECT variant_id FROM variants
            WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?
        """, (chrom, pos, ref, alt))
        result = cursor.fetchone()
        variant_id = result[0] if result else None
    except Exception as e:
        logging.error(f"Unexpected error inserting variant {chrom}:{pos}:{ref}>{alt}: {e}", exc_info=True)
        raise

    return variant_id

def insert_sample(cursor, sample_name, sample_ids):
    """
    Insert a sample into the samples table.
    """
    normalized_sample = sample_name.strip()
    if normalized_sample in sample_ids:
        return sample_ids[normalized_sample]
    try:
        cursor.execute("""
            INSERT INTO samples (sample_name)
            VALUES (?)
        """, (normalized_sample,))
        sample_id = cursor.lastrowid
        sample_ids[normalized_sample] = sample_id
    except sqlite3.IntegrityError:
        # Sample already exists
        cursor.execute("""
            SELECT sample_id FROM samples
            WHERE sample_name = ?
        """, (normalized_sample,))
        result = cursor.fetchone()
        sample_id = result[0] if result else None
        if sample_id:
            sample_ids[normalized_sample] = sample_id
        else:
            logging.error(f"Failed to retrieve sample ID for '{normalized_sample}'.")
    except Exception as e:
        logging.error(f"Unexpected error inserting sample '{normalized_sample}': {e}", exc_info=True)
        raise
    return sample_id

def insert_genotype(cursor, variant_id, sample_id, genotype):
    """
    Insert a genotype into the genotype table.
    """
    try:
        cursor.execute("""
            INSERT INTO genotype (variant_id, sample_id, genotype)
            VALUES (?, ?, ?)
        """, (variant_id, sample_id, genotype))
    except sqlite3.IntegrityError:
        # Genotype already exists
        pass
    except Exception as e:
        logging.error(f"Unexpected error inserting genotype for variant ID {variant_id}, sample ID {sample_id}: {e}", exc_info=True)
        raise

def insert_clinvar_annotation(cursor, variant_id, clinvar_info):
    """
    Insert a ClinVar annotation into the clinvar_annotations table.
    """
    def get_value(value):
        if isinstance(value, list):
            return value[0] if value else None
        elif value is not None:
            return value
        else:
            return None

    clinvar_id = get_value(clinvar_info.get('RCV'))
    clinical_significance = get_value(clinvar_info.get('CLNSIG'))
    condition = get_value(clinvar_info.get('CLNDBN'))
    review_status = get_value(clinvar_info.get('CLNREVSTAT'))
    CLNREVSTAT = get_value(clinvar_info.get('CLNREVSTAT'))
    CLNSIG = get_value(clinvar_info.get('CLNSIG'))
    CLNVC = get_value(clinvar_info.get('CLNVC'))
    CLNVCSO = get_value(clinvar_info.get('CLNVCSO'))
    GENEINFO = get_value(clinvar_info.get('GENEINFO'))
    MC = get_value(clinvar_info.get('MC'))
    ORIGIN = get_value(clinvar_info.get('ORIGIN'))
    ALLELEID = get_value(clinvar_info.get('ALLELEID'))
    CLNDISDB = get_value(clinvar_info.get('CLNDISDB'))
    CLNDN = get_value(clinvar_info.get('CLNDN'))
    CLNHGVS = get_value(clinvar_info.get('CLNHGVS'))
    AF_EXAC = get_value(clinvar_info.get('AF_EXAC'))

    # Convert lists to strings if necessary
    def convert_to_string(value):
        if isinstance(value, list):
            return ','.join(map(str, value))
        return str(value) if value is not None else None

    clinvar_id = convert_to_string(clinvar_id)
    clinical_significance = convert_to_string(clinical_significance)
    condition = convert_to_string(condition)
    review_status = convert_to_string(review_status)
    CLNREVSTAT = convert_to_string(CLNREVSTAT)
    CLNSIG = convert_to_string(CLNSIG)
    CLNVC = convert_to_string(CLNVC)
    CLNVCSO = convert_to_string(CLNVCSO)
    GENEINFO = convert_to_string(GENEINFO)
    MC = convert_to_string(MC)
    ORIGIN = convert_to_string(ORIGIN)
    CLNDISDB = convert_to_string(CLNDISDB)
    CLNDN = convert_to_string(CLNDN)
    CLNHGVS = convert_to_string(CLNHGVS)

    # Convert ALLELEID and AF_EXAC to appropriate types
    try:
        ALLELEID = int(ALLELEID) if ALLELEID is not None else None
    except ValueError:
        ALLELEID = None
    try:
        AF_EXAC = float(AF_EXAC) if AF_EXAC is not None else None
    except ValueError:
        AF_EXAC = None

    try:
        cursor.execute("""
            INSERT OR IGNORE INTO clinvar_annotations (
                variant_id, clinvar_id, clinical_significance, condition, review_status,
                CLNREVSTAT, CLNSIG, CLNVC, CLNVCSO, GENEINFO, MC, ORIGIN,
                ALLELEID, CLNDISDB, CLNDN, CLNHGVS, AF_EXAC
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            variant_id, clinvar_id, clinical_significance, condition, review_status,
            CLNREVSTAT, CLNSIG, CLNVC, CLNVCSO, GENEINFO, MC, ORIGIN,
            ALLELEID, CLNDISDB, CLNDN, CLNHGVS, AF_EXAC
        ))
    except Exception as e:
        logging.error(f"Error inserting ClinVar annotation for variant ID {variant_id}: {e}", exc_info=True)
        raise

def process_vcf(conn, vcf_path, sample_ids):
    """
    Process a single VCF or VCF.GZ file and insert its data into the database.
    """
    cursor = conn.cursor()
    try:
        conn.execute('BEGIN TRANSACTION')
        logging.info(f"Processing VCF file: {vcf_path}")
        vcf = VCF(vcf_path)

        # Insert samples first
        for sample in vcf.samples:
            insert_sample(cursor, sample, sample_ids)

        # Insert variants and genotypes
        for variant in vcf:
            variant_id = insert_variant(cursor, variant)
            if variant_id is None:
                continue  # Skip if variant ID couldn't be retrieved

            for sample_idx, sample in enumerate(vcf.samples):
                normalized_sample = sample.strip()
                if normalized_sample not in sample_ids:
                    continue  # Skip inserting genotype for this sample

                sample_id = sample_ids[normalized_sample]
                gt = variant.genotypes[sample_idx]
                if gt[0] == -1 or gt[1] == -1:
                    genotype = './.'  # Missing genotype
                else:
                    genotype = f"{gt[0]}/{gt[1]}"
                insert_genotype(cursor, variant_id, sample_id, genotype)

        conn.commit()
        logging.info(f"Successfully processed VCF file: {vcf_path}")
    except Exception as e:
        conn.rollback()
        logging.error(f"Error processing VCF file {vcf_path}: {e}", exc_info=True)
    finally:
        cursor.close()

def process_clinvar_vcf(conn, clinvar_vcf_path):
    """
    Process the ClinVar VCF file and insert annotations into the database.
    """
    cursor = conn.cursor()
    unmatched_variants = []
    try:
        conn.execute('BEGIN TRANSACTION')
        logging.info(f"Processing ClinVar VCF file: {clinvar_vcf_path}")
        clinvar_vcf = VCF(clinvar_vcf_path)

        for variant in clinvar_vcf:
            chrom = normalize_chrom(variant.CHROM)
            pos = variant.POS
            ref = variant.REF.strip()
            alt_list = [allele.strip() for allele in variant.ALT] if variant.ALT else ['.']

            for alt in alt_list:
                # Fetch variant_id from variants table
                cursor.execute("""
                    SELECT variant_id FROM variants
                    WHERE chrom = ? AND pos = ? AND ref = ? AND alt = ?
                """, (chrom, pos, ref, alt))
                result = cursor.fetchone()
                if result:
                    variant_id = result[0]
                    insert_clinvar_annotation(cursor, variant_id, variant.INFO)
                else:
                    unmatched_variants.append(f"{chrom}:{pos}:{ref}>{alt}")

        conn.commit()
        logging.info(f"Successfully processed ClinVar VCF file: {clinvar_vcf_path}")
    except Exception as e:
        conn.rollback()
        logging.error(f"Error processing ClinVar VCF file {clinvar_vcf_path}: {e}", exc_info=True)
    finally:
        cursor.close()

    # Write unmatched variants to a separate log file
    if unmatched_variants:
        with open('unmatched_variants.log', 'a') as f:
            for variant in unmatched_variants:
                f.write(f"{variant}\n")

# ---------------------------- Main Execution ---------------------------- #

def main():
    conn = connect_db()
    initialize_database(conn)

    # Process your VCF files
    vcf_files = [
        os.path.join(VCF_DIRECTORY, f)
        for f in os.listdir(VCF_DIRECTORY)
        if (f.endswith('.vcf') or f.endswith('.vcf.gz')) and 'clinvar' not in f.lower()
    ]

    sample_ids = {}

    for vcf_path in vcf_files:
        if os.path.isfile(vcf_path):
            process_vcf(conn, vcf_path, sample_ids)
        else:
            logging.warning(f"File not found: {vcf_path}")

    # Process the ClinVar VCF file
    if os.path.isfile(CLINVAR_VCF_PATH):
        process_clinvar_vcf(conn, CLINVAR_VCF_PATH)
    else:
        logging.error(f"ClinVar VCF file not found: {CLINVAR_VCF_PATH}")

    conn.close()
    logging.info("Database processing complete.")

if __name__ == "__main__":
    main()
