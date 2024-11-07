# app.py

import logging
from logging.handlers import RotatingFileHandler
from flask import Flask, render_template, request, abort, Response
from tinydb import TinyDB, Query
from flask_paginate import Pagination, get_page_parameter
import os
import sys
from datetime import datetime
import math
import csv

app = Flask(__name__)

# ---------------------------- Configuration ---------------------------- #

# Update the DB_PATH to point to the correct genomic_db.json file
DB_PATH = '/home/mohadese/Desktop/Task2/TinyDB/genomic_db.json'  # Ensure the filename and path are correct
LOG_FILE = 'flask_app.log'

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
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
console_formatter = logging.Formatter('%(levelname)s: %(message)s')
console_handler.setFormatter(console_formatter)

# Add handlers to logger
logger.addHandler(file_handler)
logger.addHandler(console_handler)

# ---------------------------- Database Setup ---------------------------- #

if not os.path.exists(DB_PATH):
    logger.error(f"Database file '{DB_PATH}' not found. Please run 'process_vcf.py' first.")
    sys.exit(1)

db = TinyDB(DB_PATH)
Variant = Query()

# ---------------------------- Helper Functions ---------------------------- #

def get_filterable_columns():
    """
    Define which columns can be filtered.

    Returns:
        list: List of dictionaries with 'name' and 'type' keys.
    """
    return [
        {'name': 'chrom', 'type': 'text'},
        {'name': 'pos', 'type': 'number'},
        {'name': 'ref', 'type': 'text'},
        {'name': 'alt', 'type': 'text'},
        {'name': 'qual', 'type': 'number'},
        {'name': 'filter', 'type': 'text'},
        {'name': 'AC', 'type': 'number'},
        {'name': 'AF', 'type': 'number'},
        {'name': 'AN', 'type': 'number'},
        {'name': 'ExcessHet', 'type': 'number'},
        {'name': 'FS', 'type': 'number'},
        {'name': 'MLEAC', 'type': 'number'},
        {'name': 'MLEAF', 'type': 'number'},
        {'name': 'MQ', 'type': 'number'},
        {'name': 'QD', 'type': 'number'},
        {'name': 'SOR', 'type': 'number'},
        {'name': 'RS', 'type': 'number'},
        # ClinVar Fields
        {'name': 'clinvar_id', 'type': 'text'},
        {'name': 'clinical_significance', 'type': 'text'},
        {'name': 'conditions', 'type': 'text'},
        {'name': 'review_status', 'type': 'text'},
        {'name': 'AF_EXAC', 'type': 'number'},
        {'name': 'CLNDISDB', 'type': 'text'},
        {'name': 'CLNDN', 'type': 'text'},
        {'name': 'CLNHGVS', 'type': 'text'},
        {'name': 'CLNVC', 'type': 'text'},
        {'name': 'CLNVCSO', 'type': 'text'},
        {'name': 'GENEINFO', 'type': 'text'},
        {'name': 'MC', 'type': 'text'},
        {'name': 'ORIGIN', 'type': 'text'},
        # Add more fields as necessary
    ]

def build_query(search_criteria, logic='and'):
    """
    Build a TinyDB query based on search criteria.

    Args:
        search_criteria (list): List of dictionaries with 'field', 'operator', and 'value'.
        logic (str): 'and' or 'or' to combine conditions.

    Returns:
        Query: TinyDB Query object or None.
    """
    if not search_criteria:
        logger.debug("No search criteria provided.")
        return None  # No filtering

    queries = []
    for idx, criterion in enumerate(search_criteria, start=1):
        field = criterion.get('field')
        operator = criterion.get('operator')
        value = criterion.get('value')

        logger.debug(f"Criterion {idx}: field={field}, operator={operator}, value={value}")

        if not field or not operator or value is None:
            logger.warning(f"Criterion {idx} is incomplete; skipping.")
            continue

        field = field.strip()
        operator = operator.strip()
        value = value.strip()

        # Determine field type
        field_type = next((col['type'] for col in get_filterable_columns() if col['name'] == field), 'text')
        logger.debug(f"Field '{field}' is of type '{field_type}'.")

        if field_type == 'number':
            try:
                numeric_value = float(value)
                logger.debug(f"Converted value '{value}' to float: {numeric_value}")
            except ValueError:
                logger.warning(f"Invalid numeric value for field '{field}': '{value}'. Skipping criterion.")
                continue

            if operator == 'equals':
                queries.append(Variant[field] == numeric_value)
                logger.debug(f"Added query: Variant['{field}'] == {numeric_value}")
            elif operator == 'greater_than':
                queries.append(Variant[field] > numeric_value)
                logger.debug(f"Added query: Variant['{field}'] > {numeric_value}")
            elif operator == 'less_than':
                queries.append(Variant[field] < numeric_value)
                logger.debug(f"Added query: Variant['{field}'] < {numeric_value}")
            elif operator == 'greater_than_or_equal':
                queries.append(Variant[field] >= numeric_value)
                logger.debug(f"Added query: Variant['{field}'] >= {numeric_value}")
            elif operator == 'less_than_or_equal':
                queries.append(Variant[field] <= numeric_value)
                logger.debug(f"Added query: Variant['{field}'] <= {numeric_value}")
            else:
                logger.warning(f"Unsupported operator '{operator}' for numeric field '{field}'.")
        else:  # Text fields
            if operator == 'equals':
                queries.append(Variant[field] == value)
                logger.debug(f"Added query: Variant['{field}'] == '{value}'")
            elif operator == 'contains':
                queries.append(Variant[field].test(lambda v: value.lower() in str(v).lower()))
                logger.debug(f"Added query: Variant['{field}'].contains('{value}')")
            else:
                logger.warning(f"Unsupported operator '{operator}' for text field '{field}'.")

    if not queries:
        logger.debug("No valid queries were built from the search criteria.")
        return None

    # Combine queries based on logic
    if logic == 'and':
        query = queries.pop()
        logger.debug("Combining queries with AND logic.")
        for q in queries:
            query &= q
            logger.debug(f"Combined query: {query}")
    else:
        query = queries.pop()
        logger.debug("Combining queries with OR logic.")
        for q in queries:
            query |= q
            logger.debug(f"Combined query: {query}")

    logger.debug(f"Final built query: {query}")
    return query

# ---------------------------- Routes ---------------------------- #

@app.route('/')
def home():
    return render_template('home.html', current_year=datetime.now().year)

@app.route('/variants', methods=['GET'])
def variants():
    # Pagination parameters
    page = request.args.get(get_page_parameter(), type=int, default=1)
    per_page = 20

    # Advanced search parameters
    fields = request.args.getlist('field[]')
    operators = request.args.getlist('operator[]')
    values = request.args.getlist('value[]')
    logic = request.args.get('logic', 'and').lower()

    search_criteria = []
    for field, operator, value in zip(fields, operators, values):
        if field and operator and value:
            search_criteria.append({'field': field, 'operator': operator, 'value': value})

    logger.debug(f"Received search criteria: {search_criteria}")
    logger.debug(f"Combine logic: {logic}")

    # Build and execute query
    query = build_query(search_criteria, logic)
    if query:
        all_variants = db.search(query)
        logger.debug(f"Applied search criteria. {len(all_variants)} variants found.")
    else:
        all_variants = db.all()
        logger.debug(f"No search criteria applied. Retrieved {len(all_variants)} variants.")

    # Total variants after filtering
    total = len(all_variants)
    logger.debug(f"Total variants after filtering: {total}")

    # Calculate total pages
    total_pages = math.ceil(total / per_page) if total > 0 else 1
    logger.debug(f"Calculated total_pages: {total_pages}")

    # Ensure current page is within range
    if page < 1:
        page = 1
    elif page > total_pages:
        page = total_pages

    # Pagination slicing
    start = (page - 1) * per_page
    end = start + per_page
    variants_paginated = all_variants[start:end]
    logger.debug(f"Start index: {start}, End index: {end}")
    logger.debug(f"Variants paginated count: {len(variants_paginated)}")

    # Prepare pagination
    pagination = Pagination(
        page=page,
        total=total,
        per_page=per_page,
        css_framework='bootstrap5',
        record_name='variants',
        bs_version=5,
        link_label='<< Previous',
        link_next='Next >>',
        # Preserve current search parameters
        query_args=request.args.to_dict(flat=False)
    )

    # Define headers based on filterable columns
    filterable_columns = get_filterable_columns()
    header_keys = [col['name'] for col in filterable_columns]

    # Log pagination details
    logger.debug(f"Pagination - Page: {page}, Per Page: {per_page}, Total: {total}, Total Pages: {pagination.total_pages}")

    return render_template('variants.html',
                           variants=variants_paginated,
                           page=page,
                           per_page=per_page,
                           pagination=pagination,
                           total=total,
                           search_criteria=search_criteria,
                           filterable_columns=filterable_columns,
                           header_keys=header_keys,
                           current_year=datetime.now().year)

@app.route('/export', methods=['GET'])
def export_variants():
    """
    Export search results as a CSV file.

    Returns:
        Response: CSV file download.
    """
    # Extract search criteria from request
    fields = request.args.getlist('field[]')
    operators = request.args.getlist('operator[]')
    values = request.args.getlist('value[]')
    logic = request.args.get('logic', 'and').lower()

    search_criteria = []
    for field, operator, value in zip(fields, operators, values):
        if field and operator and value:
            search_criteria.append({'field': field, 'operator': operator, 'value': value})

    logger.debug(f"Export - Received search criteria: {search_criteria}")
    logger.debug(f"Export - Combine logic: {logic}")

    # Build query
    query = build_query(search_criteria, logic)
    if query:
        variants = db.search(query)
        logger.debug(f"Export - Applied search criteria. {len(variants)} variants found.")
    else:
        variants = db.all()
        logger.debug(f"Export - No search criteria applied. Retrieved {len(variants)} variants.")

    # Define CSV headers
    filterable_columns = get_filterable_columns()
    headers = [col['name'] for col in filterable_columns]

    # Create CSV
    def generate():
        yield ','.join(headers) + '\n'
        for variant in variants:
            row = [str(variant.get(header, '')) for header in headers]
            # Escape quotes and commas in fields
            escaped_row = ['"{}"'.format(field.replace('"', '""')) if ',' in field or '"' in field else field for field in row]
            yield ','.join(escaped_row) + '\n'

    return Response(generate(), mimetype='text/csv',
                    headers={"Content-Disposition": "attachment;filename=variants.csv"})

# ---------------------------- Debugging Routes ---------------------------- #

@app.route('/debug')
def debug_route():
    """
    Debug route to display sample variants from the database.
    """
    sample_variants = db.all()[:5]  # Retrieve first 5 variants
    return render_template('debug.html', variants=sample_variants, current_year=datetime.now().year)

@app.route('/test-query')
def test_query():
    """
    Test route to perform a sample query and display results.
    Example: Fetch variants with AF > 0.05
    """
    sample_variants = db.search(Variant.AF > 0.05)
    return render_template('test_query.html', variants=sample_variants, current_year=datetime.now().year)

@app.route('/debug-search', methods=['GET'])
def debug_search():
    """
    Debug route to display received search parameters and built queries.
    """
    # Extract search criteria from request
    fields = request.args.getlist('field[]')
    operators = request.args.getlist('operator[]')
    values = request.args.getlist('value[]')
    logic = request.args.get('logic', 'and').lower()

    search_criteria = []
    for field, operator, value in zip(fields, operators, values):
        if field and operator and value:
            search_criteria.append({'field': field, 'operator': operator, 'value': value})

    logger.debug(f"Debug Search - Received search criteria: {search_criteria}")
    logger.debug(f"Debug Search - Combine logic: {logic}")

    # Build and execute query
    query = build_query(search_criteria, logic)
    if query:
        all_variants = db.search(query)
        logger.debug(f"Debug Search - Applied search criteria. {len(all_variants)} variants found.")
    else:
        all_variants = db.all()
        logger.debug(f"Debug Search - No search criteria applied. Retrieved {len(all_variants)} variants.")

    # Total variants after filtering
    total = len(all_variants)

    # Pagination parameters
    page = request.args.get(get_page_parameter(), type=int, default=1)
    per_page = 20
    start = (page - 1) * per_page
    end = start + per_page
    variants_paginated = all_variants[start:end]

    # Prepare pagination
    pagination = Pagination(
        page=page,
        total=total,
        per_page=per_page,
        css_framework='bootstrap5',
        record_name='variants',
        bs_version=5,
        link_label='<< Previous',
        link_next='Next >>',
        # Preserve current search parameters
        query_args=request.args.to_dict(flat=False)
    )

    return render_template('debug_search.html',
                           search_criteria=search_criteria,
                           query=str(query),
                           total=total,
                           variants=variants_paginated,
                           page=page,
                           per_page=per_page,
                           pagination=pagination,
                           current_year=datetime.now().year)

# ---------------------------- Error Handlers ---------------------------- #

@app.errorhandler(404)
def page_not_found(e):
    return render_template('404.html', current_year=datetime.now().year), 404

@app.errorhandler(500)
def internal_server_error(e):
    return render_template('500.html', current_year=datetime.now().year), 500

# ---------------------------- Entry Point ---------------------------- #

if __name__ == '__main__':
    app.run(debug=True)
