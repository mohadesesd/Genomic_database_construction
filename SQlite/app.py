# app.py
from flask import Flask, render_template, request, abort
import sqlite3
import json
import copy

app = Flask(__name__)

DATABASE = '/home/mohadese/Desktop/Task2/SQlite/genomic_variants.db' # Ensure the filename and path are correct

def dict_factory(cursor, row):
    """Convert database row objects to a dictionary keyed by column name."""
    return {col[0]: row[idx] for idx, col in enumerate(cursor.description)}

def get_db_connection():
    """Establish a connection to the SQLite database."""
    conn = sqlite3.connect(DATABASE)
    conn.row_factory = dict_factory  # Use dict_factory to get dictionaries
    return conn

def get_filterable_columns():
    """Return a list of all filterable columns with their data types."""
    # Define columns with their data types: 'text' or 'numeric'
    return [
        {'name': 'chrom', 'type': 'text'},
        {'name': 'start_pos', 'type': 'numeric'},
        {'name': 'end_pos', 'type': 'numeric'},
        {'name': 'pos', 'type': 'numeric'},
        {'name': 'ref', 'type': 'text'},
        {'name': 'alt', 'type': 'text'},
        {'name': 'qual', 'type': 'numeric'},
        {'name': 'filter', 'type': 'text'},
        {'name': 'DP', 'type': 'numeric'},
        {'name': 'AF', 'type': 'numeric'},
        {'name': 'AC', 'type': 'numeric'},
        {'name': 'AN', 'type': 'numeric'},
        {'name': 'ExcessHet', 'type': 'numeric'},
        {'name': 'FS', 'type': 'numeric'},
        {'name': 'MLEAC', 'type': 'numeric'},
        {'name': 'MLEAF', 'type': 'numeric'},
        {'name': 'MQ', 'type': 'numeric'},
        {'name': 'QD', 'type': 'numeric'},
        {'name': 'SOR', 'type': 'numeric'},
        {'name': 'RS', 'type': 'numeric'},
        {'name': 'clinvar_id', 'type': 'text'},
        {'name': 'clinical_significance', 'type': 'text'},
        {'name': 'condition', 'type': 'text'},
        {'name': 'review_status', 'type': 'text'},
        {'name': 'CLNREVSTAT', 'type': 'text'},
        {'name': 'CLNSIG', 'type': 'text'},
        {'name': 'CLNVC', 'type': 'text'},
        {'name': 'CLNVCSO', 'type': 'text'},
        {'name': 'GENEINFO', 'type': 'text'},
        {'name': 'MC', 'type': 'text'},
        {'name': 'ORIGIN', 'type': 'text'},
        {'name': 'ALLELEID', 'type': 'numeric'},
        {'name': 'CLNDISDB', 'type': 'text'},
        {'name': 'CLNDN', 'type': 'text'},
        {'name': 'CLNHGVS', 'type': 'text'},
        {'name': 'AF_EXAC', 'type': 'numeric'}
    ]

def build_where_clause(filters, logic):
    """
    Build the WHERE clause for the SQL query based on the provided filters.

    Args:
        filters (list of dict): Each dict contains 'field', 'operator', 'value'.
        logic (str): 'and' or 'or' to combine the filters.

    Returns:
        tuple: (where_clause as string, params as list)
    """
    where_clauses = []
    params = []
    for filter in filters:
        field = filter.get('field')
        operator = filter.get('operator')
        value = filter.get('value')

        if not field or not operator or value is None:
            continue  # Skip incomplete filters

        # Determine the type of the field
        column = next((col for col in get_filterable_columns() if col['name'] == field), None)
        if not column:
            continue  # Skip unknown fields

        if column['type'] == 'numeric':
            # Handle numeric operators
            try:
                if operator == 'equals':
                    where_clauses.append(f"{field} = ?")
                    params.append(float(value))
                elif operator == 'greater_than':
                    where_clauses.append(f"{field} > ?")
                    params.append(float(value))
                elif operator == 'less_than':
                    where_clauses.append(f"{field} < ?")
                    params.append(float(value))
                elif operator == 'greater_than_or_equal':
                    where_clauses.append(f"{field} >= ?")
                    params.append(float(value))
                elif operator == 'less_than_or_equal':
                    where_clauses.append(f"{field} <= ?")
                    params.append(float(value))
                elif operator == 'between':
                    # Value should be in the format "min-max"
                    min_val, max_val = value.split('-', 1)
                    where_clauses.append(f"{field} BETWEEN ? AND ?")
                    params.extend([float(min_val.strip()), float(max_val.strip())])
            except ValueError:
                continue  # Skip filters with invalid numeric values
        else:
            # Handle text operators
            if operator == 'equals':
                where_clauses.append(f"{field} = ?")
                params.append(value)
            elif operator == 'contains':
                where_clauses.append(f"{field} LIKE ?")
                params.append(f"%{value}%")
            elif operator == 'starts_with':
                where_clauses.append(f"{field} LIKE ?")
                params.append(f"{value}%")
            elif operator == 'ends_with':
                where_clauses.append(f"{field} LIKE ?")
                params.append(f"%{value}")

    if where_clauses:
        where_clause = " WHERE " + f" {logic.upper()} ".join(where_clauses)
    else:
        where_clause = ""

    return where_clause, params

@app.route('/', methods=['GET'], endpoint='variants')
def index():
    # Get all filterable columns
    filterable_columns = get_filterable_columns()

    # Retrieve advanced search criteria from the form
    criteria = []
    num_criteria = int(request.args.get('num_criteria', 0))
    logic = request.args.get('logic', 'and').lower()
    for i in range(1, num_criteria + 1):
        field = request.args.get(f'field_{i}')
        operator = request.args.get(f'operator_{i}')
        value = request.args.get(f'value_{i}')
        if field and operator and value:
            criteria.append({'field': field, 'operator': operator, 'value': value})

    # Pagination parameters
    page = request.args.get('page', 1, type=int)
    per_page = 20
    offset = (page - 1) * per_page

    # Build WHERE clause
    where_clause, params = build_where_clause(criteria, logic)

    # Base SQL query with aliases to prevent duplication
    base_query = """
        SELECT 
            variants.variant_id AS variant_variant_id, 
            variants.chrom, variants.pos, variants.ref, variants.alt,
            variants.qual, variants.filter, variants.DP, variants.AF, variants.AC, variants.AN,
            variants.ExcessHet, variants.FS, variants.MLEAC, variants.MLEAF, variants.MQ,
            variants.QD, variants.SOR, variants.RS, variants.ANN,
            clinvar_annotations.clinvar_id, clinvar_annotations.clinical_significance,
            clinvar_annotations.condition, clinvar_annotations.review_status,
            clinvar_annotations.CLNREVSTAT, clinvar_annotations.CLNSIG,
            clinvar_annotations.CLNVC, clinvar_annotations.CLNVCSO,
            clinvar_annotations.GENEINFO, clinvar_annotations.MC,
            clinvar_annotations.ORIGIN, clinvar_annotations.ALLELEID,
            clinvar_annotations.CLNDISDB, clinvar_annotations.CLNDN,
            clinvar_annotations.CLNHGVS, clinvar_annotations.AF_EXAC
        FROM variants
        LEFT JOIN clinvar_annotations ON variants.variant_id = clinvar_annotations.variant_id
    """

    # Final query for fetching variants with pagination
    final_query = base_query + where_clause + " ORDER BY chrom, pos LIMIT ? OFFSET ?"
    query_params = params.copy()
    query_params.extend([per_page, offset])

    # Fetch variants from the database
    conn = get_db_connection()
    try:
        variants = conn.execute(final_query, query_params).fetchall()
    except Exception as e:
        conn.close()
        abort(500, description=f"Database query failed: {e}")

    # Process each variant to flatten 'ANN' fields (if necessary)
    processed_variants = []
    for variant in variants:
        variant_copy = copy.deepcopy(variant)  # Create a copy to modify

        # If 'ANN' field exists and is a JSON string, parse and flatten it
        if 'ANN' in variant_copy and variant_copy['ANN']:
            try:
                ann_data = json.loads(variant_copy['ANN'])
                for ann_key, ann_value in ann_data.items():
                    new_key = f"ann_{ann_key}"
                    variant_copy[new_key] = ann_value
            except json.JSONDecodeError:
                # Handle JSON parsing errors if 'ANN' is not a valid JSON
                variant_copy['ann'] = 'N/A'
        else:
            # If 'ANN' is missing, you can set default values or leave as is
            pass  # Optionally set default values

        # Remove the original 'ANN' field as it's now flattened
        if 'ANN' in variant_copy:
            del variant_copy['ANN']

        # Append the processed variant to the list
        processed_variants.append(variant_copy)

    # Update header_keys based on the first processed variant
    if processed_variants:
        header_keys = processed_variants[0].keys()
    else:
        header_keys = []

    # Pagination logic
    count_query = "SELECT COUNT(*) AS total_variants FROM variants LEFT JOIN clinvar_annotations ON variants.variant_id = clinvar_annotations.variant_id" + where_clause
    try:
        total_variants_result = conn.execute(count_query, params).fetchone()
        total_variants = total_variants_result['total_variants'] if total_variants_result else 0
    except Exception as e:
        conn.close()
        abort(500, description=f"Database count query failed: {e}")
    conn.close()
    total_pages = max(1, (total_variants + per_page - 1) // per_page)

    # Prepare filtered_args by removing 'page' from query parameters
    filtered_args = request.args.to_dict(flat=False)
    filtered_args.pop('page', None)

    return render_template(
        'index.html',
        variants=processed_variants,
        page=page,
        total_pages=total_pages,
        criteria=criteria or [],
        header_keys=header_keys,
        filterable_columns=filterable_columns,
        filtered_args=filtered_args
    )

@app.route('/variant/<int:variant_id>')
def variant_detail(variant_id):
    conn = get_db_connection()
    variant = conn.execute("""
        SELECT 
            variants.variant_id AS variant_variant_id, 
            variants.chrom, variants.pos, variants.ref, variants.alt,
            variants.qual, variants.filter, variants.DP, variants.AF, variants.AC, variants.AN,
            variants.ExcessHet, variants.FS, variants.MLEAC, variants.MLEAF, variants.MQ,
            variants.QD, variants.SOR, variants.RS, variants.ANN,
            clinvar_annotations.clinvar_id, clinvar_annotations.clinical_significance,
            clinvar_annotations.condition, clinvar_annotations.review_status,
            clinvar_annotations.CLNREVSTAT, clinvar_annotations.CLNSIG,
            clinvar_annotations.CLNVC, clinvar_annotations.CLNVCSO,
            clinvar_annotations.GENEINFO, clinvar_annotations.MC,
            clinvar_annotations.ORIGIN, clinvar_annotations.ALLELEID,
            clinvar_annotations.CLNDISDB, clinvar_annotations.CLNDN,
            clinvar_annotations.CLNHGVS, clinvar_annotations.AF_EXAC
        FROM variants
        LEFT JOIN clinvar_annotations ON variants.variant_id = clinvar_annotations.variant_id
        WHERE variants.variant_id = ?
    """, (variant_id,)).fetchone()
    conn.close()
    if variant is None:
        abort(404, description="Variant not found")

    # Process the 'ANN' field if needed
    variant_copy = copy.deepcopy(variant)
    if 'ANN' in variant_copy and variant_copy['ANN']:
        try:
            ann_data = json.loads(variant_copy['ANN'])
            for ann_key, ann_value in ann_data.items():
                new_key = f"ann_{ann_key}"
                variant_copy[new_key] = ann_value
        except json.JSONDecodeError:
            variant_copy['ann'] = 'N/A'
    else:
        pass  # Optionally set default values

    if 'ANN' in variant_copy:
        del variant_copy['ANN']

    return render_template('variant_detail.html', variant=variant_copy)

@app.errorhandler(404)
def page_not_found(e):
    """Custom 404 error page."""
    return render_template('404.html'), 404

@app.errorhandler(500)
def internal_error(e):
    """Custom 500 error page."""
    return render_template('500.html', error=e), 500

if __name__ == '__main__':
    app.run(debug=True)