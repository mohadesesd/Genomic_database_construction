import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import sqlite3
import logging
import csv

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    filename='app_debug.log',
    filemode='w'
)

# Database path
DATABASE = '/home/mohadese/Desktop/Task2/SQlite/genomic_variants.db' # Ensure the filename and path are correct

numeric_fields = ["pos", "qual", "DP", "AF", "AC", "AN", "ExcessHet", "FS", "MLEAC", "MLEAF", "MQ", "QD", "SOR", "RS"]

current_results = []  # To store query results for export

def get_db_connection():
    try:
        conn = sqlite3.connect(DATABASE)
        conn.row_factory = sqlite3.Row
        return conn
    except sqlite3.Error as e:
        messagebox.showerror("Database Connection Error", f"Failed to connect to database: {e}")
        logging.error("Database Connection Error: %s", e)
        return None

def query_db():
    conn = get_db_connection()
    if not conn:
        return
    cursor = conn.cursor()
    
    where_clauses = []
    params = []
    for row in criteria_rows:
        field = row['field'].get()
        operator = row['operator'].get()
        value = row['value'].get().strip()
        if field and operator and value:
            if field in numeric_fields:
                try:
                    value = int(value)
                except ValueError:
                    messagebox.showerror("Input Error", f"Invalid value for {field}. Please enter an integer.")
                    return
            if operator == 'equals':
                where_clauses.append(f"{field} = ?")
                params.append(value)
            elif operator == 'contains':
                where_clauses.append(f"{field} LIKE ?")
                params.append(f"%{value}%")
            elif operator == 'greater_than':
                where_clauses.append(f"{field} > ?")
                params.append(value)
            elif operator == 'less_than':
                where_clauses.append(f"{field} < ?")
                params.append(value)

    query = """
        SELECT variants.*, clinvar_annotations.*
        FROM variants
        LEFT JOIN clinvar_annotations ON variants.variant_id = clinvar_annotations.variant_id
    """
    if where_clauses:
        query += " WHERE " + " AND ".join(where_clauses)
    
    try:
        cursor.execute(query, params)
        rows = cursor.fetchall()
        conn.close()
        show_results(rows)
    except sqlite3.Error as e:
        conn.close()
        messagebox.showerror("Query Error", f"An SQLite error occurred: {e}")
        logging.error("SQLite error: %s", e)

def show_results(rows):
    for widget in results_frame.winfo_children():
        widget.destroy()
    if not rows:
        tk.Label(results_frame, text="No results found.", font=('Arial', 14)).pack(pady=20)
        return

    headers = rows[0].keys()
    tree = ttk.Treeview(results_frame, columns=headers, show='headings', style="Custom.Treeview")
    for header in headers:
        tree.heading(header, text=header)
        tree.column(header, anchor='center', width=100)
    for row in rows:
        values = [row[header] for header in headers]
        tree.insert('', 'end', values=values)

    scrollbar = ttk.Scrollbar(results_frame, orient="vertical", command=tree.yview)
    tree.configure(yscrollcommand=scrollbar.set)
    scrollbar.pack(side='right', fill='y')
    tree.pack(expand=True, fill='both')

    global current_results
    current_results = rows

def export_results():
    if not current_results:
        messagebox.showwarning("No Data", "No data to export. Please apply a filter first.")
        return
    file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
    if file_path:
        try:
            with open(file_path, 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(current_results[0].keys())
                for row in current_results:
                    writer.writerow(row)
            messagebox.showinfo("Export Successful", f"Results exported to {file_path}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Failed to export: {e}")

def add_criteria_row():
    row_num = len(criteria_rows)
    row = {}
    row['field'] = ttk.Combobox(filter_frame, values=filterable_columns, state="readonly", width=20)
    row['field'].grid(row=row_num, column=0, padx=5, pady=5)
    row['operator'] = ttk.Combobox(filter_frame, values=["equals", "contains", "greater_than", "less_than"], state="readonly", width=15)
    row['operator'].grid(row=row_num, column=1, padx=5, pady=5)
    row['value'] = ttk.Entry(filter_frame, width=20)
    row['value'].grid(row=row_num, column=2, padx=5, pady=5)
    criteria_rows.append(row)

# Tkinter GUI setup
app = tk.Tk()
app.title("Genomic Variant Database Browser")
app.geometry("1200x700")
app.configure(bg='#fafafa')

# Set up ttk styles for the Treeview and buttons
style = ttk.Style()
style.configure("TButton", background="#005f73", foreground="white", font=('Arial', 10, 'bold'), padding=6)
style.map("TButton",
          background=[('active', '#0a9396'), ('!disabled', '#005f73')],
          foreground=[('!disabled', 'white')])
style.configure("Custom.Treeview.Heading", font=('Arial', 10, 'bold'), background="#005f73", foreground="white")
style.configure("Custom.Treeview", background="#e0f7fa", fieldbackground="#e0f7fa")

# Title label
title_label = tk.Label(app, text="Genomic Variant Database Browser", font=('Helvetica', 16, 'bold'), fg="#0a9396", bg='#fafafa')
title_label.pack(pady=10)

# Filter Frame
filter_frame = ttk.Frame(app, padding="10", style="Custom.TFrame")
filter_frame.pack(padx=10, pady=10, fill="x")

criteria_rows = []
filterable_columns = ["chrom", "pos", "ref", "alt", "qual", "filter", "DP", "AF", "AC", "AN", "ExcessHet", "FS", "MLEAC", "MLEAF", "MQ", "QD", "SOR", "RS", "clinvar_id", "clinical_significance", "condition", "review_status", "CLNREVSTAT", "CLNSIG", "CLNVC", "CLNVCSO", "GENEINFO", "MC", "ORIGIN", "ALLELEID", "CLNDISDB", "CLNDN", "CLNHGVS", "AF_EXAC"]
add_criteria_row()

# Button frame for adding and applying filters
button_frame = ttk.Frame(filter_frame)
button_frame.grid(row=0, column=3, padx=10, pady=5)

# Create and style buttons
add_criteria_button = ttk.Button(button_frame, text="Add Criteria", command=add_criteria_row, style="TButton")
add_criteria_button.grid(row=0, column=0, padx=5)

apply_button = ttk.Button(button_frame, text="Apply Filter", command=query_db, style="TButton")
apply_button.grid(row=0, column=1, padx=5)

export_button = ttk.Button(button_frame, text="Export Results", command=export_results, style="TButton")
export_button.grid(row=0, column=2, padx=5)

# Results Frame with Scrollbar
results_container = ttk.Frame(app, style="Custom.TFrame")
results_container.pack(padx=10, pady=10, fill="both", expand=True)
results_frame = ttk.Frame(results_container)
results_frame.pack(fill="both", expand=True)

app.mainloop()
