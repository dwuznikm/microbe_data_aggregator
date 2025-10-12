import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import webbrowser
from api_client import get_genome_summary, fetch_ncbi_taxonomy
import csv


class GenomeApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Genome Data Explorer")
        self.geometry("1100x650")
        self.create_widgets()
        self.genome_data = []

    def create_widgets(self):
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(expand=True, fill="both")

        # --- Search Tab ---
        self.search_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.search_frame, text="Search")

        # TaxID Entry
        tk.Label(self.search_frame, text="Enter Taxonomy ID:").pack(pady=5)
        self.taxid_entry = ttk.Entry(self.search_frame, width=30)
        self.taxid_entry.pack()

        # Database checkboxes
        db_frame = ttk.Frame(self.search_frame)
        db_frame.pack(pady=5)
        self.ncbi_var = tk.BooleanVar(value=True)
        self.ensembl_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(db_frame, text="NCBI", variable=self.ncbi_var).pack(side="left", padx=5)
        ttk.Checkbutton(db_frame, text="Ensembl", variable=self.ensembl_var).pack(side="left", padx=5)

        # Progress label
        self.progress_label = ttk.Label(self.search_frame, text="")
        self.progress_label.pack(pady=5)

        # Search button
        search_button = ttk.Button(self.search_frame, text="Search", command=self.search_taxid)
        search_button.pack(pady=5)

        # --- Taxonomy Tab ---
        self.taxonomy_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.taxonomy_frame, text="Taxonomy")

        self.tax_tree = ttk.Treeview(self.taxonomy_frame)
        self.tax_tree.pack(expand=True, fill="both", padx=10, pady=10)

        # --- Genome Tab ---
        self.genome_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.genome_frame, text="Genome Info")

        # --- Filters ---
        filter_frame = ttk.Frame(self.genome_frame)
        filter_frame.pack(fill="x", padx=10, pady=5)

        self.ref_var = tk.BooleanVar()
        ttk.Checkbutton(filter_frame, text="Reference genomes only", variable=self.ref_var).grid(row=0, column=0, padx=5)
        
        self.remove_dups_var = tk.BooleanVar()
        ttk.Checkbutton(filter_frame, text="Remove duplicates", variable=self.remove_dups_var).grid(row=0, column=8, padx=5)

        tk.Label(filter_frame, text="Seq Length min:").grid(row=0, column=1, padx=5)
        self.seq_min = ttk.Entry(filter_frame, width=10)
        self.seq_min.grid(row=0, column=2, padx=5)

        tk.Label(filter_frame, text="Seq Length max:").grid(row=0, column=3, padx=5)
        self.seq_max = ttk.Entry(filter_frame, width=10)
        self.seq_max.grid(row=0, column=4, padx=5)

        tk.Label(filter_frame, text="Limit rows:").grid(row=0, column=5, padx=5)
        self.limit_rows = ttk.Entry(filter_frame, width=10)
        self.limit_rows.grid(row=0, column=6, padx=5)
        export_btn = ttk.Button(filter_frame, text="Export to CSV", command=self.export_to_csv)
        export_btn.grid(row=0, column=9, padx=10)
        apply_btn = ttk.Button(filter_frame, text="Apply Filters", command=self.apply_filters)
        apply_btn.grid(row=0, column=7, padx=10)

        # --- Genome Table ---
        columns = ("Accession", "Assembly Level", "Seq Length", "GC Content", "# Genes", "Source", "Reference", "Link", "Duplicate")
        self.genome_table = ttk.Treeview(self.genome_frame, columns=columns, show="headings")
        
        for col in columns:
            self.genome_table.heading(col, text=col, command=lambda c=col: self.sort_treeview_column(self.genome_table, c, False))
            self.genome_table.column(col, width=120, anchor="w")
        self.genome_table.pack(expand=True, fill="both", padx=10, pady=10)

        # Bind double-click to open link
        self.genome_table.bind("<Double-1>", self.open_genome_link)

    # ---------------- Genome Link ----------------
    def open_genome_link(self, event):
        selected_item = self.genome_table.selection()
        if not selected_item:
            return
        row = self.genome_table.item(selected_item)
        col = self.genome_table.identify_column(event.x)
        col_index = int(col.replace("#", "")) - 1
        values = row["values"]

        # Open Link column
        if col_index == 7:
            url = values[col_index]
            if url:
                webbrowser.open(url)

    # ---------------- Search ----------------
    def search_taxid(self):
        taxid_str = self.taxid_entry.get().strip()
        if not taxid_str.isdigit():
            messagebox.showerror("Invalid Input", "Please enter a valid numeric Tax ID.")
            return

        taxid = int(taxid_str)

        # Determine which databases to fetch
        selected_sources = []
        if self.ncbi_var.get():
            selected_sources.append("NCBI")
        if self.ensembl_var.get():
            selected_sources.append("Ensembl")
        if not selected_sources:
            messagebox.showerror("No Database Selected", "Select at least one database to search.")
            return

        # Clear previous results
        self.tax_tree.delete(*self.tax_tree.get_children())
        for row in self.genome_table.get_children():
            self.genome_table.delete(row)
        self.genome_data.clear()

        # --- Show progress ---
        self.progress_label.config(text=f"Fetching genomes from {', '.join(selected_sources)}...")
        self.update_idletasks()

        # --- Fetch genome summary ---
        summary = get_genome_summary(taxid, sources=selected_sources)
        if not summary:
            messagebox.showerror("Error", f"Failed to fetch data for Tax ID {taxid}")
            self.progress_label.config(text="")
            return

        self.genome_data = summary.get("genomes", [])

        # Populate genome table
        self.populate_genome_table(self.genome_data)

        # --- Fetch Taxonomy Tree ---
        self.progress_label.config(text="Fetching taxonomy data...")
        self.update_idletasks()
        tax_data = fetch_ncbi_taxonomy(taxid)
        if not tax_data:
            messagebox.showwarning("Warning", "Failed to fetch taxonomy data")
            self.progress_label.config(text="")
            return
        self.build_taxonomy_tree(tax_data)

        # Update label to done
        self.progress_label.config(text="Done")

    # ---------------- Populate Genome Table ----------------
    def populate_genome_table(self, genomes):
        
        for row in self.genome_table.get_children():
            self.genome_table.delete(row)

        accession_counts = {}
        for g in genomes:
            acc = g.get("Accession")
            if acc:
                accession_counts[acc] = accession_counts.get(acc, 0) + 1

        for genome in genomes:
            acc = genome.get("Accession")
            duplicate_flag = "Yes" if acc and accession_counts.get(acc, 0) > 1 else "No"
            genome["Duplicate"] = duplicate_flag

        if self.remove_dups_var.get():
            seen = set()
            filtered_genomes = []
            for g in genomes:
                acc = g.get("Accession")
                if acc:
                    if acc not in seen:
                        filtered_genomes.append(g)
                        seen.add(acc)
                else:
                    filtered_genomes.append(g)
            genomes = filtered_genomes


        for genome in genomes:
            values = (
                genome.get("Accession", ""),
                genome.get("Assembly Level", ""),
                genome.get("Seq Length", ""),
                genome.get("GC Content", ""),
                genome.get("# Genes", ""),
                genome.get("Source", ""),
                genome.get("Reference", ""),
                genome.get("Link", ""),
                genome.get("Duplicate", "")
            )
            self.genome_table.insert("", "end", values=values)

        self.adjust_column_widths(self.genome_table)

    # ---------------- Apply Filters ----------------
    def apply_filters(self):
        filtered = self.genome_data

        # Reference genome filter
        if self.ref_var.get():
            filtered = [g for g in filtered if g.get("Reference") == "Yes"]

        # Sequence length min/max
        try:
            seq_min = int(self.seq_min.get()) if self.seq_min.get() else None
            seq_max = int(self.seq_max.get()) if self.seq_max.get() else None
        except ValueError:
            messagebox.showerror("Invalid Input", "Sequence length filters must be numeric")
            return

        if seq_min is not None:
            filtered = [g for g in filtered if g.get("Seq Length") and int(g.get("Seq Length")) >= seq_min]
        if seq_max is not None:
            filtered = [g for g in filtered if g.get("Seq Length") and int(g.get("Seq Length")) <= seq_max]

        # Limit rows
        try:
            limit = int(self.limit_rows.get()) if self.limit_rows.get() else None
        except ValueError:
            messagebox.showerror("Invalid Input", "Limit must be numeric")
            return

        if limit is not None:
            filtered = filtered[:limit]

        self.populate_genome_table(filtered)

    def sort_treeview_column(self, tree, col, reverse=False):
        """
        Sort the Treeview column when header is clicked.
        tree: ttk.Treeview instance
        col: column name (string)
        reverse: bool, sort descending if True
        """
        data_list = [(tree.set(k, col), k) for k in tree.get_children('')]

        try:
            data_list.sort(key=lambda t: float(t[0]) if t[0] else float('-inf'), reverse=reverse)
        except ValueError:
            data_list.sort(key=lambda t: t[0], reverse=reverse)

        for index, (val, k) in enumerate(data_list):
            tree.move(k, '', index)

        tree.heading(col, command=lambda: self.sort_treeview_column(tree, col, not reverse))

    def export_to_csv(self):
        file_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")],
            title="Save genome table as CSV"
        )
        if not file_path:
            return

        columns = self.genome_table["columns"]

        try:
            with open(file_path, mode="w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(columns)  # header
                for row_id in self.genome_table.get_children():
                    row = self.genome_table.item(row_id)["values"]
                    writer.writerow(row)
            messagebox.showinfo("Export Successful", f"Genome table saved to {file_path}")
        except Exception as e:
            messagebox.showerror("Export Failed", f"Could not save CSV:\n{e}")


    # ---------------- Adjust Column Widths ----------------
    @staticmethod
    def adjust_column_widths(tree):
        max_widths = []
        for col in tree["columns"]:
            max_width = max(
                [len(str(tree.set(item, col))) for item in tree.get_children()] + [len(col)]
            )
            pixel_width = min(max_width * 7 + 20, 500)  # cap at 500px
            max_widths.append(pixel_width)
        for i, col in enumerate(tree["columns"]):
            tree.column(col, width=max_widths[i])

    # ---------------- Taxonomy Tree ----------------
    def build_taxonomy_tree(self, tax_data):
        try:
            self.tax_tree.delete(*self.tax_tree.get_children())
            tax = tax_data["reports"][0]["taxonomy"]
            classification = tax.get("classification", {})

            ranks = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
            indent = 0

            for rank in ranks:
                node = classification.get(rank)
                if node:
                    name = node.get("name", "Unknown")
                    tid = node.get("id", "")
                    display_text = " " * (indent * 4) + f"{rank.capitalize()}: {name} (TaxID: {tid})"
                    item_id = self.tax_tree.insert("", "end", text=display_text)
                    indent += 1

            species_name = classification.get("species", {}).get("name")
            current_name = tax.get("current_scientific_name", {}).get("name")

            # Only add strain info if current_name != species_name (since strain is not included in the api update)
            if species_name and current_name and current_name != species_name:
                display_text = " " * (indent * 4) + f"Strain/variant: {current_name}"
                self.tax_tree.insert("", "end", text=display_text)

            # Expand all items
            for item in self.tax_tree.get_children():
                self.tax_tree.item(item, open=True)

        except Exception:
            self.tax_tree.insert("", "end", text="Could not parse taxonomy tree")



if __name__ == "__main__":
    app = GenomeApp()
    app.mainloop()
