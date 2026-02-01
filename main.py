import tkinter as tk
from tkinter import ttk, messagebox, simpledialog, filedialog
import webbrowser
import threading

import api_client


class GenomeApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Genome Data Explorer")
        self.geometry("1100x650")

        # --- session state ---
        self.user_email = None
        self.email_provided = False
        self.max_ncbi_genome_count = None
        self.max_ncbi_gene_count = None
        self.genome_data = []

        # --- build UI ---
        self.create_widgets()

    def create_widgets(self):
        # --- Tabs ---
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(expand=True, fill="both")

        # --- Search tab ---
        self.search_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.search_frame, text="Search")

        ttk.Label(self.search_frame, text="Enter Taxonomy ID:").pack(pady=(12, 6))
        self.taxid_entry = ttk.Entry(self.search_frame, width=30)
        self.taxid_entry.pack()

        sources_frame = ttk.Frame(self.search_frame)
        sources_frame.pack(pady=10)
        self.ncbi_var = tk.BooleanVar(value=True)
        self.ensembl_var = tk.BooleanVar(value=True)
        self.ena_var = tk.BooleanVar(value=True)
        self.bvbrc_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(sources_frame, text="NCBI", variable=self.ncbi_var).pack(
            side="left", padx=6
        )
        ttk.Checkbutton(sources_frame, text="Ensembl", variable=self.ensembl_var).pack(
            side="left", padx=6
        )
        ttk.Checkbutton(sources_frame, text="ENA", variable=self.ena_var).pack(
            side="left", padx=6
        )
        ttk.Checkbutton(sources_frame, text="BV-BRC", variable=self.bvbrc_var).pack(
            side="left", padx=6
        )

        # Search button
        self.search_btn = ttk.Button(
            self.search_frame, text="Search", command=self.on_search_button
        )
        self.search_btn.pack(pady=16)

        # --- Progress frame ---
        self.progress_frame = ttk.Frame(self.search_frame)
        self.progress_label = ttk.Label(self.progress_frame, text="")
        self.progress_label.pack()
        self.progress_bar = ttk.Progressbar(
            self.progress_frame, orient="horizontal", length=400, mode="determinate"
        )
        self.progress_percent = ttk.Label(self.progress_frame, text="0%")
        self.progress_bar.pack(pady=2)
        self.progress_percent.pack(pady=2)
        self.progress_frame.pack(pady=10)
        self.progress_frame.pack_forget()  # hide initially
        self.progress_text_only = ttk.Label(
            self.search_frame, text="", foreground="blue"
        )
        self.progress_text_only.pack(pady=5)
        self.progress_text_only.pack_forget()  # hide initially

        # --- Taxonomy tab ---
        self.taxonomy_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.taxonomy_frame, text="Taxonomy")
        self.tax_tree = ttk.Treeview(self.taxonomy_frame)
        self.tax_tree.pack(expand=True, fill="both", padx=10, pady=10)

        # --- Genome tab ---
        self.genome_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.genome_frame, text="Genome Info")
        self._build_genome_table()
        self._build_genome_filters()

        # --- Gene tab ---
        self.gene_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.gene_frame, text="Gene Info")
        self._build_gene_table()

    # ---------- UI Table Builders ----------
    def _build_genome_table(self):
        columns = (
            "Accession",
            "Assembly Level",
            "Seq Length",
            "GC Content",
            "# Genes",
            "Source",
            "Reference",
            "Link",
            "Duplicate",
        )
        self.genome_table = ttk.Treeview(
            self.genome_frame, columns=columns, show="headings"
        )
        for col in columns:
            self.genome_table.heading(
                col,
                text=col,
                command=lambda c=col: self.sort_treeview_column(
                    self.genome_table, c, False
                ),
            )
            self.genome_table.column(col, width=120, anchor="w")
        self.genome_table.pack(expand=True, fill="both", padx=10, pady=10)
        self.genome_table.bind("<Double-1>", self.open_genome_link)

    def _build_genome_filters(self):
        filter_frame = ttk.Frame(self.genome_frame)
        filter_frame.pack(fill="x", padx=10, pady=5)
        self.ref_var = tk.BooleanVar()
        self.remove_dups_var = tk.BooleanVar()
        ttk.Checkbutton(
            filter_frame, text="Reference genomes only", variable=self.ref_var
        ).grid(row=0, column=0, padx=5)
        ttk.Checkbutton(
            filter_frame, text="Remove duplicates", variable=self.remove_dups_var
        ).grid(row=0, column=1, padx=5)
        ttk.Label(filter_frame, text="Seq Length min:").grid(row=0, column=2, padx=5)
        self.seq_min = ttk.Entry(filter_frame, width=10)
        self.seq_min.grid(row=0, column=3, padx=5)
        ttk.Label(filter_frame, text="Seq Length max:").grid(row=0, column=4, padx=5)
        self.seq_max = ttk.Entry(filter_frame, width=10)
        self.seq_max.grid(row=0, column=5, padx=5)
        ttk.Label(filter_frame, text="Limit rows:").grid(row=0, column=6, padx=5)
        self.limit_rows = ttk.Entry(filter_frame, width=10)
        self.limit_rows.grid(row=0, column=7, padx=5)
        ttk.Button(filter_frame, text="Apply Filters", command=self.apply_filters).grid(
            row=0, column=8, padx=10
        )
        ttk.Button(filter_frame, text="Export to CSV", command=self.export_to_csv).grid(
            row=0, column=9, padx=10
        )

    def _build_gene_table(self):
        columns = (
            "Gene ID",
            "Type",
            "Description",
            "Symbol",
            "Locus",
            "Assemblies",
            "Link",
            "Duplicate",
        )
        self.gene_table = ttk.Treeview(
            self.gene_frame, columns=columns, show="headings"
        )
        for col in columns:
            self.gene_table.heading(
                col,
                text=col,
                command=lambda c=col: self.sort_treeview_column(
                    self.gene_table, c, False
                ),
            )
            self.gene_table.column(col, width=150, anchor="w")
        self.gene_table.pack(expand=True, fill="both", padx=10, pady=10)
        self.gene_table.bind("<Double-1>", self.on_gene_double_click)

    # ---------- Initial Search ----------
    def on_search_button(self):
        taxid_str = self.taxid_entry.get().strip()
        if not taxid_str.isdigit():
            messagebox.showerror("Invalid Input", "Please enter a numeric Taxonomy ID.")
            return
        taxid = int(taxid_str)

        sources = []
        if self.ncbi_var.get():
            sources.append("NCBI")
        if self.ensembl_var.get():
            sources.append("Ensembl")
        if self.ena_var.get():
            sources.append("ENA")
        if self.bvbrc_var.get():
            sources.append("BV-BRC")
        if not sources:
            messagebox.showerror(
                "No Database Selected", "Select at least one database."
            )
            return

        email = None
        if "NCBI" in sources and not self.email_provided:
            email = simpledialog.askstring(
                "NCBI Email", "Enter your email for NCBI Entrez:", parent=self
            )
            if not email:
                return
            self.user_email = email
            self.email_provided = True
        elif "NCBI" in sources:
            email = self.user_email

        # --- Use proper popup ---
        if "NCBI" in sources:
            threading.Thread(
                target=self.fetch_counts_thread,
                args=(taxid, email, sources),
                daemon=True,
            ).start()
        else:
            # Show text-only progress
            self.progress_text_only.config(text="Starting search...")
            self.progress_text_only.pack()
            self.update_idletasks()

            threading.Thread(
                target=self.full_search_thread,
                args=(taxid, None, None, sources, False),
                daemon=True,
            ).start()

    # ---------- NCBI Max Count Popup ----------
    def fetch_counts_thread(self, taxid, email, sources):
        try:
            genome_count, gene_count = api_client.get_total_ncbi_genome_count(
                taxid, email
            )
            self.max_ncbi_genome_count = genome_count
            self.max_ncbi_gene_count = gene_count
            self.after(
                0, lambda: self.show_max_popup(genome_count, gene_count, taxid, sources)
            )
        except Exception as e:
            self.after(
                0,
                lambda: messagebox.showerror(
                    "Count Failed", f"Could not fetch counts:\n{e}"
                ),
            )

    def show_max_popup(self, genome_count, gene_count, taxid, sources):
        popup = tk.Toplevel(self)
        popup.title("NCBI Max Records")
        popup.grab_set()

        ttk.Label(popup, text="Max genomic records:").grid(
            row=0, column=0, padx=6, pady=6, sticky="e"
        )
        max_genome_entry = ttk.Entry(popup, width=12)
        max_genome_entry.grid(row=0, column=1, sticky="w")
        max_genome_entry.insert(0, str(genome_count))
        ttk.Label(popup, text=f"Max: {genome_count}").grid(row=0, column=2, sticky="w")

        ttk.Label(popup, text="Max genetic records:").grid(
            row=1, column=0, padx=6, pady=6, sticky="e"
        )
        max_gene_entry = ttk.Entry(popup, width=12)
        max_gene_entry.grid(row=1, column=1, sticky="w")
        max_gene_entry.insert(0, str(gene_count))
        ttk.Label(popup, text=f"Max: {gene_count}").grid(row=1, column=2, sticky="w")

        def start_search():
            try:
                max_genomes_val = int(max_genome_entry.get())
                max_genes_val = int(max_gene_entry.get())
            except ValueError:
                messagebox.showerror("Invalid Input", "Max records must be numeric")
                return
            popup.destroy()
            self.on_start_search(taxid, sources, max_genomes_val, max_genes_val)

        ttk.Button(popup, text="Start Search", command=start_search).grid(
            row=2, column=0, columnspan=3, pady=10
        )

    def on_start_search(self, taxid, sources, max_genomes=None, max_genes=None):
        self.progress_frame.pack()
        self.progress_bar["value"] = 0
        self.progress_percent.config(text="0%")
        self.progress_label.config(text="Starting search...")
        self.update_idletasks()

        threading.Thread(
            target=self.full_search_thread,
            args=(taxid, max_genomes, max_genes, sources),
            daemon=True,
        ).start()

    def full_search_thread(self, taxid, max_genomes, max_genes, sources, show_bar=True):
        # Reset data at start of search
        self.genome_data = []

        def update_progress_label(text):
            if show_bar:
                self.after(0, lambda: self.progress_label.config(text=text))
            else:
                self.after(0, lambda: self.progress_text_only.config(text=text))

        def update_progress_bar(percent):
            if show_bar:
                self.after(0, lambda: self.progress_bar.config(value=percent))
                self.after(0, lambda: self.progress_percent.config(text=f"{percent}%"))

        # --- Taxonomy ---
        update_progress_label("Fetching taxonomy data...")
        try:
            tax_data = api_client.fetch_ncbi_taxonomy(taxid)
            self.build_taxonomy_tree(tax_data)
        except Exception:
            pass

        # --- Genomes ---
        for source in sources:
            update_progress_label(f"Fetching {source} genomic data...")

            def genome_progress(current, total):
                pct = int(current / total * 100) if total else 0
                update_progress_bar(pct)

            try:
                summary = api_client.get_genome_summary(
                    taxid,
                    max_genomes or 0,
                    sources=[source],
                    progress_callback=genome_progress,
                )
                genomes = summary.get("genomes", [])
                self.genome_data.extend(genomes)
                self.after(
                    0, lambda g=genomes: self.populate_genome_table(self.genome_data)
                )
            except Exception as e:
                self.after(
                    0,
                    lambda e=e: messagebox.showerror(
                        "Error", f"Failed to fetch genomic data from {source}:\n{e}"
                    ),
                )

        # --- Genes (NCBI only) ---
        if "NCBI" in sources:
            update_progress_label("Fetching NCBI genetic data...")

            def gene_progress(current, total):
                pct = int(current / total * 100) if total else 0
                update_progress_bar(pct)

            try:
                genes = api_client.get_gene_summary(
                    taxid, max_genes or 0, progress_callback=gene_progress
                )
                self.after(0, lambda: self.populate_gene_table(genes))
            except Exception as e:
                self.after(
                    0,
                    lambda e=e: messagebox.showerror(
                        "Error", f"Failed to fetch genetic data:\n{e}"
                    ),
                )

            update_progress_bar(100)

        update_progress_label("Done")
        if show_bar:
            self.after(2000, self.progress_frame.pack_forget)
        else:
            self.after(2000, self.progress_text_only.pack_forget)

    # ---------- Populate Tables ----------
    def populate_genome_table(self, genomes):
        for r in self.genome_table.get_children():
            self.genome_table.delete(r)
        accession_counts = {}
        for g in genomes:
            acc = g.get("Accession")
            if acc:
                accession_counts[acc] = accession_counts.get(acc, 0) + 1
        for g in genomes:
            g["Duplicate"] = (
                "Yes" if accession_counts.get(g.get("Accession"), 0) > 1 else "No"
            )
        if self.remove_dups_var.get():
            seen = set()
            genomes = [
                x
                for x in genomes
                if x.get("Accession") not in seen
                and (seen.add(x.get("Accession")) or True)
            ]
        for g in genomes:
            self.genome_table.insert(
                "",
                "end",
                values=(
                    g.get("Accession", ""),
                    g.get("Assembly Level", ""),
                    g.get("Seq Length", ""),
                    g.get("GC Content", ""),
                    g.get("# Genes", ""),
                    g.get("Source", ""),
                    g.get("Reference", ""),
                    g.get("Link", ""),
                    g.get("Duplicate", ""),
                ),
            )
        self.adjust_column_widths(self.genome_table)

    def populate_gene_table(self, genes):
        for r in self.gene_table.get_children():
            self.gene_table.delete(r)
        gene_id_counts = {}
        for g in genes:
            gid = g.get("Gene ID")
            if gid:
                gene_id_counts[gid] = gene_id_counts.get(gid, 0) + 1
        for g in genes:
            g["Duplicate"] = (
                "Yes" if gene_id_counts.get(g.get("Gene ID"), 0) > 1 else "No"
            )
            self.gene_table.insert(
                "",
                "end",
                values=(
                    g.get("Gene ID", ""),
                    g.get("Type", ""),
                    g.get("Description", ""),
                    g.get("Symbol", ""),
                    g.get("Locus", ""),
                    g.get("Assemblies", ""),
                    g.get("Link", ""),
                    g.get("Duplicate", ""),
                ),
            )
        self.adjust_column_widths(self.gene_table)

    # ---------- Filters ----------
    def apply_filters(self):
        filtered = self.genome_data
        if self.ref_var.get():
            filtered = [g for g in filtered if g.get("Reference") == "Yes"]
        try:
            seq_min = int(self.seq_min.get()) if self.seq_min.get() else None
            seq_max = int(self.seq_max.get()) if self.seq_max.get() else None
        except:
            messagebox.showerror("Invalid Input", "Seq filters must be numeric")
            return
        if seq_min:
            filtered = [
                g
                for g in filtered
                if g.get("Seq Length") and int(g.get("Seq Length")) >= seq_min
            ]
        if seq_max:
            filtered = [
                g
                for g in filtered
                if g.get("Seq Length") and int(g.get("Seq Length")) <= seq_max
            ]
        try:
            limit = int(self.limit_rows.get()) if self.limit_rows.get() else None
        except:
            messagebox.showerror("Invalid Input", "Limit must be numeric")
            return
        if limit:
            filtered = filtered[:limit]
        self.populate_genome_table(filtered)

    # ---------- Utilities ----------
    def sort_treeview_column(self, tree, col, reverse=False):
        data = [(tree.set(k, col), k) for k in tree.get_children()]
        try:
            data.sort(
                key=lambda t: float(t[0]) if t[0] else float("-inf"), reverse=reverse
            )
        except:
            data.sort(key=lambda t: t[0], reverse=reverse)
        for i, (v, k) in enumerate(data):
            tree.move(k, "", i)
        tree.heading(
            col, command=lambda: self.sort_treeview_column(tree, col, not reverse)
        )

    def export_to_csv(self):
        path = filedialog.asksaveasfilename(
            defaultextension=".csv", filetypes=[("CSV", "*.csv")]
        )
        if not path:
            return
        cols = self.genome_table["columns"]
        with open(path, "w", newline="", encoding="utf-8") as f:
            import csv

            w = csv.writer(f)
            w.writerow(cols)
            for r in self.genome_table.get_children():
                w.writerow(self.genome_table.item(r)["values"])
        messagebox.showinfo("Exported", f"Saved to {path}")

    def open_genome_link(self, event):
        sel = self.genome_table.selection()
        if not sel:
            return
        row = self.genome_table.item(sel[0])
        col = int(self.genome_table.identify_column(event.x).replace("#", "")) - 1
        url = row["values"][col]
        if url:
            webbrowser.open(url)

    def on_gene_double_click(self, event):
        sel = self.gene_table.selection()
        if not sel:
            return
        row = self.gene_table.item(sel[0])
        col = int(self.gene_table.identify_column(event.x).replace("#", "")) - 1
        url = row["values"][col]
        if url:
            webbrowser.open(url)

    def adjust_column_widths(self, tree):
        for col in tree["columns"]:
            max_len = max(
                [len(str(tree.set(item, col))) for item in tree.get_children()]
                + [len(col)]
            )
            tree.column(col, width=min(max_len * 6 + 20, 500))

    def build_taxonomy_tree(self, tax_data):
        try:
            self.tax_tree.delete(*self.tax_tree.get_children())
            tax = tax_data["reports"][0]["taxonomy"]
            cls = tax.get("classification", {})
            indent = 0
            ranks = [
                "domain",
                "kingdom",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
            ]
            for r in ranks:
                node = cls.get(r)
                if node:
                    self.tax_tree.insert(
                        "",
                        "end",
                        text=" " * indent * 4
                        + f"{r.capitalize()}: {node.get('name')} (TaxID:{node.get('id')})",
                    )
                    indent += 1
            species_name = cls.get("species", {}).get("name")
            current_name = tax.get("current_scientific_name", {}).get("name")
            if species_name and current_name and current_name != species_name:
                self.tax_tree.insert(
                    "", "end", text=" " * indent * 4 + f"Strain/variant: {current_name}"
                )
            for item in self.tax_tree.get_children():
                self.tax_tree.item(item, open=True)
        except:
            self.tax_tree.insert("", "end", text="Could not parse taxonomy tree")


if __name__ == "__main__":
    app = GenomeApp()
    app.mainloop()
