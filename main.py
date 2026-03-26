import logging
import base64
import os
from io import BytesIO

import tkinter as tk
from tkinter import ttk, messagebox, simpledialog, filedialog
import webbrowser
import threading
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from collections import Counter, defaultdict
import numpy as np
from datetime import datetime

import api_client
import reports

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


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
        self.max_bvbrc_genome_count = None
        self.genome_data = []
        self.taxonomy_data = None

        # --- runtime state ---
        self._search_running = False

        # --- build UI ---
        self.create_widgets()

    def _build_summary_tab(self):
        self.summary_stats_frame = ttk.Frame(self.summary_frame)
        self.summary_stats_frame.pack(fill="x", padx=10, pady=10)

        self.summary_label = ttk.Label(
            self.summary_stats_frame, text="", justify="left"
        )
        self.summary_label.pack(anchor="w")

        self.export_html_btn = ttk.Button(
            self.summary_stats_frame,
            text="Export HTML report",
            command=self.export_html_report,
        )
        self.export_html_btn.pack(anchor="e", pady=(6, 0))

        self.summary_charts_frame = ttk.Frame(self.summary_frame)
        self.summary_charts_frame.pack(expand=True, fill="both", padx=10, pady=10)

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

        # Auto export mode (faster run, no interactive tables)
        self.auto_export_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            self.search_frame,
            text="Auto-export (CSV + HTML) only",
            variable=self.auto_export_var,
        ).pack(pady=(0, 16))

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
        self.notebook.add(self.genome_frame, text="Genome data")
        self._build_genome_table()
        self._build_genome_filters()

        # --- Gene tab ---
        self.gene_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.gene_frame, text="Genetic data")
        self._build_gene_table()

        # --- Summary tab ---
        self.summary_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.summary_frame, text="Summary")
        self._build_summary_tab()

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

        auto_export = self.auto_export_var.get()
        export_dir = None
        if auto_export:
            export_dir = filedialog.askdirectory(
                title="Select folder to save CSV + HTML exports"
            )
            if not export_dir:
                return

        sources = []
        if self.ensembl_var.get():
            sources.append("Ensembl")
        if self.ena_var.get():
            sources.append("ENA")
        if self.bvbrc_var.get():
            sources.append("BV-BRC")
        if self.ncbi_var.get():
            sources.append("NCBI")

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

        # Auto-export mode does not populate interactive data tabs.
        self.set_data_tabs_enabled(not auto_export)

        # --- Use proper popup ---
        show_bar = not auto_export and "NCBI" in sources

        self.progress_text_only.config(text="Starting...")
        self.progress_text_only.pack()
        self.update_idletasks()

        if "NCBI" in sources:
            threading.Thread(
                target=self.fetch_counts_thread,
                args=(taxid, email, sources, export_dir, show_bar),
                daemon=True,
            ).start()
        elif "BV-BRC" in sources:
            threading.Thread(
                target=self.fetch_bvbrc_count_thread,
                args=(taxid, sources, export_dir, show_bar),
                daemon=True,
            ).start()
        else:
            threading.Thread(
                target=self.full_search_thread,
                args=(taxid, None, None, None, sources, show_bar, export_dir),
                daemon=True,
            ).start()

    # ---------- NCBI Max Count Popup ----------
    def fetch_counts_thread(
        self,
        taxid,
        email,
        sources,
        export_dir=None,
        show_bar=True,
    ):
        try:
            genome_count, gene_count = api_client.get_total_ncbi_genome_count(
                taxid, email
            )
            bvbrc_count = None
            if "BV-BRC" in sources:
                bvbrc_count = api_client.get_total_bvbrc_genome_count(taxid)

            self.max_ncbi_genome_count = genome_count
            self.max_ncbi_gene_count = gene_count
            self.max_bvbrc_genome_count = bvbrc_count
            self.after(
                0,
                lambda: self.show_max_popup(
                    genome_count,
                    gene_count,
                    bvbrc_count,
                    taxid,
                    sources,
                    export_dir,
                    show_bar,
                ),
            )
        except Exception as e:
            self.after(
                0,
                lambda: messagebox.showerror(
                    "Count Failed", f"Could not fetch counts:\n{e}"
                ),
            )

    def fetch_bvbrc_count_thread(self, taxid, sources, export_dir=None, show_bar=True):
        try:
            bvbrc_count = api_client.get_total_bvbrc_genome_count(taxid)
            self.max_bvbrc_genome_count = bvbrc_count
            self.after(
                0,
                lambda: self.show_bvbrc_max_popup(
                    bvbrc_count, taxid, sources, export_dir, show_bar
                ),
            )
        except Exception as e:
            self.after(
                0,
                lambda: messagebox.showerror(
                    "Count Failed", f"Could not fetch BV-BRC count:\n{e}"
                ),
            )

    def show_max_popup(
        self,
        genome_count,
        gene_count,
        bvbrc_count,
        taxid,
        sources,
        export_dir=None,
        show_bar=True,
    ):
        popup = tk.Toplevel(self)
        popup.title("NCBI Max Records")
        popup.grab_set()

        ttk.Label(popup, text="Max NCBI genomic records:").grid(
            row=0, column=0, padx=6, pady=6, sticky="e"
        )
        max_genome_entry = ttk.Entry(popup, width=12)
        max_genome_entry.grid(row=0, column=1, sticky="w")
        max_genome_entry.insert(0, str(genome_count))
        ttk.Label(popup, text=f"Max: {genome_count}").grid(row=0, column=2, sticky="w")

        ttk.Label(popup, text="Max NCBI genetic records:").grid(
            row=1, column=0, padx=6, pady=6, sticky="e"
        )
        max_gene_entry = ttk.Entry(popup, width=12)
        max_gene_entry.grid(row=1, column=1, sticky="w")
        max_gene_entry.insert(0, str(gene_count))
        ttk.Label(popup, text=f"Max: {gene_count}").grid(row=1, column=2, sticky="w")

        max_bvbrc_entry = None
        button_row = 2
        if "BV-BRC" in sources:
            ttk.Label(popup, text="Max BV-BRC genomic records:").grid(
                row=2, column=0, padx=6, pady=6, sticky="e"
            )
            max_bvbrc_entry = ttk.Entry(popup, width=12)
            max_bvbrc_entry.grid(row=2, column=1, sticky="w")
            initial_bvbrc = bvbrc_count if bvbrc_count is not None else 0
            max_bvbrc_entry.insert(0, str(initial_bvbrc))
            max_text = f"Max: {bvbrc_count}" if bvbrc_count is not None else "0 = all"
            ttk.Label(popup, text=max_text).grid(row=2, column=2, sticky="w")
            button_row = 3

        def start_search():
            try:
                max_genomes_val = int(max_genome_entry.get())
                max_genes_val = int(max_gene_entry.get())
                max_bvbrc_val = (
                    int(max_bvbrc_entry.get()) if max_bvbrc_entry is not None else None
                )
            except ValueError:
                messagebox.showerror("Invalid Input", "Max records must be numeric")
                return
            popup.destroy()
            self.on_start_search(
                taxid,
                sources,
                max_genomes_val,
                max_genes_val,
                max_bvbrc_val,
                export_dir,
                show_bar,
            )

        ttk.Button(popup, text="Start Search", command=start_search).grid(
            row=button_row, column=0, columnspan=3, pady=10
        )

    def show_bvbrc_max_popup(
        self,
        bvbrc_count,
        taxid,
        sources,
        export_dir=None,
        show_bar=True,
    ):
        popup = tk.Toplevel(self)
        popup.title("BV-BRC Max Records")
        popup.grab_set()

        ttk.Label(popup, text="Max BV-BRC genomic records:").grid(
            row=0, column=0, padx=6, pady=6, sticky="e"
        )
        max_bvbrc_entry = ttk.Entry(popup, width=12)
        max_bvbrc_entry.grid(row=0, column=1, sticky="w")
        max_bvbrc_entry.insert(0, str(bvbrc_count if bvbrc_count is not None else 0))
        max_text = f"Max: {bvbrc_count}" if bvbrc_count is not None else "0 = all"
        ttk.Label(popup, text=max_text).grid(row=0, column=2, sticky="w")

        def start_search():
            try:
                max_bvbrc_val = int(max_bvbrc_entry.get())
            except ValueError:
                messagebox.showerror("Invalid Input", "Max records must be numeric")
                return

            popup.destroy()
            self.on_start_search(
                taxid,
                sources,
                None,
                None,
                max_bvbrc_val,
                export_dir,
                show_bar,
            )

        ttk.Button(popup, text="Start Search", command=start_search).grid(
            row=1, column=0, columnspan=3, pady=10
        )

    def on_start_search(
        self,
        taxid,
        sources,
        max_genomes=None,
        max_genes=None,
        max_bvbrc_genomes=None,
        export_dir=None,
        show_bar=True,
    ):
        show_bar = show_bar and "NCBI" in sources
        self.set_search_running(True)
        if show_bar:
            self.progress_bar["value"] = 0
            self.progress_percent.config(text="0%")
            self.progress_label.config(text="Starting search...")
        else:
            self.progress_text_only.config(text="Starting search...")
            self.progress_text_only.pack()
        self.update_idletasks()

        threading.Thread(
            target=self.full_search_thread,
            args=(
                taxid,
                max_genomes,
                max_genes,
                max_bvbrc_genomes,
                sources,
                show_bar,
                export_dir,
            ),
            daemon=True,
        ).start()

    def set_search_running(self, running: bool):
        """Enable/disable UI controls during a search."""
        self._search_running = running
        state = "disabled" if running else "normal"
        self.search_btn.config(state=state)
        self.taxid_entry.config(state=state)

    def set_data_tabs_enabled(self, enabled: bool):
        state = "normal" if enabled else "disabled"
        for tab_frame in (
            self.taxonomy_frame,
            self.genome_frame,
            self.gene_frame,
            self.summary_frame,
        ):
            self.notebook.tab(tab_frame, state=state)

        if not enabled:
            selected = self.notebook.select()
            disabled_tabs = {
                str(self.taxonomy_frame),
                str(self.genome_frame),
                str(self.gene_frame),
                str(self.summary_frame),
            }
            if selected in disabled_tabs:
                self.notebook.select(self.search_frame)

    def full_search_thread(
        self,
        taxid,
        max_genomes,
        max_genes,
        max_bvbrc_genomes,
        sources,
        show_bar=True,
        export_dir=None,
    ):
        # Reset data at start of search
        self.genome_data = []
        self.taxonomy_data = None

        def update_progress_label(text):
            self.after(0, lambda: self.progress_text_only.config(text=text))

        def update_progress_bar(percent):
            if show_bar:
                self.after(0, lambda: self.progress_bar.config(value=percent))
                self.after(0, lambda: self.progress_percent.config(text=f"{percent}%"))

        # --- Taxonomy ---
        update_progress_label("Fetching taxonomy data...")
        try:
            tax_data = api_client.fetch_ncbi_taxonomy(taxid)
            self.taxonomy_data = tax_data
            self.build_taxonomy_tree(tax_data)
        except Exception:
            logging.exception("Failed to fetch taxonomy")

        # --- Genomes ---
        update_progress_label("Fetching genomic data from multiple sources...")

        if "NCBI" in sources:
            self.after(0, lambda: self.progress_text_only.pack_forget())
            self.after(0, lambda: self.progress_frame.pack())
            self.after(0, lambda: self.progress_label.config(text="Fetching genomic data from multiple sources..."))

        def genome_progress(current, total):
            pct = min(100, int(current / total * 100)) if total else 0
            update_progress_bar(pct)

        try:
            summary = api_client.get_genome_summary(
                taxid,
                max_genomes or 0,
                max_bvbrc_records=max_bvbrc_genomes or 0,
                sources=sources,
                progress_callback=genome_progress,
            )
            genomes = summary.get("genomes", [])
            if not export_dir:
                self.after(
                    0, lambda g=genomes: self.append_genomes(g)
                )
            else:
                self.genome_data.extend(genomes)
        except Exception as e:
            logging.exception("Failed to fetch genomic data")
            self.after(
                0,
                lambda e=e: messagebox.showerror(
                    "Error", f"Failed to fetch genomic data:\n{e}"
                ),
            )

        # --- Genes (NCBI only) ---
        if "NCBI" in sources:
            update_progress_label("Fetching NCBI genetic data...")
            self.after(0, lambda: self.progress_label.config(text="Fetching NCBI genetic data..."))

            def gene_progress(current, total):
                pct = min(100, int(current / total * 100)) if total else 0
                update_progress_bar(pct)

            try:
                genes = api_client.get_gene_summary(
                    taxid, max_genes or 0, progress_callback=gene_progress
                )
                # store genes even in export mode (optional in future)
                if not export_dir:
                    self.after(0, lambda: self.populate_gene_table(genes))
            except Exception as e:
                logging.exception("Failed to fetch genetic data")
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

        if export_dir:
            try:
                csv_path, html_path = reports._export_results(
                    export_dir,
                    taxid,
                    self.genome_data,
                    taxonomy_data=self.taxonomy_data,
                )
                self.after(
                    0,
                    lambda: messagebox.showinfo(
                        "Exported",
                        f"Saved CSV and HTML to:\n{csv_path}\n{html_path}",
                    ),
                )
            except Exception as e:
                logging.exception("Failed to export results")
                self.after(
                    0,
                    lambda e=e: messagebox.showerror(
                        "Error", f"Could not export results:\n{e}"
                    ),
                )
        else:
            self.after(0, self.update_summary)

        self.after(0, lambda: self.set_search_running(False))

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

    def append_genomes(self, genomes):
        """Add newly fetched genomes and refresh the table."""
        self.genome_data.extend(genomes)
        self.populate_genome_table(self.genome_data)

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

    def normalize_assembly(self, level):
        if not level:
            return "Unknown"

        level = level.strip().lower().replace("_", " ")

        mapping = {
            "complete genome": "Complete Genome",
            "chromosome": "Chromosome",
            "scaffold": "Scaffold",
            "contig": "Contig",
            "primary assembly": "Primary Assembly",
        }

    def update_summary(self):
        if not self.genome_data:
            return

        # Basic stats
        total = len(self.genome_data)

        accessions = [
            g.get("Accession") for g in self.genome_data if g.get("Accession")
        ]
        unique_accessions = len(set(accessions))
        duplicates = total - unique_accessions

        source_counts = Counter(g.get("Source") or "Unknown" for g in self.genome_data)

        stats_text = (
            f"Total genome records: {total}\n"
            f"Unique accessions: {unique_accessions}\n"
            f"Duplicate records (cross-database overlap): {duplicates}\n\n"
            f"Records per database:\n"
        )

        for source, count in source_counts.items():
            stats_text += f"  - {source}: {count}\n"

        self.summary_label.config(text=stats_text)

        # Clear previous charts
        for widget in self.summary_charts_frame.winfo_children():
            widget.destroy()

        fig = plt.Figure(figsize=(12, 8))

        # Color scheme
        colors = {
            "NCBI": "#d62728",
            "ENA": "#1f77b4",
            "Ensembl": "#2ca02c",
            "BV-BRC": "#9467bd",
            "Unknown": "gray",
        }

        sources = sorted(source_counts.keys())

        # Records per Database
        ax1 = fig.add_subplot(221)
        ax1.bar(
            source_counts.keys(),
            source_counts.values(),
            color=[colors.get(s, "gray") for s in source_counts.keys()],
        )
        ax1.set_title("Records per Database")
        ax1.tick_params(axis="x", rotation=45)

        # Assembly Level

        assembly_db_counts = defaultdict(lambda: defaultdict(int))

        for g in self.genome_data:
            assembly = reports.normalize_assembly(g.get("Assembly Level"))
            source = g.get("Source") or "Unknown"
            assembly_db_counts[assembly][source] += 1

        assemblies = sorted(assembly_db_counts.keys())

        ax2 = fig.add_subplot(222)

        bottom = np.zeros(len(assemblies))

        for source in sources:
            values = [
                assembly_db_counts[assembly].get(source, 0) for assembly in assemblies
            ]

            ax2.bar(
                assemblies,
                values,
                bottom=bottom,
                label=source,
                color=colors.get(source, "gray"),
            )

            bottom += np.array(values)

        ax2.set_title("Assembly Level Distribution (by Database)")
        ax2.tick_params(axis="x", rotation=30)
        ax2.legend()

        # Sequence Length
        seq_data = []

        for g in self.genome_data:
            try:
                length = int(g.get("Seq Length"))
                source = g.get("Source") or "Unknown"
                seq_data.append((length, source))
            except:
                continue

        if seq_data:
            ax3 = fig.add_subplot(212)

            all_lengths = [x[0] for x in seq_data]
            bins = np.histogram_bin_edges(all_lengths, bins="auto")

            for source in sources:
                source_lengths = [length for length, s in seq_data if s == source]

                ax3.hist(
                    source_lengths,
                    bins=bins,
                    stacked=True,
                    label=source,
                    color=colors.get(source, "gray"),
                )

            ax3.set_title("Sequence Length Distribution (by Database)")
            ax3.set_xlabel("Sequence Length (bp)")
            ax3.set_ylabel("Count")
            ax3.legend()

        fig.tight_layout(pad=3.0)

        canvas = FigureCanvasTkAgg(fig, master=self.summary_charts_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(expand=True, fill="both")

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
        reports._write_genome_csv(path, self.genome_data)
        messagebox.showinfo("Exported", f"Saved to {path}")

    def export_html_report(self):
        if not self.genome_data:
            messagebox.showinfo("No data", "No genome data to export.")
            return

        path = filedialog.asksaveasfilename(
            defaultextension=".html", filetypes=[("HTML", "*.html")]
        )
        if not path:
            return

        try:
            reports._write_html_report(
                path,
                self.genome_data,
                taxonomy_data=self.taxonomy_data,
                open_after=True,
            )
            messagebox.showinfo("Exported", f"Saved report to {path}")
        except Exception as e:
            logging.exception("Failed to export HTML report")
            messagebox.showerror("Error", f"Could not export HTML report:\n{e}")

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
