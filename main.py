import logging
import base64
import os
from io import BytesIO

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import webbrowser
import threading
from concurrent.futures import ThreadPoolExecutor
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
        self.max_ncbi_genome_count = None
        self.max_ncbi_gene_count = None
        self.max_bvbrc_genome_count = None
        self.genome_data = []
        self.taxonomy_data = None

        # --- runtime state ---
        self._search_running = False
        self._cancel_event = None

        # --- build UI ---
        self.create_widgets()

    def _build_summary_tab(self):
        # Create a frame for the header/button
        header_frame = ttk.Frame(self.summary_frame)
        header_frame.pack(fill="x", padx=10, pady=10)

        self.export_html_btn = ttk.Button(
            header_frame,
            text="Export HTML report",
            command=self.export_html_report,
        )
        self.export_html_btn.pack(side="right")

        # Create scrollable frame for stats and charts
        canvas = tk.Canvas(self.summary_frame, highlightthickness=0)
        scrollbar = ttk.Scrollbar(self.summary_frame, orient="vertical", command=canvas.yview)
        self.summary_stats_frame = ttk.Frame(canvas)
        
        self.summary_stats_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=self.summary_stats_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True, padx=10, pady=10)
        scrollbar.pack(side="right", fill="y")

        # Bind mousewheel to canvas
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas.bind_all("<MouseWheel>", _on_mousewheel)

        self.summary_label = ttk.Label(
            self.summary_stats_frame, text="", justify="left"
        )
        self.summary_label.pack(anchor="w")

        self.summary_charts_frame = ttk.Frame(self.summary_stats_frame)
        self.summary_charts_frame.pack(expand=True, fill="both", pady=10)

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

        # Search controls
        controls_frame = ttk.Frame(self.search_frame)
        controls_frame.pack(pady=16)
        self.search_btn = ttk.Button(
            controls_frame, text="Search", command=self.on_search_button
        )
        self.search_btn.pack(side="left", padx=(0, 8))
        self.quick_summary_search_btn = ttk.Button(
            controls_frame,
            text="Quick Summary",
            command=self.show_quick_summary,
        )
        self.quick_summary_search_btn.pack(side="left", padx=(0, 8))
        self.cancel_btn = ttk.Button(
            controls_frame,
            text="Cancel",
            command=self.cancel_search,
            state="disabled",
        )
        self.cancel_btn.pack(side="left")

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
                args=(taxid, sources, export_dir, show_bar),
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
        sources,
        export_dir=None,
        show_bar=True,
    ):
        try:
            genome_count, gene_count = api_client.get_total_ncbi_genome_count(taxid)
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
                lambda e=e: messagebox.showerror(
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
                lambda e=e: messagebox.showerror(
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
            max_text = f"Max: {bvbrc_count}" if bvbrc_count is not None else "0 = skip"
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
        max_text = f"Max: {bvbrc_count}" if bvbrc_count is not None else "0 = skip"
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
        select_summary_on_complete=False,
    ):
        ncbi_genome_enabled = "NCBI" in sources and (max_genomes is None or max_genomes > 0)
        bvbrc_genome_enabled = "BV-BRC" in sources and (
            max_bvbrc_genomes is None or max_bvbrc_genomes > 0
        )
        show_bar = show_bar and (ncbi_genome_enabled or bvbrc_genome_enabled)
        self._cancel_event = threading.Event()
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
                select_summary_on_complete,
            ),
            daemon=True,
        ).start()

    def set_search_running(self, running: bool):
        """Enable/disable UI controls during a search."""
        self._search_running = running
        if not running:
            self._cancel_event = None
        search_state = "disabled" if running else "normal"
        self.search_btn.config(state=search_state)
        self.taxid_entry.config(state=search_state)
        self.cancel_btn.config(state="normal" if running else "disabled")

    def cancel_search(self):
        if not self._search_running or not self._cancel_event:
            return
        self._cancel_event.set()
        self.progress_text_only.config(text="Cancelling... waiting for active requests to stop")
        self.progress_text_only.pack()
        if self.progress_frame.winfo_ismapped():
            self.progress_label.config(text="Cancelling... waiting for active requests to stop")

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
        select_summary_on_complete=False,
    ):
        cancel_event = self._cancel_event

        def is_cancelled():
            return bool(cancel_event and cancel_event.is_set())

        # Reset data at start of search
        self.genome_data = []
        self.taxonomy_data = None

        def update_progress_label(text):
            self.after(0, lambda: self.progress_text_only.config(text=text))

        def update_progress_bar(percent):
            if show_bar:
                capped = min(int(percent), 90)
                self.after(0, lambda: self.progress_bar.config(value=capped))
                self.after(0, lambda: self.progress_percent.config(text=f"{capped}%"))

        genome_sources = list(sources)
        if "NCBI" in genome_sources and (max_genomes is not None and max_genomes <= 0):
            genome_sources.remove("NCBI")
        if "BV-BRC" in genome_sources and (
            max_bvbrc_genomes is not None and max_bvbrc_genomes <= 0
        ):
            genome_sources.remove("BV-BRC")
        fetch_genes = "NCBI" in sources and (max_genes is not None and max_genes > 0)

        # --- Taxonomy ---
        update_progress_label("Fetching taxonomy data...")
        try:
            tax_data = api_client.fetch_ncbi_taxonomy(taxid, stop_event=cancel_event)
            self.taxonomy_data = tax_data
            if not is_cancelled():
                self.build_taxonomy_tree(tax_data)
        except Exception:
            logging.exception("Failed to fetch taxonomy")

        if is_cancelled():
            self.after(0, lambda: self.progress_text_only.config(text="Cancelled"))
            if show_bar:
                self.after(0, self.progress_frame.pack_forget)
            else:
                self.after(0, self.progress_text_only.pack_forget)
            self.after(0, lambda: self.set_search_running(False))
            return

        # --- Genomes + Genes (parallel when enabled) ---
        update_progress_label("Fetching genomic data from multiple sources...")

        if "NCBI" in genome_sources or "BV-BRC" in genome_sources:
            self.after(0, lambda: self.progress_text_only.pack_forget())
            self.after(0, lambda: self.progress_frame.pack())
            self.after(
                0,
                lambda: self.progress_label.config(
                    text="Fetching genomic and genetic data..."
                ),
            )

        tracked_sources = []
        if "NCBI" in genome_sources:
            tracked_sources.append("NCBI")
        if "BV-BRC" in genome_sources:
            tracked_sources.append("BV-BRC")
        source_progress = {source: 0 for source in tracked_sources}

        def genome_progress(source, current, total):
            if source not in source_progress:
                return
            pct = min(100, int(current / total * 100)) if total else 0
            source_progress[source] = pct
            # Progress follows the slowest tracked source.
            combined_pct = min(source_progress.values()) if source_progress else pct
            update_progress_bar(combined_pct)

        def fetch_genomes_task():
            if not genome_sources:
                return {"genomes": []}
            return api_client.get_genome_summary(
                taxid,
                max_genomes or 0,
                max_bvbrc_records=max_bvbrc_genomes or 0,
                sources=genome_sources,
                progress_callback=genome_progress,
                stop_event=cancel_event,
            )

        def fetch_genes_task():
            if not fetch_genes:
                return []
            return api_client.get_gene_summary(
                taxid,
                max_genes or 0,
                stop_event=cancel_event,
            )

        summary = None
        genes = None

        with ThreadPoolExecutor(max_workers=2) as executor:
            genome_future = executor.submit(fetch_genomes_task)
            gene_future = executor.submit(fetch_genes_task) if fetch_genes else None

            try:
                summary = genome_future.result()
            except Exception as e:
                logging.exception("Failed to fetch genomic data")
                self.after(
                    0,
                    lambda e=e: messagebox.showerror(
                        "Error", f"Failed to fetch genomic data:\n{e}"
                    ),
                )

            if gene_future is not None:
                try:
                    genes = gene_future.result()
                except Exception as e:
                    logging.exception("Failed to fetch genetic data")
                    self.after(
                        0,
                        lambda e=e: messagebox.showerror(
                            "Error", f"Failed to fetch genetic data:\n{e}"
                        ),
                    )

        if summary:
            genomes = summary.get("genomes", [])
            if not export_dir and not is_cancelled():
                self.after(0, lambda g=genomes: self.append_genomes(g))
            elif not is_cancelled():
                self.genome_data.extend(genomes)

        if genes is not None and not export_dir and not is_cancelled():
            self.after(0, lambda g=genes: self.populate_gene_table(g))

        if "NCBI" in sources and not fetch_genes and not export_dir:
            self.after(0, lambda: self.populate_gene_table([]))

        if show_bar and not is_cancelled():
            update_progress_bar(100)
            self.after(0, lambda: self.progress_bar.config(value=100))
            self.after(0, lambda: self.progress_percent.config(text="100%"))

        if is_cancelled():
            update_progress_label("Cancelled")
        else:
            update_progress_label("Done")
        if show_bar:
            self.after(2000, self.progress_frame.pack_forget)
        else:
            self.after(2000, self.progress_text_only.pack_forget)

        if export_dir and not is_cancelled():
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
        elif not is_cancelled():
            self.after(0, self.update_summary)
            if select_summary_on_complete:
                self.after(0, lambda: self.notebook.select(self.summary_frame))

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

    def show_quick_summary(self):
        """Run a fast pre-scan for current taxon and open the Summary tab."""
        taxid_str = self.taxid_entry.get().strip()
        if not taxid_str.isdigit():
            messagebox.showerror("Invalid Input", "Please enter a numeric Taxonomy ID.")
            return

        if self._search_running:
            messagebox.showinfo("Quick Summary", "A search is already running.")
            return

        taxid = int(taxid_str)

        sources = []
        if self.ensembl_var.get():
            sources.append("Ensembl")
        if self.ena_var.get():
            sources.append("ENA")
        if self.bvbrc_var.get():
            sources.append("BV-BRC")
        if self.ncbi_var.get():
            sources.append("NCBI")

        if not sources:
            messagebox.showerror("No Source Selected", "Select at least one data source.")
            return

        # Fast pre-scan limits: enough data for a meaningful summary, but quick to fetch.
        quick_max_ncbi_genomes = 5000
        quick_max_ncbi_genes = 100
        quick_max_bvbrc_genomes = 1500

        self.set_data_tabs_enabled(True)
        self.progress_text_only.config(text="Starting quick summary scan...")
        self.progress_text_only.pack()
        self.update_idletasks()

        self.on_start_search(
            taxid,
            sources,
            max_genomes=quick_max_ncbi_genomes,
            max_genes=quick_max_ncbi_genes,
            max_bvbrc_genomes=quick_max_bvbrc_genomes,
            export_dir=None,
            show_bar=True,
            select_summary_on_complete=True,
        )

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
        geography = reports.summarize_geography(self.genome_data)
        all_regions = geography["top_regions"]
        top_regions = all_regions[:10]
        other_count = sum(count for _, count in all_regions[10:])
        display_regions = list(top_regions)
        if other_count > 0:
            display_regions.append(("Other", other_count))
        region_denom = geography["region_covered"] or 1
        missing_region_pct = (geography["missing_region"] / total * 100) if total else 0.0
        region_coverage_pct = (geography["region_covered"] / total * 100) if total else 0.0

        stats_text = (
            f"Total genome records: {total}\n"
            f"Unique accessions: {unique_accessions}\n"
            f"Duplicate records (cross-database overlap): {duplicates}\n\n"
            f"Records per database:\n"
        )

        for source, count in source_counts.items():
            stats_text += f"  - {source}: {count}\n"

        stats_text += (
            "\nGeographic metadata coverage:\n"
            f"  - Without region metadata: {geography['missing_region']} / {geography['total']} ({missing_region_pct:.1f}%)\n"
            f"  - With region metadata: {geography['region_covered']} / {geography['total']} ({region_coverage_pct:.1f}%)\n"
            "Top 10 collection countries + Other (share among records that have region metadata):\n"
        )

        if display_regions:
            for region, count in display_regions:
                stats_text += f"  - {region}: {count / region_denom * 100:.1f}%\n"
        else:
            stats_text += "  - No geographic metadata available\n"

        self.summary_label.config(text=stats_text)

        # Clear previous charts
        for widget in self.summary_charts_frame.winfo_children():
            widget.destroy()

        fig = plt.Figure(figsize=(12, 8))

        # Color scheme
        colors = {
            "NCBI": "#d62728",
            "REFSEQ": "#d62728",
            "GENBANK": "#ff7f0e",
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
            ax3 = fig.add_subplot(223)

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

        # Geographic Distribution
        if display_regions:
            ax4 = fig.add_subplot(224)
            regions_list = display_regions
            region_names = [r[0] for r in regions_list]
            region_counts = [r[1] for r in regions_list]
            
            # Reverse order so largest is at top
            region_names = region_names[::-1]
            region_counts = region_counts[::-1]
            
            bars = ax4.barh(region_names, region_counts, color="#2ca02c")
            
            # Add numeric labels on bars
            for i, (bar, count) in enumerate(zip(bars, region_counts)):
                ax4.text(count, i, f" {count}", va="center")
            
            ax4.set_title("Top 10 Countries + Other")
            ax4.set_xlabel("Count")

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
