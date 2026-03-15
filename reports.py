import os
import base64
from io import BytesIO
import csv
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import numpy as np
import webbrowser


def normalize_assembly(level):
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

    return mapping.get(level, level.title())


def _write_genome_csv(path, genomes):
    cols = (
        "Accession",
        "Assembly Level",
        "Seq Length",
        "GC Content",
        "# Genes",
        "Source",
        "Reference",
        "Link",
    )
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for g in genomes:
            w.writerow(
                [
                    g.get("Accession", ""),
                    g.get("Assembly Level", ""),
                    g.get("Seq Length", ""),
                    g.get("GC Content", ""),
                    g.get("# Genes", ""),
                    g.get("Source", ""),
                    g.get("Reference", ""),
                    g.get("Link", ""),
                ]
            )


def _write_html_report(path, genomes, taxid=None, open_after=True):
    # Filter to reference genomes only for the report
    ref_genomes = [g for g in genomes if g.get("Reference") == "Yes"]

    def _figure_to_base64(fig):
        buf = BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight")
        buf.seek(0)
        data = base64.b64encode(buf.read()).decode("ascii")
        plt.close(fig)
        return data

    # build charts
    source_counts = Counter(g.get("Source") or "Unknown" for g in ref_genomes)
    colors = {
        "NCBI": "#d62728",
        "ENA": "#1f77b4",
        "Ensembl": "#2ca02c",
        "BV-BRC": "#9467bd",
        "Unknown": "gray",
    }
    sources = sorted(source_counts.keys())

    fig = plt.Figure(figsize=(12, 8))
    ax1 = fig.add_subplot(221)
    ax1.bar(
        source_counts.keys(),
        source_counts.values(),
        color=[colors.get(s, "gray") for s in source_counts.keys()],
    )
    ax1.set_title("Records per Database")
    ax1.tick_params(axis="x", rotation=45)

    assembly_db_counts = defaultdict(lambda: defaultdict(int))
    for g in ref_genomes:
        assembly = normalize_assembly(g.get("Assembly Level"))
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

    seq_data = []
    for g in ref_genomes:
        try:
            length = int(g.get("Seq Length"))
            source = g.get("Source") or "Unknown"
            seq_data.append((length, source))
        except Exception:
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

    img_data = _figure_to_base64(fig)

    total = len(ref_genomes)
    accessions = [g.get("Accession") for g in ref_genomes if g.get("Accession")]
    unique_accessions = len(set(accessions))
    duplicates = total - unique_accessions
    source_counts = Counter(g.get("Source") or "Unknown" for g in ref_genomes)

    rows = []
    for g in ref_genomes:
        rows.append(
            """
                <tr>
                  <td>{acc}</td>
                  <td>{level}</td>
                  <td>{length}</td>
                  <td>{gc}</td>
                  <td>{genes}</td>
                  <td>{source}</td>
                  <td>{ref}</td>
                  <td><a href=\"{link}\" target=\"_blank\">{link}</a></td>
                </tr>
                """.format(
                acc=g.get("Accession", ""),
                level=g.get("Assembly Level", ""),
                length=g.get("Seq Length", ""),
                gc=g.get("GC Content", ""),
                genes=g.get("# Genes", ""),
                source=g.get("Source", ""),
                ref=g.get("Reference", ""),
                link=g.get("Link", ""),
            )
        )

    html = f"""
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>Genome Data Report</title>
  <style>
    body {{ font-family: system-ui, sans-serif; margin: 18px; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 16px; }}
    th, td {{ border: 1px solid #ccc; padding: 6px 8px; text-align: left; }}
    th {{ background: #f2f2f2; }}
    img {{ max-width: 100%; height: auto; }}
    .section {{ margin-top: 24px; }}
  </style>
</head>
<body>
  <h1>Genome Data Report</h1>
  <div class="section">
    <h2>Summary</h2>
    <p>Total reference records: {total}</p>
    <p>Unique accessions: {unique_accessions}</p>
    <p>Duplicate records: {duplicates}</p>
    <p>Records per database:</p>
    <ul>
      {''.join(f'<li>{src}: {cnt}</li>' for src, cnt in source_counts.items())}
    </ul>
  </div>
  <div class="section">
    <h2>Charts</h2>
    <img src="data:image/png;base64,{img_data}" alt="Summary charts" />
  </div>
  <div class="section">
    <h2>Genome Table (Reference Genomes Only)</h2>
    <table>
      <thead>
        <tr>
          <th>Accession</th>
          <th>Assembly Level</th>
          <th>Seq Length</th>
          <th>GC Content</th>
          <th># Genes</th>
          <th>Source</th>
          <th>Reference</th>
          <th>Link</th>
        </tr>
      </thead>
      <tbody>
        {''.join(rows)}
      </tbody>
    </table>
  </div>
</body>
</html>
"""

    with open(path, "w", encoding="utf-8") as f:
        f.write(html)

    if open_after:
        webbrowser.open(f"file://{path}")


def _export_results(export_dir, taxid, genomes):
    from datetime import datetime
    timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    base_name = f"genome_report_{taxid or 'unknown'}_{timestamp}"
    csv_path = os.path.join(export_dir, f"{base_name}.csv")
    html_path = os.path.join(export_dir, f"{base_name}.html")

    # Filter to reference genomes for export
    ref_genomes = [g for g in genomes if g.get("Reference") == "Yes"]
    _write_genome_csv(csv_path, ref_genomes)
    _write_html_report(html_path, ref_genomes, taxid=taxid, open_after=False)

    return csv_path, html_path