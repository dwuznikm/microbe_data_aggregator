import requests

BASE_URL_NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2"
BASE_URL_ENSEMBL = "https://rest.ensembl.org"

def fetch_json_from_api(url: str):
    headers = {
        "Accept": "application/json",
        "User-Agent": "GenomicDataExplorer/0.1"
    }
    try:
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data from API: {e}")
        return None


# --- NCBI ---

def fetch_ncbi_genomes(tax_id: int):
    url = f"{BASE_URL_NCBI}/genome/taxon/{tax_id}/dataset_report"
    return fetch_json_from_api(url)


def fetch_ncbi_taxonomy(tax_id: int):
    url = f"{BASE_URL_NCBI}/taxonomy/taxon/{tax_id}/dataset_report"
    return fetch_json_from_api(url)


def extract_ncbi_genome_metadata(report: dict):
    try:
        accession = report.get("accession")
        assembly_info = report.get("assembly_info", {})
        assembly_stats = report.get("assembly_stats", {})
        annotation_info = report.get("annotation_info", {})

        assembly_level = assembly_info.get("assembly_level", "Unknown")
        seq_len = assembly_stats.get("total_sequence_length")
        gc_content = assembly_stats.get("gc_percent")
        gene_count = annotation_info.get("stats", {}).get("gene_counts", {}).get("total")

        ref_genome = "Yes" if assembly_info.get("refseq_category") == "reference genome" else "No"
        link = f"https://www.ncbi.nlm.nih.gov/assembly/{accession}" if accession else ""

        return {
            "Accession": accession,
            "Assembly Level": assembly_level,
            "Seq Length": seq_len,
            "GC Content": gc_content,
            "# Genes": gene_count,
            "Source": "NCBI",
            "Reference": ref_genome,
            "Link": link,
        }
    except Exception:
        return None


# --- Ensembl ---

def fetch_ensembl_genomes(tax_id: int):
    url = f"{BASE_URL_ENSEMBL}/info/genomes/taxonomy/{tax_id}"
    return fetch_json_from_api(url)


def extract_ensembl_genome_metadata(genome: dict):
    try:
        accession = genome.get("assembly_accession")
        assembly_level = genome.get("assembly_level", "Unknown")
        seq_len = genome.get("base_count")
        link_name = genome.get("url_name") or genome.get("name")
        link = f"https://bacteria.ensembl.org/{link_name}/Info/Index" if link_name else ""
        ref_genome = "Yes" if genome.get("reference") else "No"

        return {
            "Accession": accession,
            "Assembly Level": assembly_level,
            "Seq Length": seq_len,
            "GC Content": "",
            "# Genes": "",
            "Source": "Ensembl",
            "Reference": ref_genome,
            "Link": link,
        }
    except Exception:
        return None


# --- Combined ---

def get_genome_summary(tax_id: int, sources=None):
    """
    Fetches NCBI + Ensembl genome data, returns unified list of genome metadata dicts.
    sources: list of "NCBI" and/or "Ensembl". If None, fetch both.
    """
    if sources is None:
        sources = ["NCBI", "Ensembl"]

    taxonomy_data = fetch_ncbi_taxonomy(tax_id)
    summary = {"tax_id": tax_id, "organism": "Unknown", "genomes": []}

    # Organism name
    try:
        summary["organism"] = taxonomy_data["reports"][0]["taxonomy"]["current_scientific_name"]["name"]
    except Exception:
        pass

    # --- NCBI genomes ---
    if "NCBI" in sources:
        ncbi_data = fetch_ncbi_genomes(tax_id)
        if ncbi_data:
            reports = ncbi_data.get("reports", [])
            for r in reports:
                metadata = extract_ncbi_genome_metadata(r)
                if metadata:
                    summary["genomes"].append(metadata)

    # --- Ensembl genomes ---
    if "Ensembl" in sources:
        ensembl_data = fetch_ensembl_genomes(tax_id)
        if ensembl_data and isinstance(ensembl_data, list):
            for g in ensembl_data:
                metadata = extract_ensembl_genome_metadata(g)
                if metadata:
                    summary["genomes"].append(metadata)

    return summary
