import requests
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from queue import Queue, Empty
import threading
import time
from Bio import Entrez
import ssl

ssl._create_default_https_context = ssl._create_unverified_context

BASE_URL_NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2"
BASE_URL_ENA = "https://www.ebi.ac.uk/ena/portal/api/search"

BASE_URL_ENSEMBL = "https://rest.ensembl.org"
NCBI_API_KEY = "ad05b3160b82232233a302e84033a1eb8307"


def fetch_json_from_api(url: str, params=None, retries=25, delay=0.1):
    """
    Fetch JSON from API with retries and constant delay.
    Supports NCBI API key for increased rate limits (10 req/sec).
    """
    headers = {"Accept": "application/json", "User-Agent": "MicrobeDataAggregator/0.1"}

    params = dict(params or {})
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    for attempt in range(retries):
        try:
            response = requests.get(url, headers=headers, params=params, timeout=60)
            response.raise_for_status()
            return response.json()

        except requests.exceptions.HTTPError as e:
            status = response.status_code if "response" in locals() else None

            # Retry for transient server or rate-limit errors
            if status in (429, 500, 502, 503, 504):
                print(f"Server error {status}, retrying ({attempt+1}/{retries})...")
                time.sleep(delay)
                continue

            print(f"Non-retryable HTTP error {status}: {e}")
            break

        except requests.exceptions.RequestException as e:
            print(f"Request failed: {e}")
            time.sleep(delay)
            continue

    print(f"Giving up after {retries} attempts for URL: {url}")
    return None


def fetch_ncbi_genomes(tax_id: int, max_records, max_workers=6):
    url = f"{BASE_URL_NCBI}/genome/taxon/{tax_id}/dataset_report"
    all_reports = []
    token_queue = Queue()
    lock = threading.Lock()
    stop_flag = threading.Event()

    print(
        f"🚀 Starting concurrent genome fetch for tax_id={tax_id} ({max_workers} threads)..."
    )

    # --- Initial request (bootstrap) ---
    first = fetch_json_from_api(url)
    if not first:
        print("❌ Failed initial fetch.")
        return {"reports": []}

    reports = first.get("reports", [])
    with lock:
        all_reports.extend(reports)
    print(f"✅ Page 1 fetched: {len(reports)} records (total: {len(all_reports)})")

    next_token = first.get("next_page_token")
    if not next_token:
        print("Only one page found.")
        return {"reports": all_reports}

    token_queue.put(next_token)

    def worker():
        nonlocal all_reports
        while not stop_flag.is_set():
            try:
                token = token_queue.get(timeout=2)
            except Empty:
                # If the queue is empty and stop_flag is set, exit gracefully
                if stop_flag.is_set():
                    break
                continue

            # Stop worker if sentinel is encountered
            if token is None:
                token_queue.task_done()
                break

            # Fetch next page
            data = fetch_json_from_api(url, params={"page_token": token})
            if not data:
                token_queue.task_done()
                continue

            new_reports = data.get("reports", [])
            next_token = data.get("next_page_token")

            with lock:
                all_reports.extend(new_reports)
                total = len(all_reports)
                print(
                    f"🧩 {threading.current_thread().name}: +{len(new_reports)} (total={total})"
                )
                if total >= max_records:
                    stop_flag.set()

            # Queue the next token if allowed
            if next_token and not stop_flag.is_set():
                token_queue.put(next_token)

            token_queue.task_done()

    # Start workers
    threads = [
        threading.Thread(target=worker, name=f"Worker-{i}", daemon=True)
        for i in range(max_workers)
    ]
    for t in threads:
        t.start()

    # --- Main thread watcher loop ---
    while True:
        time.sleep(0.5)
        with lock:
            total = len(all_reports)
        if total >= max_records:
            stop_flag.set()
        # Stop when both conditions hold: queue empty + all workers idle
        if token_queue.empty() and all(
            not t.is_alive() or stop_flag.is_set() for t in threads
        ):
            break
        # Safety timeout to avoid infinite wait
        if stop_flag.is_set() and token_queue.empty():
            break

    # --- Clean shutdown ---
    while not token_queue.empty():
        try:
            token_queue.get_nowait()
            token_queue.task_done()
        except Empty:
            break
    for _ in threads:
        token_queue.put(None)
    for t in threads:
        t.join(timeout=0.5)

    print(f"✅ Completed: {len(all_reports)} genome records fetched across threads.\n")
    return {"reports": all_reports[:max_records]}


def fetch_ncbi_taxonomy(tax_id: int):
    """
    Fetch taxonomy report for a given tax_id.
    """
    url = f"{BASE_URL_NCBI}/taxonomy/taxon/{tax_id}/dataset_report"
    return fetch_json_from_api(url)


def extract_ncbi_genome_metadata(report: dict):
    """
    Extract and normalize genome metadata from an NCBI genome report.
    """
    try:
        accession = report.get("accession")
        source_database = report.get("source_database")
        assembly_info = report.get("assembly_info", {})
        assembly_stats = report.get("assembly_stats", {})
        annotation_info = report.get("annotation_info", {})

        assembly_level = assembly_info.get("assembly_level", "Unknown")
        seq_len = assembly_stats.get("total_sequence_length")
        gc_content = assembly_stats.get("gc_percent")
        gene_count = (
            annotation_info.get("stats", {}).get("gene_counts", {}).get("total")
        )

        ref_genome = (
            "Yes"
            if assembly_info.get("refseq_category") == "reference genome"
            else "No"
        )
        link = f"https://www.ncbi.nlm.nih.gov/assembly/{accession}" if accession else ""

        return {
            "Accession": accession,
            "Assembly Level": assembly_level,
            "Seq Length": seq_len,
            "GC Content": gc_content,
            "# Genes": gene_count,
            "Source": source_database,
            "Reference": ref_genome,
            "Link": link,
        }
    except Exception:
        return None


# --- NCBI Gene ---


def get_total_ncbi_genome_count(tax_id: int, email):
    """
    Fetch all gene reports for a given tax_id from NCBI Datasets API with pagination.
    Displays progress with page numbers and total records fetched.
    """
    Entrez.email = email
    handle_genome = Entrez.esearch(db="assembly", term=f"txid{tax_id}[Organism:exp]")
    record_genome = Entrez.read(handle_genome)
    handle_gene = Entrez.esearch(db="gene", term=f"txid{tax_id}[Organism:exp]")
    record_gene = Entrez.read(handle_gene)
    return record_genome["Count"], record_gene["Count"]


def fetch_ncbi_genes(tax_id: int, max_records):
    """
    Fetch all gene reports for a given tax_id from NCBI Datasets API with pagination.
    Displays progress with page numbers and total records fetched.
    """
    url = f"{BASE_URL_NCBI}/gene/taxon/{tax_id}"
    all_reports = []
    next_page_token = None
    page_count = 0

    while True:
        params = {}
        if next_page_token:
            params["page_token"] = next_page_token

        page_count += 1
        print(f"Fetching gene page {page_count}...")

        data = fetch_json_from_api(url, params)
        if not data:
            print("No data returned or request failed.")
            break

        reports = data.get("reports", [])
        all_reports.extend(reports)

        print(
            f"→ Page {page_count} fetched: {len(reports)} records (total: {len(all_reports)})."
        )
        if len(all_reports) >= max_records:
            break
        next_page_token = data.get("next_page_token")
        if not next_page_token:
            print(
                f"Completed fetching {page_count} pages ({len(all_reports)} total gene reports)."
            )
            break
    return {"reports": all_reports}


def extract_ncbi_gene_metadata(report: dict):
    """
    Extract and normalize gene metadata from an NCBI gene report.
    """
    try:
        gene = report.get("gene", {})
        gene_id = (
            gene.get("gene_id")
            or report.get("gene_id")
            or report.get("accession")
            or report.get("uid")
        )
        symbol = gene.get("symbol", "")
        gene_type = gene.get("type", "")
        description = gene.get("description", "")
        locus_tag = gene.get("locus_tag", "")

        assemblies = []
        annotations = gene.get("annotations", [])
        for ann in annotations:
            acc = ann.get("assembly_accession")
            date = ann.get("annotation_release_date", "")
            if acc:
                assemblies.append(f"{acc} ({date})")

        assemblies_str = "; ".join(assemblies) if assemblies else "-"

        link = f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}" if gene_id else ""

        return {
            "Gene ID": gene_id,
            "Type": gene_type,
            "Symbol": symbol,
            "Locus": locus_tag,
            "Description": description,
            "Assemblies": assemblies_str,
            "Link": link,
        }
    except Exception as e:
        print(f"Error extracting gene metadata: {e}")
        return None


def get_gene_summary(tax_id: int, max_records):
    """Fetches and returns all NCBI gene metadata for a given tax_id."""
    data = fetch_ncbi_genes(tax_id, max_records)
    if not data:
        return []

    genes = []
    for r in data.get("reports", []):
        meta = extract_ncbi_gene_metadata(r)
        if meta:
            genes.append(meta)
    return genes


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
        link = (
            f"https://bacteria.ensembl.org/{link_name}/Info/Index" if link_name else ""
        )
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


# --- ENA ---


def fetch_ena_genomes(tax_id: int):
    url = f"{BASE_URL_ENA}?result=assembly&query=tax_eq({tax_id})&fields=accession,assembly_level,base_count&format=json"
    return fetch_json_from_api(url)


def extract_ena_genome_metadata(genome: dict):
    try:
        accession = genome.get("accession")
        assembly_level = genome.get("assembly_level", "Unknown")
        seq_len = genome.get("base_count")

        link = (
            f"https://www.ebi.ac.uk/ena/browser/view/{accession}" if accession else ""
        )

        return {
            "Accession": accession,
            "Assembly Level": assembly_level,
            "Seq Length": seq_len,
            "GC Content": "",
            "# Genes": "",
            "Source": "ENA",
            "Reference": "No",
            "Link": link,
        }
    except Exception:
        return None


# --- Combined ---


def get_genome_summary(tax_id: int, max_records, sources=None):
    """
    Fetches NCBI + Ensembl genome data, returns unified list of genome metadata dicts.
    sources: list of "NCBI" and/or "Ensembl". If None, fetch both.
    """
    if sources is None:
        sources = ["NCBI", "Ensembl", "ENA"]

    taxonomy_data = fetch_ncbi_taxonomy(tax_id)
    summary = {"tax_id": tax_id, "organism": "Unknown", "genomes": []}

    try:
        summary["organism"] = taxonomy_data["reports"][0]["taxonomy"][
            "current_scientific_name"
        ]["name"]
    except Exception:
        pass

    # --- NCBI genomes ---
    if "NCBI" in sources:
        ncbi_data = fetch_ncbi_genomes(tax_id, max_records)
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

    # --- ENA genomes ---
    if "ENA" in sources:
        ena_data = fetch_ena_genomes(tax_id)

        if ena_data and isinstance(ena_data, list):
            for g in ena_data:
                metadata = extract_ena_genome_metadata(g)
                if metadata:
                    summary["genomes"].append(metadata)

    return summary
