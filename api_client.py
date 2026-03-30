import requests
import time
import os
import logging
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from queue import Queue, Empty
import threading
import time
import ssl

ssl._create_default_https_context = ssl._create_unverified_context

BASE_URL_NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2"
BASE_URL_ENA = "https://www.ebi.ac.uk/ena/portal/api/search"
BASE_URL_BVBRC = "https://www.bv-brc.org/api"

BASE_URL_ENSEMBL = "https://rest.ensembl.org"
NCBI_API_KEY = "ad05b3160b82232233a302e84033a1eb8307"
NCBI_PAGE_SIZE = 1000

LOG_DIR = os.path.join(os.path.dirname(__file__), "logs")
os.makedirs(LOG_DIR, exist_ok=True)
FETCH_TIMING_LOG_PATH = os.path.join(LOG_DIR, "fetch_timing.log")

fetch_timing_logger = logging.getLogger("microbe_data_aggregator.fetch_timing")
if not fetch_timing_logger.handlers:
    handler = logging.FileHandler(FETCH_TIMING_LOG_PATH, encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s"))
    fetch_timing_logger.addHandler(handler)
fetch_timing_logger.setLevel(logging.INFO)
fetch_timing_logger.propagate = False


def _log_fetch_timing(operation, tax_id, start_time, status="ok", records=None, details=None):
    elapsed = time.perf_counter() - start_time
    parts = [
        f"operation={operation}",
        f"tax_id={tax_id}",
        f"status={status}",
        f"elapsed_s={elapsed:.3f}",
    ]
    if records is not None:
        parts.append(f"records={records}")
    if details:
        parts.append(f"details={details}")
    fetch_timing_logger.info(" | ".join(parts))


def fetch_json_from_api(url: str, database=None, params=None, retries=None, delay=0.1):
    """
    Fetch JSON from API with retries and constant delay.
    Supports NCBI API key for increased rate limits (10 req/sec).
    """
    if retries is None:
        retries = 25 if database == "NCBI" else 2
    headers = {"Accept": "application/json", "User-Agent": "MicrobeDataAggregator/0.1"}

    params = dict(params or {})
    if NCBI_API_KEY and database == "NCBI":
        params["api_key"] = NCBI_API_KEY

    for attempt in range(retries):
        try:
            response = requests.get(url, headers=headers, params=params, timeout=30)
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


def fetch_json_and_headers_from_api(
    url: str, database=None, params=None, retries=None, delay=0.1
):
    """Fetch JSON payload with response headers."""
    if retries is None:
        retries = 25 if database == "NCBI" else 2
    headers = {"Accept": "application/json", "User-Agent": "MicrobeDataAggregator/0.1"}

    params = dict(params or {})
    if NCBI_API_KEY and database == "NCBI":
        params["api_key"] = NCBI_API_KEY

    for attempt in range(retries):
        try:
            response = requests.get(url, headers=headers, params=params, timeout=30)
            response.raise_for_status()
            return response.json(), dict(response.headers)

        except requests.exceptions.HTTPError as e:
            status = response.status_code if "response" in locals() else None

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
    return None, {}


def _parse_content_range_total(content_range):
    """Parse BV-BRC Content-Range like 'items 0-24/1234' and return total as int."""
    if not content_range:
        return None
    try:
        total_part = str(content_range).split("/")[-1].strip()
        if total_part == "*":
            return None
        return int(total_part)
    except Exception:
        return None


def get_total_bvbrc_genome_count(tax_id: int):
    """Get total BV-BRC genome records for a taxon using Content-Range headers."""
    url = (
        f"{BASE_URL_BVBRC}/genome/?eq(taxon_id,{tax_id})"
        f"&limit(1,0)&http_accept=application/json"
    )
    data, headers = fetch_json_and_headers_from_api(url)
    if data is None:
        raise RuntimeError("Could not fetch BV-BRC genome count")

    total = _parse_content_range_total(
        headers.get("Content-Range") or headers.get("content-range")
    )

    if total is not None:
        return total

    if isinstance(data, list):
        return len(data)

    return 0


def fetch_ncbi_genomes(tax_id: int, max_records, max_workers=6, progress_callback=None):
    url = f"{BASE_URL_NCBI}/genome/taxon/{tax_id}/dataset_report"
    start_time = time.perf_counter()
    limit = int(max_records or 0)
    all_reports = []
    token_queue = Queue()
    lock = threading.Lock()
    stop_flag = threading.Event()

    print(
        f"Starting concurrent genome fetch for tax_id={tax_id} ({max_workers} threads)..."
    )

    # --- Initial request (bootstrap) ---
    first = fetch_json_from_api(url, "NCBI", params={"page_size": NCBI_PAGE_SIZE})
    if not first:
        print("Failed initial fetch.")
        _log_fetch_timing(
            "genome_fetch_ncbi",
            tax_id,
            start_time,
            status="no_data",
            records=0,
            details="initial request failed",
        )
        return {"reports": []}

    reports = first.get("reports", [])
    with lock:
        all_reports.extend(reports)
        if limit > 0 and len(all_reports) > limit:
            del all_reports[limit:]

    if progress_callback:
        target = limit if limit > 0 else len(all_reports)
        progress_callback(len(all_reports), target)

    if limit > 0 and len(all_reports) >= limit:
        print(f"Reached requested max genome records ({limit}) on first page.")
        _log_fetch_timing(
            "genome_fetch_ncbi",
            tax_id,
            start_time,
            status="max_reached",
            records=len(all_reports[:limit]),
            details="reached limit on first page",
        )
        return {"reports": all_reports[:limit]}

    next_token = first.get("next_page_token")
    if not next_token:
        print("Only one page found.")
        final_reports = all_reports[:limit] if limit > 0 else all_reports
        _log_fetch_timing(
            "genome_fetch_ncbi",
            tax_id,
            start_time,
            status="single_page",
            records=len(final_reports),
        )
        return {"reports": final_reports}

    token_queue.put(next_token)

    def worker():
        nonlocal all_reports
        while not stop_flag.is_set():
            try:
                token = token_queue.get(timeout=2)
            except Empty:
                if stop_flag.is_set():
                    break
                continue

            if token is None:
                token_queue.task_done()
                break

            data = fetch_json_from_api(
                url,
                "NCBI",
                params={"page_size": NCBI_PAGE_SIZE, "page_token": token},
            )
            if not data:
                token_queue.task_done()
                continue

            new_reports = data.get("reports", [])
            next_token = data.get("next_page_token")

            with lock:
                all_reports.extend(new_reports)
                if limit > 0 and len(all_reports) > limit:
                    del all_reports[limit:]
                total = len(all_reports)
                if progress_callback:
                    target = limit if limit > 0 else total
                    progress_callback(total, target)
                if limit > 0 and total >= limit:
                    stop_flag.set()

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
        if limit > 0 and total >= limit:
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

    print(f"Completed: {len(all_reports)} genome records fetched across threads.\n")
    final_reports = all_reports[:limit] if limit > 0 else all_reports
    status = "max_reached" if limit > 0 and len(final_reports) >= limit else "completed"
    _log_fetch_timing(
        "genome_fetch_ncbi",
        tax_id,
        start_time,
        status=status,
        records=len(final_reports),
    )
    return {"reports": final_reports}


def fetch_ncbi_taxonomy(tax_id: int):
    """
    Fetch taxonomy report for a given tax_id.
    """
    start_time = time.perf_counter()
    url = f"{BASE_URL_NCBI}/taxonomy/taxon/{tax_id}/dataset_report"
    try:
        data = fetch_json_from_api(url, "NCBI")
        records = len(data.get("reports", [])) if isinstance(data, dict) else 0
        status = "ok" if data else "no_data"
        _log_fetch_timing(
            "taxonomy_fetch_ncbi",
            tax_id,
            start_time,
            status=status,
            records=records,
        )
        return data
    except Exception as e:
        _log_fetch_timing(
            "taxonomy_fetch_ncbi",
            tax_id,
            start_time,
            status="error",
            details=str(e),
        )
        raise


def extract_ncbi_genome_metadata(report: dict):
    """
    Extract and normalize genome metadata from an NCBI genome report.
    """
    try:
        accession = report.get("accession")
        source_database = report.get("source_database")
        source_database = (
            "REFSEQ"
            if report.get("source_database") == "SOURCE_DATABASE_REFSEQ"
            else "GENBANK"
        )
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


def get_total_ncbi_genome_count(tax_id: int):
    """Fetch genome and gene total_count values from NCBI Datasets endpoints."""
    genome_url = f"{BASE_URL_NCBI}/genome/taxon/{tax_id}/dataset_report"
    gene_url = f"{BASE_URL_NCBI}/gene/taxon/{tax_id}"

    genome_data = fetch_json_from_api(genome_url, "NCBI")
    gene_data = fetch_json_from_api(gene_url, "NCBI")

    if not genome_data:
        raise RuntimeError("Could not fetch NCBI genome total_count")
    if not gene_data:
        raise RuntimeError("Could not fetch NCBI gene total_count")

    return int(genome_data.get("total_count", 0)), int(gene_data.get("total_count", 0))


def fetch_ncbi_genes(tax_id: int, max_records, progress_callback=None):
    """
    Fetch all gene reports for a given tax_id from NCBI Datasets API with pagination.
    Displays progress with page numbers and total records fetched.
    """
    start_time = time.perf_counter()
    url = f"{BASE_URL_NCBI}/gene/taxon/{tax_id}"
    all_reports = []
    next_page_token = None
    page_count = 0

    limit = int(max_records or 0)
    end_status = "completed"
    end_details = None

    while True:
        params = {"page_size": NCBI_PAGE_SIZE}
        if next_page_token:
            params["page_token"] = next_page_token

        page_count += 1
        print(f"Fetching gene page {page_count}...")

        data = fetch_json_from_api(url, "NCBI", params)
        if not data:
            print("No data returned or request failed.")
            end_status = "no_data"
            end_details = "request failed or empty response"
            break

        reports = data.get("reports", [])
        all_reports.extend(reports)

        if limit > 0 and len(all_reports) >= limit:
            all_reports = all_reports[:limit]

        if progress_callback:
            target = limit if limit > 0 else len(all_reports)
            progress_callback(len(all_reports), target)

        print(
            f"→ Page {page_count} fetched: {len(reports)} records (total: {len(all_reports)})."
        )

        if limit > 0 and len(all_reports) >= limit:
            print(
                f"Reached requested max gene records ({limit}). "
                f"Completed after {page_count} pages."
            )
            end_status = "max_reached"
            end_details = f"completed after {page_count} pages"
            break

        next_page_token = data.get("next_page_token")
        if not next_page_token:
            print(
                f"Completed fetching {page_count} pages ({len(all_reports)} total gene reports)."
            )
            end_status = "completed"
            end_details = f"completed after {page_count} pages"
            break

    if page_count == 0:
        print("Completed fetching 0 pages (0 total gene reports).")

    _log_fetch_timing(
        "gene_fetch_ncbi",
        tax_id,
        start_time,
        status=end_status,
        records=len(all_reports),
        details=end_details,
    )

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


def get_gene_summary(tax_id: int, max_records, progress_callback=None):
    """Fetches and returns all NCBI gene metadata for a given tax_id."""
    data = fetch_ncbi_genes(tax_id, max_records, progress_callback=progress_callback)
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
    start_time = time.perf_counter()
    url = f"{BASE_URL_ENSEMBL}/info/genomes/taxonomy/{tax_id}"
    try:
        data = fetch_json_from_api(url)
        records = len(data) if isinstance(data, list) else 0
        status = "ok" if data is not None else "no_data"
        _log_fetch_timing(
            "genome_fetch_ensembl",
            tax_id,
            start_time,
            status=status,
            records=records,
        )
        return data
    except Exception as e:
        _log_fetch_timing(
            "genome_fetch_ensembl",
            tax_id,
            start_time,
            status="error",
            details=str(e),
        )
        raise


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
    start_time = time.perf_counter()
    url = f"{BASE_URL_ENA}?result=assembly&query=tax_eq({tax_id})&fields=accession,assembly_level,base_count&format=json"
    try:
        data = fetch_json_from_api(url)
        records = len(data) if isinstance(data, list) else 0
        status = "ok" if data is not None else "no_data"
        _log_fetch_timing(
            "genome_fetch_ena",
            tax_id,
            start_time,
            status=status,
            records=records,
        )
        return data
    except Exception as e:
        _log_fetch_timing(
            "genome_fetch_ena",
            tax_id,
            start_time,
            status="error",
            details=str(e),
        )
        raise


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


# --- BV-BRC ---


def fetch_bvbrc_genomes(tax_id: int, max_records=0, page_size=250):
    start_time = time.perf_counter()
    limit = int(max_records or 0)
    start = 0
    pages = 0
    all_records = []
    end_status = "completed"
    end_details = None
    total_available = None
    planned_records = None
    planned_pages = None

    try:
        first_url = (
            f"{BASE_URL_BVBRC}/genome/?eq(taxon_id,{tax_id})"
            f"&limit({page_size},0)&http_accept=application/json"
        )
        first_page, first_headers = fetch_json_and_headers_from_api(first_url)

        if first_page is None:
            end_status = "no_data"
            _log_fetch_timing(
                "genome_fetch_bvbrc",
                tax_id,
                start_time,
                status=end_status,
                records=0,
                details="initial request failed",
            )
            return []

        if not isinstance(first_page, list):
            end_status = "unexpected_response"
            end_details = f"type={type(first_page).__name__}"
            _log_fetch_timing(
                "genome_fetch_bvbrc",
                tax_id,
                start_time,
                status=end_status,
                records=0,
                details=end_details,
            )
            return []

        total_available = _parse_content_range_total(
            first_headers.get("Content-Range") or first_headers.get("content-range")
        )

        if total_available is not None:
            planned_records = min(limit, total_available) if limit > 0 else total_available
            planned_pages = max(1, math.ceil(planned_records / page_size)) if planned_records > 0 else 0
            print(
                f"BV-BRC total available: {total_available}. "
                f"Planned fetch: {planned_records} records across ~{planned_pages} pages."
            )

        pages = 1
        all_records.extend(first_page)

        if limit > 0 and len(all_records) >= limit:
            all_records = all_records[:limit]
            end_status = "max_reached"
            end_details = "reached limit on first page"
        else:
            start = page_size

        while end_status == "completed":
            url = (
                f"{BASE_URL_BVBRC}/genome/?eq(taxon_id,{tax_id})"
                f"&limit({page_size},{start})&http_accept=application/json"
            )
            page_data = fetch_json_from_api(url)

            if not page_data:
                break

            if not isinstance(page_data, list):
                end_status = "unexpected_response"
                end_details = f"type={type(page_data).__name__}"
                break

            pages += 1
            all_records.extend(page_data)

            if limit > 0 and len(all_records) >= limit:
                all_records = all_records[:limit]
                end_status = "max_reached"
                end_details = f"completed after {pages} pages"
                break

            if len(page_data) < page_size:
                break

            start += page_size

            if planned_pages is not None and pages >= planned_pages:
                break

        detail_parts = []
        if total_available is not None:
            detail_parts.append(f"total_available={total_available}")
        if planned_pages is not None:
            detail_parts.append(f"planned_pages={planned_pages}")
        detail_parts.append(f"fetched_pages={pages}")
        if end_details:
            detail_parts.append(end_details)
        end_details = "; ".join(detail_parts)

        _log_fetch_timing(
            "genome_fetch_bvbrc",
            tax_id,
            start_time,
            status=end_status,
            records=len(all_records),
            details=end_details,
        )
        return all_records
    except Exception as e:
        _log_fetch_timing(
            "genome_fetch_bvbrc",
            tax_id,
            start_time,
            status="error",
            details=str(e),
        )
        raise


def extract_bvbrc_genome_metadata(genome: dict):
    try:
        accession = genome.get("assembly_accession")

        cds = genome.get("cds", 0)
        rrna = genome.get("rrna", 0)
        trna = genome.get("trna", 0)
        gene_count = cds + rrna + trna

        return {
            "Accession": accession,
            "Assembly Level": "",
            "Seq Length": genome.get("genome_length"),
            "GC Content": genome.get("gc_content"),
            "# Genes": gene_count,
            "Source": "BV-BRC",
            "Reference": "No",
            "Link": (
                f"https://www.bv-brc.org/view/Genome/{accession}" if accession else ""
            ),
        }
    except Exception:
        return None


# --- Combined ---


def get_genome_summary(
    tax_id: int,
    max_records,
    max_bvbrc_records=0,
    sources=None,
    progress_callback=None,
):
    """
    Fetches NCBI + Ensembl genome data, returns unified list of genome metadata dicts.
    sources: list of "NCBI" and/or "Ensembl". If None, fetch both.
    """
    if sources is None:
        sources = ["NCBI", "Ensembl", "ENA", "BV-BRC"]

    summary = {"tax_id": tax_id, "organism": "Unknown", "genomes": []}

    # Fetch all sources concurrently
    fetch_funcs = {
        "NCBI": lambda: fetch_ncbi_genomes(
            tax_id, max_records, progress_callback=progress_callback
        ),
        "Ensembl": lambda: fetch_ensembl_genomes(tax_id),
        "ENA": lambda: fetch_ena_genomes(tax_id),
        "BV-BRC": lambda: fetch_bvbrc_genomes(
            tax_id, max_records=max_bvbrc_records
        ),
    }

    with ThreadPoolExecutor(max_workers=len(sources)) as executor:
        futures = {source: executor.submit(fetch_funcs[source]) for source in sources if source in fetch_funcs}
        results = {source: future.result() for source, future in futures.items()}

    # Process results
    if "NCBI" in results and results["NCBI"]:
        reports = results["NCBI"].get("reports", [])
        for r in reports:
            metadata = extract_ncbi_genome_metadata(r)
            if metadata:
                summary["genomes"].append(metadata)

    if "Ensembl" in results and results["Ensembl"] and isinstance(results["Ensembl"], list):
        for g in results["Ensembl"]:
            metadata = extract_ensembl_genome_metadata(g)
            if metadata:
                summary["genomes"].append(metadata)

    if "ENA" in results and results["ENA"] and isinstance(results["ENA"], list):
        for g in results["ENA"]:
            metadata = extract_ena_genome_metadata(g)
            if metadata:
                summary["genomes"].append(metadata)

    if "BV-BRC" in results and results["BV-BRC"] and isinstance(results["BV-BRC"], list):
        for g in results["BV-BRC"]:
            meta = extract_bvbrc_genome_metadata(g)
            if meta:
                summary["genomes"].append(meta)

    return summary
