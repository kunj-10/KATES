from typing import Dict, List, Optional
from collections import Counter
from itertools import product
from Bio import SeqIO
import csv
import re

K = 4
BASES = ["A", "T", "C", "G"]
KMERS = ["".join(p) for p in product(BASES, repeat=K)]

def kmer_count_record(seq: str, k: int = K, kmers: List[str] = KMERS) -> List[int]:
    seq = str(seq).upper()
    seq = re.sub(r"[^ATCG]", "", seq)

    if len(seq) < k:
        return [0] * len(kmers)

    counts = Counter([seq[i : i + k] for i in range(len(seq) - k + 1)])
    return [counts.get(kmer, 0) for kmer in kmers]


def parse_fasta_description(description: str) -> Dict[str, Optional[str]]:
    description = description.strip()

    if ">" in description:
        parts = [part.strip() for part in description.split(">") if part.strip()]
        if parts:
            description = parts[-1]

    parts = description.split(maxsplit=1)
    accession = parts[0] if parts else ""
    rest = parts[1] if len(parts) > 1 else ""

    if not rest:
        return {
            "accession": accession,
            "species": "Unknown species",
            "strain": None,
            "gene_type": None,
        }

    strain = None
    strain_patterns = [
        r"strain\s+([^\s]+(?:\s+[^\s]+)*?)(?=\s+\d+S|\s+ribosomal|\s+RNA|$)", 
        r"DSM\s+(\d+)",
        r"ATCC\s+(\d+)",
        r"NCFB\s+(\d+)", 
        r"strain\s+([A-Za-z0-9\-_]+)", 
    ]

    full_strain_text = None

    for pattern in strain_patterns:
        match = re.search(pattern, rest, re.IGNORECASE)
        if match:
            if "strain" in pattern:
                full_strain_text = match.group(0)  
                strain = match.group(1)
            else:
                full_strain_text = match.group(0)  
                strain = full_strain_text
            break

    species_end_markers = [
        r"\s+strain\s+",
        r"\s+DSM\s+",
        r"\s+ATCC\s+",
        r"\s+NCFB\s+",
        r"\s+\d+S\s+ribosomal",
        r"\s+ribosomal",
        r"\s+\d+$",
    ]

    species_part = rest

    earliest_pos = len(rest)
    for marker_pattern in species_end_markers:
        match = re.search(marker_pattern, rest, re.IGNORECASE)
        if match:
            earliest_pos = min(earliest_pos, match.start())

    if earliest_pos < len(rest):
        species_part = rest[:earliest_pos].strip()

    species = species_part.strip() if species_part.strip() else "Unknown species"

    gene_type_patterns = [
        r"(\d+S\s+ribosomal\s+RNA.*?)(?:\s*>|$)",
        r"(ribosomal\s+RNA.*?)(?:\s*>|$)",
        r"(\d+S.*?)(?:\s*>|$)",
    ]

    gene_type = None
    for pattern in gene_type_patterns:
        match = re.search(pattern, rest, re.IGNORECASE)
        if match:
            gene_type = match.group(1).strip()
            gene_type = re.sub(r"\s+", " ", gene_type).strip()
            break

    if not gene_type and "RNA" in rest:
        rna_match = re.search(r"RNA.*", rest, re.IGNORECASE)
        if rna_match:
            gene_type = rna_match.group(0).strip()

    return {
        "accession": accession,
        "species": species,
        "strain": strain,
        "gene_type": gene_type,
    }

def fasta_to_csv_with_kmers(fasta_path: str, csv_path: str):
    records: list[Dict[str, any]] = []

    with open(fasta_path, "r") as handle:
        fasta_text = handle.read()
        fasta_text_clean = re.sub(r"(?<!\n)>", "_GT_", fasta_text)

    from io import StringIO

    handle = StringIO(fasta_text_clean)
    for record in SeqIO.parse(handle, "fasta"):
        record.description = record.description.replace("_GT_", ">")
        kmer_features = kmer_count_record(record.seq)
        meta = parse_fasta_description(record.description)

        record_dict = {
            "accession": meta["accession"],
        }

        for i, kmer in enumerate(KMERS):
            record_dict[f"{kmer}"] = kmer_features[i]

        record_dict["species"] = meta["species"]
        records.append(record_dict)

    fieldnames = ["accession"] + [f"{kmer}" for kmer in KMERS] + ["species"]

    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for r in records:
            writer.writerow(r)

    print(f"CSV with k-mer features saved to {csv_path}")
    print(f"Total records: {len(records)}")
    print(f"Features per record: {len(KMERS)} k-mers + metadata")

def process_single_record(record, include_metadata: bool = True) -> Dict[str, any]:
    kmer_features = kmer_count_record(record.seq)

    result = {}

    if include_metadata:
        meta = parse_fasta_description(record.description)
        result.update(
            {
                "accession": meta["accession"],
                "species": meta["species"],
                "strain": meta["strain"],
                "gene_type": meta["gene_type"],
            }
        )

    for i, kmer in enumerate(KMERS):
        result[f"kmer_{kmer}"] = kmer_features[i]

    return result
