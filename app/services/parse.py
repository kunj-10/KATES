import csv
import re
from collections import Counter
from itertools import product
from typing import Dict, List, Optional

from Bio import SeqIO

# Pre-compute k-mers for k=4 (256 total combinations)
K = 4
BASES = ["A", "T", "C", "G"]
KMERS = ["".join(p) for p in product(BASES, repeat=K)]


def kmer_count_record(seq: str, k: int = K, kmers: List[str] = KMERS) -> List[int]:
    """
    Calculate k-mer counts for a single sequence.

    Args:
        seq: DNA sequence string
        k: k-mer length (default: 4)
        kmers: pre-computed list of all possible k-mers

    Returns:
        List of k-mer counts in the same order as kmers list
    """
    seq = str(seq).upper()
    # Remove any non-ATCG characters and replace with empty string
    seq = re.sub(r"[^ATCG]", "", seq)

    if len(seq) < k:
        # If sequence is shorter than k, return all zeros
        return [0] * len(kmers)

    # Count k-mers in the sequence
    counts = Counter([seq[i : i + k] for i in range(len(seq) - k + 1)])

    # Return counts in the same order as kmers list
    return [counts.get(kmer, 0) for kmer in kmers]


def parse_fasta_description(description: str) -> Dict[str, Optional[str]]:
    """
    Parse FASTA header robustly.

    Handles cases like:
    - Multiple '>' in description (takes last one)
    - Missing strain information
    - Various strain formats (strain X, DSM X, ATCC X, NCFB X, or just numbers)
    - Gene type extraction

    Returns dict with: accession, species, strain, gene_type
    """
    description = description.strip()

    # Handle multiple '>' characters - take the last complete entry
    if ">" in description:
        # Split by '>' and take the last non-empty part
        parts = [part.strip() for part in description.split(">") if part.strip()]
        if parts:
            description = parts[-1]

    # Split into accession and rest
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

    # Extract strain information with various patterns
    strain = None
    strain_patterns = [
        r"strain\s+([^\s]+(?:\s+[^\s]+)*?)(?=\s+\d+S|\s+ribosomal|\s+RNA|$)",  # strain NRRL 11412
        r"DSM\s+(\d+)",  # DSM 4304
        r"ATCC\s+(\d+)",  # ATCC 31181
        r"NCFB\s+(\d+)",  # NCFB 2751
        r"strain\s+([A-Za-z0-9\-_]+)",  # strain VC-16
    ]

    # strain_match = None
    full_strain_text = None

    for pattern in strain_patterns:
        match = re.search(pattern, rest, re.IGNORECASE)
        if match:
            if "strain" in pattern:
                full_strain_text = match.group(0)  # Full "strain XXX" text
                strain = match.group(1)
            else:
                full_strain_text = match.group(0)  # Full "DSM XXX" etc
                strain = full_strain_text
            # strain_match = match
            break

    # Extract species name
    species_end_markers = [
        r"\s+strain\s+",
        r"\s+DSM\s+",
        r"\s+ATCC\s+",
        r"\s+NCFB\s+",
        r"\s+\d+S\s+ribosomal",  # 16S ribosomal
        r"\s+ribosomal",
        r"\s+\d+$",  # Just a number at the end
    ]

    species_part = rest

    # Find the earliest occurrence of any end marker
    earliest_pos = len(rest)
    for marker_pattern in species_end_markers:
        match = re.search(marker_pattern, rest, re.IGNORECASE)
        if match:
            earliest_pos = min(earliest_pos, match.start())

    if earliest_pos < len(rest):
        species_part = rest[:earliest_pos].strip()

    # Clean up species name
    species = species_part.strip() if species_part.strip() else "Unknown species"

    # Extract gene type (everything from 16S onwards, or ribosomal RNA info)
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
            # Clean up any trailing punctuation or extra whitespace
            gene_type = re.sub(r"\s+", " ", gene_type).strip()
            break

    # If no specific gene type found, try to extract any remaining meaningful info
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
    """
    Read a FASTA file and write a CSV with columns: accession, species, and 256 k-mer features.
    Robust to internal `>` characters in descriptions.
    """
    records: list[Dict[str, any]] = []

    # First, read the FASTA safely using SeqIO
    with open(fasta_path, "r") as handle:
        # preprocess lines to ensure '>' at start of line signals new entry
        fasta_text = handle.read()
        # replace any '>' that are not at the start of a line with a placeholder
        fasta_text_clean = re.sub(r"(?<!\n)>", "_GT_", fasta_text)

    # Now parse with SeqIO from a string
    from io import StringIO

    handle = StringIO(fasta_text_clean)
    for record in SeqIO.parse(handle, "fasta"):
        # replace placeholder back to original
        record.description = record.description.replace("_GT_", ">")

        # Get k-mer features instead of raw sequence
        kmer_features = kmer_count_record(record.seq)

        # Parse metadata
        meta = parse_fasta_description(record.description)

        # Create record with k-mer features
        record_dict = {
            "accession": meta["accession"],
        }

        # Add k-mer features as separate columns
        for i, kmer in enumerate(KMERS):
            record_dict[f"{kmer}"] = kmer_features[i]

        record_dict["species"] = meta["species"]

        records.append(record_dict)

    # Prepare fieldnames: metadata + k-mer columns
    fieldnames = ["accession"] + [f"{kmer}" for kmer in KMERS] + ["species"]

    # Write CSV
    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for r in records:
            writer.writerow(r)

    print(f"CSV with k-mer features saved to {csv_path}")
    print(f"Total records: {len(records)}")
    print(f"Features per record: {len(KMERS)} k-mers + metadata")


def process_single_record(record, include_metadata: bool = True) -> Dict[str, any]:
    """
    Process a single SeqIO record and return k-mer features with optional metadata.

    Args:
        record: BioPython SeqRecord object
        include_metadata: Whether to include parsed metadata

    Returns:
        Dictionary with k-mer features and optional metadata
    """
    # Get k-mer features
    kmer_features = kmer_count_record(record.seq)

    result = {}

    # Add metadata if requested
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

    # Add k-mer features
    for i, kmer in enumerate(KMERS):
        result[f"kmer_{kmer}"] = kmer_features[i]

    return result


if __name__ == "__main__":
    # Example usage
    fasta_file = r"D:\Hackathon\SIH 2025\KATES\final.fasta"  # Path to your input FASTA
    csv_file = "output_kmers.csv"  # Path to your output CSV
    fasta_to_csv_with_kmers(fasta_file, csv_file)

    # # Example of processing a single record
    # print("\nExample of processing single record:")
    # with open(fasta_file, "r") as handle:
    #     fasta_text = handle.read()
    #     fasta_text_clean = re.sub(r"(?<!\n)>", "_GT_", fasta_text)

    # from io import StringIO

    # handle = StringIO(fasta_text_clean)

    # # Process first record as example
    # for record in SeqIO.parse(handle, "fasta"):
    #     record.description = record.description.replace("_GT_", ">")
    #     single_record_result = process_single_record(record)
    #     print(f"Accession: {single_record_result['accession']}")
    #     print(f"Species: {single_record_result['species']}")
    #     print(
    #         f"First 5 k-mer features: {[single_record_result[f'kmer_{kmer}'] for kmer in KMERS[:5]]}"
    #     )
    #     break
