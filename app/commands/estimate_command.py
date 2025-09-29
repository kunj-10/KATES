from pathlib import Path

import typer
from rich.console import Console
from rich.panel import Panel

from app.services import HDBSCANService

console = Console()
estimate_app = typer.Typer()


@estimate_app.command()
def abundance(
    edna_fasta: str = typer.Argument(..., help="Path to eDNA FASTA file"),
    reference_fasta: str = typer.Argument(..., help="Path to reference FASTA file"),
    output_dir: str = typer.Option("results", help="Directory to save CSVs"),
    k: int = typer.Option(4, help="k-mer size for annotations"),
):
    """
    Compute sequence-level taxonomy, cluster sequences, and generate abundance CSVs.
    """
    service = HDBSCANService(k=k)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    console.rule("[bold green]Starting eDNA Abundance Pipeline")

    steps = [
        "Embedding sequences",
        "Clustering sequences with UMAP + HDBSCAN",
        "Building reference mapping (species/genus)",
        "Annotating sequences by k-mer overlap",
        "Building sequence-level taxonomy table",
        "Computing abundance per species, genus, cluster",
        "Saving results to CSVs",
    ]

    results = []

    for step in steps:
        with console.status(f"[cyan]{step}...", spinner="dots"):
            if step == "Embedding sequences":
                service.embed_fasta(edna_fasta)
                results.append(f"✅ {step} ({len(service.seq_ids)} sequences embedded)")
            elif step == "Clustering sequences with UMAP + HDBSCAN":
                df_clusters = service.cluster_sequences()
                results.append(
                    f"✅ {step} ({len(df_clusters['Cluster'].unique())} clusters found)"
                )
            elif step == "Building reference mapping (species/genus)":
                df_mapping = service.build_reference_mapping(reference_fasta)
                results.append(
                    f"✅ {step} ({len(df_mapping)} reference sequences mapped)"
                )
            elif step == "Annotating sequences by k-mer overlap":
                seq_annotation = service.annotate_sequences(
                    edna_fasta, reference_fasta, k=k
                )
                results.append(f"✅ {step} ({len(seq_annotation)} sequences annotated)")
            elif step == "Building sequence-level taxonomy table":
                df_seq_taxonomy = service.build_sequence_taxonomy(
                    df_mapping, df_clusters, seq_annotation
                )
                results.append(f"✅ {step} ({len(df_seq_taxonomy)} rows built)")
            elif step == "Computing abundance per species, genus, cluster":
                species_counts, genus_counts, cluster_counts = (
                    service.compute_abundance(df_seq_taxonomy)
                )
                results.append(
                    f"✅ {step} (Species: {len(species_counts)}, Genus: {len(genus_counts)})"
                )
            elif step == "Saving results to CSVs":
                service.export_results(output_dir)
                results.append(f"✅ {step} (All CSVs saved in '{output_dir}')")

    # Show all steps with ticks
    console.print(
        Panel(
            "\n".join(results), title="[bold green]Pipeline Complete ✅", expand=False
        )
    )
