import os
import subprocess
import tarfile

import typer


class PreprocessingService:
    @staticmethod
    def extract_tarfile(file_path: str, extract_dir: str = "extracted_files") -> str:
        """Extracts a .tar.gz archive into a folder"""
        os.makedirs(extract_dir, exist_ok=True)
        with tarfile.open(file_path, "r:gz") as tar:
            tar.extractall(path=extract_dir)

        typer.echo(f"Tar file: {file_path} extracted to {extract_dir}")
        return extract_dir

    @staticmethod
    def make_fasta_file(file_path: str, output_fasta: str = "sequences.fasta"):
        """Convert extracted BLAST database to FASTA using blastdbcmd"""
        extract_dir = PreprocessingService.extract_tarfile(file_path)

        blast_db = None
        for file in os.listdir(extract_dir):
            if file.endswith(".nsq"):
                blast_db = os.path.splitext(file)[0]
                break

        if not blast_db:
            typer.echo("No BLAST database found in the archive.")
            return

        db_path = os.path.join(extract_dir, blast_db)
        typer.echo(f"Found BLAST database: {db_path}")

        try:
            subprocess.run(
                [
                    "blastdbcmd",
                    "-db",
                    db_path,
                    "-entry",
                    "all",
                    "-out",
                    os.path.join(extract_dir, output_fasta),
                ],
                check=True,
            )
            typer.echo(f"FASTA file created: {output_fasta}")

        except FileNotFoundError:
            typer.echo(
                "Error: blastdbcmd not found. Please install NCBI BLAST+ and ensure it’s in PATH."
            )
        except subprocess.CalledProcessError as e:
            typer.echo(f"Error while running blastdbcmd: {e}")
