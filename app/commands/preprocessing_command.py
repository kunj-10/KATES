import typer

from app.services import PreprocessingService

preprocessing_app = typer.Typer()


@preprocessing_app.command()
def load(file_path: str, output_fasta: str = "final.fasta"):
    """
    Takes the file path to the .tar.gz files and get fasta files from them
    """
    PreprocessingService.make_fasta_file(file_path, output_fasta)
