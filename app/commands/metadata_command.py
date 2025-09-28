import typer

from app.services import MetadataService

metadata_app = typer.Typer()

@metadata_app.command()
def download(url: str, out_dir: str):
    """
    Downloads all the metadata.json files from NCBI database url and stores it in the out_dir directory.
    """
    MetadataService.download_metadata_jsons(url, out_dir)