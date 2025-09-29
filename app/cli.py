import warnings

import typer

from app.commands import estimate_app, metadata_app, preprocessing_app

warnings.filterwarnings("ignore")

app = typer.Typer()


@app.command()
def hello(name: str):
    """Say hello"""
    typer.echo(f"Hello {name}!")


@app.command()
def goodbbye(name: str):
    """GoodBye"""
    typer.echo("goodbye")


app.add_typer(metadata_app, name="metadata")
app.add_typer(preprocessing_app, name="preprocess")
app.add_typer(estimate_app, name="estimate")

if __name__ == "__main__":
    app()
