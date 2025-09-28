import typer

from app.commands import metadata_app

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

if __name__ == "__main__":
    app()
