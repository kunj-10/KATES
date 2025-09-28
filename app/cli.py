import typer

app = typer.Typer()


@app.command()
def hello(name: str):
    """Say hello"""
    typer.echo(f"Hello {name}!")

@app.command()  
def goodbbye(name: str):
    """GoodBye"""
    typer.echo("goodbye")

if __name__ == "__main__":
    app()
