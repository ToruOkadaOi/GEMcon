import typer
from typing import Optional
from .flow import main_flow

app = typer.Typer()


@app.command()
def main(
    branch: str = typer.Option(..., "--branch", help="[transcriptomic] or [proteomic] Branch"), # typer.Argument() insteade?
    input: Optional[str] = typer.Option(None, "--input", help="Input file eg - .loom/.h5ad"),
    task: Optional[str] = typer.Option(None, "--task", help="Task type eg - annotate/metabolic"),
    algo: Optional[str] = typer.Option(
        None, 
        "--algo", 
        help="Algorithm for metabolic modeling (default - gimme for transcriptomic, geckopy for proteomic; Choices - [gimme, tinit, fastcore(in beta), imat(in beta)])")
):
    if branch not in ["transcriptomic", "proteomic"]:
        raise typer.BadParameter("branch must be 'transcriptomic' or 'proteomic'")
    
    if task and task not in ["annotate", "metabolic"]:
        raise typer.BadParameter("task must be 'annotate' or 'metabolic'")
    
    if branch == "transcriptomic" and task is None:
        raise typer.BadParameter("--task is required when --branch is transcriptomic")

    main_flow(branch, input, task, algo)


if __name__ == "__main__":
    app()