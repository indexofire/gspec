"""The `gspec typing` command: designation from screen result tables.

Lives outside cli.py to keep that module small; cli.py registers it via
``register_typing_command``. The second stage of the designation pipeline:
``gspec screen -o result.tsv --db NAME`` detects genes (contig or reads
mode), then ``gspec typing result.tsv`` resolves NAME from the table's
DATABASE column and designates (the use-case lives in gspec.typing.ops).
The table also arrives on stdin: ``gspec screen 1.fna --db NAME | gspec
typing`` types the piped output, and the explicit ``-`` marker reads stdin
even under a terminal (until EOF).
"""

from pathlib import Path

import click

from gspec.dispatch import (
    EnumChoice,
    datadir_option,
    dispatch,
    output_target,
    quiet_option,
    stdin_is_tty,
    stdin_table,
)
from gspec.engines import OutputFormat
from gspec.errors import usage_fail
from gspec.typing.ops import run_typing


def _stdin_is_tty() -> bool:
    """Whether stdin is a terminal — the injectable probe behind the
    bare-help guard (a CliRunner stdin is never a tty)."""
    return stdin_is_tty()


@click.argument(
    "files",
    nargs=-1,
    required=False,
    type=click.Path(path_type=Path),
    help=(
        "gspec/abricate screen result TSV/CSV file(s) to type — contig"
        " (FILE/%IDENTITY/%COVERAGE) and reads (#SAMPLE/BREADTH%) tables are"
        " accepted, plus cluster JSON (gspec.cluster/2 or gspec.locus_map/1);"
        " '-' = stdin."
    ),
)
@datadir_option()
@click.option(
    "--format",
    "-f",
    "output_format",
    type=EnumChoice(OutputFormat),
    default=OutputFormat.tsv,
    help="Output format.",
)
@quiet_option()
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    help=(
        "Write the typing results to PATH instead of stdout (truncates; stdout"
        " then carries no data)."
    ),
)
@click.pass_context
def typing_command(
    ctx: click.Context,
    files: tuple[Path, ...],
    datadir: Path | None,
    output_format: OutputFormat,
    quiet: bool,
    output: Path | None,
) -> None:
    """Designate phenotypes from screen result tables or cluster JSON.

    Gene databases evaluate the table; cluster JSON is already designated by
    ``gspec screen`` and is rendered into the same concise typing result.
    """
    args = list(files or [])
    # Bare invocation at a terminal keeps the help-on-bare policy; with a
    # pipe attached the table arrives on stdin. Runs before dispatch so the
    # help path exits cleanly — never through the error-envelope handler.
    if not args and _stdin_is_tty():
        click.echo(ctx.get_help())
        ctx.exit(2)

    def run() -> None:
        if output_format is OutputFormat.csv:
            usage_fail("--format csv is not available for typing (use tsv, json, or md)")
        stdin_text = stdin_table(
            args, code="SCREEN_TABLE_MALFORMED", noun="screen result", verb="type"
        )
        paths = args if stdin_text is None else []
        with output_target(output) as deliver:
            # run_typing emits its output through `deliver`; echoing the
            # return value here would duplicate every byte.
            run_typing(paths, datadir, output_format, quiet, emit=deliver, stdin_text=stdin_text)

    dispatch(run)


def register_typing_command(app: click.Group) -> None:
    """Attach the typing command to the CLI app (bare invocation at a
    terminal prints help; piped stdin types without arguments)."""
    app.command("typing")(typing_command)
