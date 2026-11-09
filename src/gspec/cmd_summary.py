"""The `gspec summary` command: matrix over abricate-format report tables.

Lives outside cli.py to keep that module small; cli.py registers it via
``register_summary_command``. The report table also arrives on stdin:
``gspec screen -d eco_dec *.fna | gspec summary`` summarizes the piped
output (a batched multi-FILE table is one input, so dutch mode applies),
and the explicit ``-`` marker reads stdin even under a terminal (until
EOF). A bare invocation at a terminal keeps the help-on-bare policy.
"""

from datetime import UTC, datetime
from pathlib import Path

import click

from gspec.dispatch import (
    EnumChoice,
    dispatch,
    quiet_option,
    stdin_is_tty,
    stdin_table,
)
from gspec.formats.summary import format_summary_tsv, render_summary_json, render_summary_md
from gspec.screening.contigs import OutputFormat
from gspec.summary import SummaryParams, build_summary


def _stdin_is_tty() -> bool:
    """Whether stdin is a terminal — the injectable probe behind the
    bare-help guard (a CliRunner stdin is never a tty)."""
    return stdin_is_tty()


@click.argument(
    "files",
    nargs=-1,
    required=False,
    type=click.Path(path_type=Path),
    help="Abricate-format report file(s) to summarize ('-' = stdin).",
)
@click.option(
    "--identity", "-i", is_flag=True, help="Cells show %IDENTITY values (default: +/- presence)."
)
@click.option(
    "--coverage", "-c", is_flag=True, help="Cells show %COVERAGE values (default: +/- presence)."
)
@click.option(
    "--nopath", "-p", is_flag=True, help="Basename row keys (FILE values / input filenames)."
)
@quiet_option()
@click.option(
    "--format",
    "-f",
    "output_format",
    type=EnumChoice(OutputFormat),
    default=OutputFormat.tsv,
    help="Output format.",
)
@click.pass_context
def summary_command(
    ctx: click.Context,
    files: tuple[Path, ...],
    identity: bool,
    coverage: bool,
    nopath: bool,
    quiet: bool,
    output_format: OutputFormat,
) -> None:
    """Summarize report table(s) into a gene presence/absence matrix."""
    args = list(files or [])
    # Bare invocation at a terminal keeps the help-on-bare policy; with a
    # pipe attached the table arrives on stdin. Runs before dispatch so the
    # help path exits cleanly — never through the error-envelope handler.
    if not args and _stdin_is_tty():
        click.echo(ctx.get_help())
        ctx.exit(2)

    def run() -> None:
        stdin_text = stdin_table(args, code="SUMMARY_MALFORMED", noun="report", verb="summarize")
        paths = args if stdin_text is None else []
        if identity and coverage and not quiet:
            click.echo("Using %IDENTITY/%COVERAGE per hit for the summary table", err=True)
        elif identity and not quiet:
            click.echo("Using %IDENTITY for the summary table instead of +/- presence", err=True)
        elif coverage and not quiet:
            click.echo("Using %COVERAGE for the summary table instead of +/- presence", err=True)
        params = SummaryParams(identity=identity, coverage=coverage, nopath=nopath)

        def warn(message: str) -> None:
            if not quiet:
                click.echo(f"WARNING: {message}", err=True)

        matrix = build_summary(paths, params, warn=warn, stdin_text=stdin_text)
        now = datetime.now(UTC)
        if output_format is OutputFormat.json:
            click.echo(render_summary_json(matrix, now=now), nl=False)
        elif output_format is OutputFormat.md:
            click.echo(render_summary_md(matrix, now=now), nl=False)
        else:
            click.echo(format_summary_tsv(matrix, csv=output_format is OutputFormat.csv), nl=False)

    dispatch(run)


def register_summary_command(app: click.Group) -> None:
    """Attach the summary command to the CLI app (bare invocation at a
    terminal prints help; piped stdin summarizes without arguments)."""
    app.command("summary")(summary_command)
