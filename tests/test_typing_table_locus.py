"""typing.table: the frozen contig table carries its locus into typing.

The 15-column screen table has always had SEQUENCE/START/END/STRAND, but the
typing boundary used to drop them. The antigen engine needs the locus to keep
one physical flagellin gene out of both phase families; without it a
cross-reacting `fliC` allele wins the `fljB` family and the real phase is lost
(real genome: `9,12:l,v:1,7` came back as `9,12:l,v:-`). These tests pin the
parse so the table and the JSON screen document type identically.
"""

from __future__ import annotations

import pytest

from gspec.errors import InputError
from gspec.typing.table import parse_table

HEADER = (
    "#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS"
    "\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
)
ROW = (
    "genome.fna\tctg1\t466021\t467523\t-\th1_l_v__22\t1-1503/1503\t==============="
    "\t0/0\t100.00\t99.73\tsen_sero\tupstream:l,v\tantigen lv\tantigen\n"
)


def test_contig_table_carries_the_locus() -> None:
    """Given the frozen 15-column table, When it is parsed for typing, Then
    every row carries the contig/interval and strand the phase guard needs."""
    (row,) = parse_table("screen.tsv", HEADER + ROW)
    assert (row.contig, row.start, row.end, row.strand) == ("ctg1", 466021, 467523, "-")


def test_table_without_locus_columns_keeps_the_locus_unknown() -> None:
    """A hand-written table with only the required columns still types: the
    locus stays unknown instead of being invented."""
    text = (
        "FILE\tGENE\tDATABASE\t%COVERAGE\t%IDENTITY\n"
        "genome.fna\th1_l_v__22\tsen_sero\t100.00\t99.73\n"
    )
    (row,) = parse_table("screen.tsv", text)
    assert (row.contig, row.start, row.end, row.strand) == (None, None, None, None)


def test_blank_locus_cells_keep_the_locus_unknown() -> None:
    """A row whose coordinate cells are empty has no locus, not coordinate 0."""
    text = (
        "FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tDATABASE\t%COVERAGE\t%IDENTITY\n"
        "genome.fna\tctg1\t\t\t-\th1_l_v__22\tsen_sero\t100.00\t99.73\n"
    )
    (row,) = parse_table("screen.tsv", text)
    assert (row.contig, row.start, row.end) == (None, None, None)


def test_non_numeric_coordinates_are_rejected() -> None:
    """Coordinates are parsed, never silently ignored: a malformed START is a
    typed input error like any other malformed table cell."""
    text = (
        "FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tDATABASE\t%COVERAGE\t%IDENTITY\n"
        "genome.fna\tctg1\tstart\t1600\t-\th1_l_v__22\tsen_sero\t100.00\t99.73\n"
    )
    with pytest.raises(InputError) as excinfo:
        parse_table("screen.tsv", text)
    assert excinfo.value.code == "SCREEN_TABLE_MALFORMED"
