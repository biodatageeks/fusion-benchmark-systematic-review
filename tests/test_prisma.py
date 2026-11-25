"""Figure 2 (PRISMA): the hard-coded counts must add up and agree with the benchmarks in the analysis."""
from __future__ import annotations

import plot_prisma_flow as pf

# Benchmarks found by backward and forward citation searching (Methods 2.1, plot_prisma_flow docstring)
CITATION_BENCHMARKS = {("Creason", 2021), ("Liu", 2016), ("Tan", 2015), ("Rodriguez-Martin", 2017), ("Jia", 2013),
                       ("Benelli", 2012), ("Nicorici", 2014)}
# Of these, also retrieved by the database queries (Supp. Table S2) but not selected from the database records
ALSO_IN_DATABASE = {("Liu", 2016), ("Tan", 2015), ("Rodriguez-Martin", 2017), ("Jia", 2013)}


def test_database_branch_adds_up():
    assert pf.PUBMED + pf.SCOPUS == 1449  # Supp. Table S2: total before deduplication
    assert pf.PUBMED + pf.SCOPUS - pf.DUPLICATES == pf.SCREENED
    assert pf.SCREENED - pf.SCREEN_EXCLUDED == pf.SOUGHT
    assert pf.SOUGHT - pf.NOT_RETRIEVED == pf.ASSESSED
    assert pf.ASSESSED - pf.FULLTEXT_EXCLUDED == pf.DATABASE_INCLUDED
    assert sum(n for _, n in pf.EXCLUSION_REASONS) == pf.FULLTEXT_EXCLUDED
    assert all(n > 0 for _, n in pf.EXCLUSION_REASONS)


def test_citation_branch_and_total(benchmarks):
    assert pf.CITATION_IDENTIFIED >= pf.CITATION_ASSESSED >= pf.CITATION_INCLUDED == len(CITATION_BENCHMARKS)
    assert pf.DATABASE_INCLUDED + pf.CITATION_INCLUDED == len(benchmarks) == 18
    present = set(zip(benchmarks["first_author"], benchmarks["year_x"]))
    assert CITATION_BENCHMARKS <= present
    assert len(present - CITATION_BENCHMARKS) == pf.DATABASE_INCLUDED
    assert ALSO_IN_DATABASE <= CITATION_BENCHMARKS and len(ALSO_IN_DATABASE) == pf.CITATION_ALSO_IN_DATABASE


def test_supplementary_table_s2_matches_figure():
    """Supp. Table S2 rows (supplement.tex) use the same numbers as Figure 2."""
    expected = [647, 802, 1449, 472, 977, 826, 151, 0, 151, 140, 11, 7, 4, 7, 7, 18]
    figure = [pf.PUBMED, pf.SCOPUS, pf.PUBMED + pf.SCOPUS, pf.DUPLICATES, pf.SCREENED, pf.SCREEN_EXCLUDED, pf.SOUGHT,
              pf.NOT_RETRIEVED, pf.ASSESSED, pf.FULLTEXT_EXCLUDED, pf.DATABASE_INCLUDED, pf.CITATION_IDENTIFIED,
              pf.CITATION_ALSO_IN_DATABASE, pf.CITATION_ASSESSED, pf.CITATION_INCLUDED,
              pf.DATABASE_INCLUDED + pf.CITATION_INCLUDED]
    assert figure == expected
