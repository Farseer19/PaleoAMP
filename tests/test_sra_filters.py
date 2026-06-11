"""Tests for SRA search filters."""

import pandas as pd
import pytest

from paleoamp.data.sra import filter_for_ancient_indicators, filter_for_microbial


def _make_df(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


class TestFilterForMicrobial:
    def test_passes_metagenomic_with_microbial_keyword(self):
        df = _make_df([{
            "library_source": "METAGENOMIC",
            "study_title": "Ancient permafrost microbiome study",
        }])
        result = filter_for_microbial(df)
        assert len(result) == 1

    def test_drops_non_metagenomic_library_source(self):
        df = _make_df([{
            "library_source": "GENOMIC",
            "study_title": "Ancient microbial metagenome permafrost",
        }])
        result = filter_for_microbial(df)
        assert len(result) == 0

    def test_drops_metagenomic_without_microbial_keyword(self):
        df = _make_df([{
            "library_source": "METAGENOMIC",
            "study_title": "Ancient arctic plant and animal eDNA shotgun sequencing",
        }])
        result = filter_for_microbial(df)
        assert len(result) == 0

    def test_drops_excluded_plant_animal_study(self):
        df = _make_df([{
            "library_source": "METAGENOMIC",
            "study_title": "Ancient oral microbiome plant and animal eDNA",
        }])
        result = filter_for_microbial(df)
        assert len(result) == 0

    def test_passes_dental_calculus(self):
        df = _make_df([{
            "library_source": "METAGENOMIC",
            "study_title": "Medieval dental calculus metagenome",
        }])
        result = filter_for_microbial(df)
        assert len(result) == 1

    def test_passes_gut_metagenome(self):
        df = _make_df([{
            "library_source": "METAGENOMIC",
            "study_title": "Ancient human gut metagenome from archaeological site",
        }])
        result = filter_for_microbial(df)
        assert len(result) == 1

    def test_passes_archaea_keyword(self):
        df = _make_df([{
            "library_source": "METAGENOMIC",
            "study_title": "Permafrost archaea diversity ancient sediment",
        }])
        result = filter_for_microbial(df)
        assert len(result) == 1

    def test_keyword_match_across_columns(self):
        df = _make_df([{
            "library_source": "METAGENOMIC",
            "study_title": "Ancient arctic sediment sequencing",
            "environment (material)": "soil metagenome",
        }])
        result = filter_for_microbial(df)
        assert len(result) == 1

    def test_empty_dataframe(self):
        result = filter_for_microbial(pd.DataFrame())
        assert result.empty

    def test_multiple_rows_mixed(self):
        df = _make_df([
            {
                "library_source": "METAGENOMIC",
                "study_title": "Ancient oral microbiome from medieval burial",
            },
            {
                "library_source": "METAGENOMIC",
                "study_title": "Ancient plant and animal eDNA arctic",
            },
            {
                "library_source": "AMPLICON",
                "study_title": "Ancient gut bacteria 16S survey",
            },
            {
                "library_source": "METAGENOMIC",
                "study_title": "Permafrost soil metagenome bacteria ancient",
            },
        ])
        result = filter_for_microbial(df)
        assert len(result) == 2
        assert set(result["study_title"]) == {
            "Ancient oral microbiome from medieval burial",
            "Permafrost soil metagenome bacteria ancient",
        }

    def test_library_source_case_insensitive(self):
        df = _make_df([{
            "library_source": "metagenomic",
            "study_title": "Ancient oral microbiome study",
        }])
        result = filter_for_microbial(df)
        assert len(result) == 1

    def test_no_library_source_column(self):
        df = _make_df([{
            "study_title": "Ancient oral microbiome study",
        }])
        # No library_source column — criterion 1 is skipped, falls through to keyword check
        result = filter_for_microbial(df)
        assert len(result) == 1


class TestFilterForAncientIndicators:
    def test_passes_ancient_keyword(self):
        df = _make_df([{"study_title": "Ancient permafrost microbiome"}])
        assert len(filter_for_ancient_indicators(df)) == 1

    def test_passes_permafrost(self):
        df = _make_df([{"study_title": "Permafrost microbial diversity study"}])
        assert len(filter_for_ancient_indicators(df)) == 1

    def test_passes_coprolite(self):
        df = _make_df([{"study_title": "Microbiome of a medieval coprolite"}])
        assert len(filter_for_ancient_indicators(df)) == 1

    def test_passes_keyword_in_organism_name(self):
        # keyword in organism_name column, not study_title
        df = _make_df([{"study_title": "Shotgun metagenomics", "organism_name": "ancient metagenome"}])
        assert len(filter_for_ancient_indicators(df)) == 1

    def test_drops_modern_sample(self):
        df = _make_df([{"study_title": "Modern soil microbiome survey"}])
        assert len(filter_for_ancient_indicators(df)) == 0

    def test_sediment_no_longer_passes(self):
        df = _make_df([{"study_title": "Marine sediment metagenome survey"}])
        assert len(filter_for_ancient_indicators(df)) == 0

    def test_edna_no_longer_passes(self):
        df = _make_df([{"study_title": "Freshwater eDNA fish diversity"}])
        assert len(filter_for_ancient_indicators(df)) == 0

    def test_empty_dataframe(self):
        assert filter_for_ancient_indicators(pd.DataFrame()).empty

    def test_case_insensitive(self):
        df = _make_df([{"study_title": "ANCIENT SIBERIAN PERMAFROST METAGENOME"}])
        assert len(filter_for_ancient_indicators(df)) == 1
