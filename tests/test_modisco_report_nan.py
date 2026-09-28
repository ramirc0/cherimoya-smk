"""modisco 2.5.2 report vs pandas 3: a pattern with fewer TomTom matches than
n_matches crashes unless pandas' string inference is off (modisco_report.py)."""

import pandas
import pytest

MEME = "MEME version 4\n\nALPHABET= ACGT\n\nMOTIF ZBTB_0\nletter-probability matrix: alength= 4 w= 2\n" \
    "1 0 0 0\n0 1 0 0\n\n"


def _report_logos(tmp_path):
    # The match table and dict as generate_descriptive_report builds them.
    from modiscolite.descriptive_report import create_tomtom_match_logos

    db = tmp_path / "db.meme"
    db.write_text(MEME)
    # Other patterns fill match1, so the column holds strings and the gap.
    df = pandas.DataFrame({"match0": ["ZBTB_0", "ZBTB_0"], "qval0": [0.04, 0.01],
        "match1": [None, "ZBTB_0"], "qval1": [None, 0.02]})
    data = {"pos_patterns.pattern_0": {f"match_{j}": df.iloc[0][f"match{j}"] for j in range(2)}}
    return create_tomtom_match_logos(data, str(tmp_path), str(db), 2)


def test_missing_match_crashes_with_string_inference(tmp_path):
    with pytest.raises(AttributeError, match="strip"):
        _report_logos(tmp_path)


def test_missing_match_skipped_without_string_inference(tmp_path):
    with pandas.option_context("future.infer_string", False):
        logos = _report_logos(tmp_path)
    assert list(logos["pos_patterns.pattern_0"]) == ["match_0_logo", "match_0_base64"]
