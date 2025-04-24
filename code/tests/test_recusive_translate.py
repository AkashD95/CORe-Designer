import pytest
import pandas as pd
from unittest.mock import patch
from app_uniprot_id_to_genomic_loci_to_guide_generation import recursive_translate  # Replace with your actual module name


@pytest.mark.parametrize("genomic_loci_exons, mock_translation, expected_protein, expected_result, loci_exons_expected", [
    # Case 1: Translation matches the first ATG
    ("GGGATGAAATGCCCT", "MKCP", "MKCP", "Pass", "ATGAAATGCCCT"),

    # Case 2: Translation does not match
    ("GGGATGAAATGCCCTAG", "XYZ", "MA*", "Fail", None),

    # Case 3: No ATG present
    ("GGGCCCCTTT", "MA*", "MA*", "Fail", None),
])
@patch("app_utils.translate")
def test_recursive_translate_param(mock_translate, genomic_loci_exons, mock_translation, expected_protein, expected_result, loci_exons_expected):
    mock_translate.return_value = mock_translation
    df_input = pd.DataFrame({"sequence": [expected_protein]}).copy(deep=True)

    result = recursive_translate(genomic_loci_exons, df_input.copy())

    assert result["translation_check"].iloc[0] == expected_result
    if loci_exons_expected:
        assert result["loci_exons"].iloc[0] == loci_exons_expected
    else:
        assert "loci_exons" not in result.columns or pd.isna(result["loci_exons"].iloc[0])
