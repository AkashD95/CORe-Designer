import pytest
from app_utils import reverse_complement

@pytest.mark.parametrize(
    "input_seq,expected_output",
    [
        ("ATCG", "CGAT"),
        ("atcg", "cgat"),
        ("AaTtCcGg", "cCgGaAtT"),  # preserves case
        ("", ""),
        ("ATXB", "BXAT"),  # unknown characters passed through
        ("GAATTC", "GAATTC"),  # palindrome
        ("NNNN", "NNNN"),  # ambiguous bases unchanged
    ]
)
def test_reverse_complement(input_seq, expected_output):
    assert reverse_complement(input_seq) == expected_output
