import pytest
from app_uniprot_id_to_genomic_loci_to_guide_generation import get_GC_content  # Replace 'your_module' with your actual module name

@pytest.mark.parametrize(
    "sequence, expected_gc",
    [
        ("ATATAT", 0.0),
        ("GGGG", 100.0),
        ("CCCC", 100.0),
        ("GCGCGC", 100.0),
        ("ATGC", 50.0),
        ("atgc", 50.0),
        ("AaTtGgCc", 50.0),
        ("A1T2G3C4", 50.0),
        ("!@#ATGC$%^", 50.0),
    ]
)
def test_get_GC_content_valid(sequence, expected_gc):
    assert get_GC_content(sequence) == expected_gc


@pytest.mark.parametrize("invalid_input", [
    "",    # Empty string
    None,  # None input
])
def test_get_GC_content_invalid_input_raises_assertion(invalid_input):
    with pytest.raises(AssertionError):
        get_GC_content(invalid_input)