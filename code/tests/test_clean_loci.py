import pytest
from app_utils import clean_loci

@pytest.mark.parametrize(
    "input_loci, expected_loci, expected_exon",
    [   
        # Case 1: Mixed case with spaces
        ("  ATGcgtACGT  ", "ATGcgtACGT", "ATGACGT"),

        # Case 2: all introns     
        ("acgtacgt", "acgtacgt", ""),

        # Case 3: all exons
        ("ACGTACGT", "ACGTACGT", "ACGTACGT"),
                    
        # Case 4: all empty
        ("", "", ""),                                 

    ]
)
def test_clean_loci(input_loci, expected_loci, expected_exon):
    cleaned_loci, loci_exon = clean_loci(input_loci)
    assert cleaned_loci == expected_loci
    assert loci_exon == expected_exon
