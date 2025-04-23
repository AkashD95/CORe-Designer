import pytest
from app_utils import clean_loci

@pytest.mark.parametrize(
    "input_loci, expected_loci, expected_exon",
    [
        ("  ATGcgtACGT  ", "ATGcgtACGT", "ATGACGT"),     # Mixed case with spaces
        ("acgtacgt", "acgtacgt", ""),                    # All lowercase (introns only)
        ("ACGTACGT", "ACGTACGT", "ACGTACGT"),            # All uppercase (exons only)
        ("", "", ""),                                    # Empty string
        ("aTgC", "aTgC", "TC"),                           # Single lowercase at front
        (" tTgGcC ", "tTgGcC", "TGC"),                   # Mixed and whitespace
        ("ATGcatCAT", "ATGcatCAT", "ATGCAT"),            # Mixed introns in the middle
        ("   ACGTacgtACGT   ", "ACGTacgtACGT", "ACGTACGT") # Padded whitespace both sides
    ]
)
def test_clean_loci(input_loci, expected_loci, expected_exon):
    cleaned_loci, loci_exon = clean_loci(input_loci)
    assert cleaned_loci == expected_loci
    assert loci_exon == expected_exon
