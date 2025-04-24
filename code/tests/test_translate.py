import pytest
from app_utils import translate  # Replace with the actual name of your module

@pytest.mark.parametrize(
    "dna_seq, expected_protein",
    [
        ("ATGGCC", "MA"),                     # Valid codons: Met, Ala
        ("atggcc", "xx"),                     # Lowercase introns → should yield 'x'
        ("ATGatgGCC", "MxA"),                 # Mixed case: second codon lowercase
        ("ATGTAA", "M*"),                     # Contains stop codon
        ("ATGXXX", "M?"),                     # Unknown codon (XXX is not valid)
        ("", ""),                             # Empty string → empty protein
        ("ATGGC", "M"),                       # Incomplete final codon is ignored
        ("GCTGCCGCA", "AAA"),                 # Synonymous codons for Alanine
        ("tgcTGT", "xC"),                     # First lowercase → 'x', second uppercase → Cysteine
        ("AAA", "K"),                         # Single codon for Lysine
        ("ATGTAGTGA", "M**"),                 # Multiple stop codons
    ]
)
def test_translate(dna_seq, expected_protein):
    assert translate(dna_seq) == expected_protein