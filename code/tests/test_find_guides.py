import pytest
from unittest.mock import patch
from app_uniprot_id_to_genomic_loci_to_guide_generation import find_guides

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

@pytest.mark.parametrize(
    "loci, pam, expected_output",
    [   
        # Case 1: guide on forward strand
        ("AACTGACGTACGTACGTACGTACCGG", "NGG", (["CGG"], [23], ["TGACGTACGTACGTACGTAC"], ["forward"])),
        
        # Case 2: guide on reverse strand
        ("CCGGTACGTACGTACGTACGTCAGTT", "NGG", (["CGG"], [3], ["TGACGTACGTACGTACGTAC"], ["reverse"])),
        # Case 3: no guides found
        ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "NGG", ([], [], [], []))
    ]
)
def test_find_guides(loci, pam, expected_output):
    with patch("app_uniprot_id_to_genomic_loci_to_guide_generation.reverse_complement", side_effect=reverse_complement):
        result = find_guides(loci, pam)
        assert result == expected_output
