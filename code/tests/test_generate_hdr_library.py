import pytest
from app_selected_guides_to_hdr_template import generate_hdr_library


# Helper function to mock the translation function used in the original code
def mock_translate(codon):
    codon_map = {
        'ATC': 'I', 'ATG': 'M', 'ACC': 'T', 'AAC': 'N', 'AAG': 'K', 'AGC': 'S',
        'CGG': 'R', 'CTG': 'L', 'CCC': 'P', 'CAC': 'H', 'CAG': 'Q', 'GTG': 'V',
        'GCC': 'A', 'GAC': 'D', 'GAG': 'E', 'GGC': 'G', 'TTC': 'F', 'TAC': 'Y',
        'TGC': 'C', 'TGG': 'W', 'TAA': '*'
    }
    return codon_map.get(codon, 'X')


# Parameterized test 1: Test basic functionality with different codons and homology arms
@pytest.mark.parametrize(
    "recodonised_hdr_template, spec_aa_codon_start, upstream_homology_arm, downstream_homology_arm, expected_mutated_sequences_length",
    [
        ('ATGACCAGC', 0, 'AAGCTG', 'CTGAAG', 21),  # 21 mutations in the library
        ('ATGACCAGC', 3, 'GTCGTA', 'AGTCAA', 21),  # Check with different homology arms
    ]
)
def test_generate_hdr_library_basic(recodonised_hdr_template, spec_aa_codon_start, upstream_homology_arm, downstream_homology_arm, expected_mutated_sequences_length):
    result = generate_hdr_library(recodonised_hdr_template, spec_aa_codon_start, upstream_homology_arm, downstream_homology_arm)

    # Assert correct number of mutated sequences generated
    assert len(result) == expected_mutated_sequences_length
    # Assert homology arms are correctly included in mutated sequences
    assert upstream_homology_arm in result[0]
    assert downstream_homology_arm in result[0]


# Parameterized test 2: Test with a specific variant codon
@pytest.mark.parametrize(
    "recodonised_hdr_template, spec_aa_codon_start, upstream_homology_arm, downstream_homology_arm, variant_specific, codon_to_check",
    [
        ('ATGACCAGC', 0, 'AAGCTG', 'CTGAAG', 'GTC', 'GTC'),  # Ensure variant codon GTC is included
        ('ATGACCAGC', 3, 'AAGCTG', 'CTGAAG', 'GGC', 'GGC'),  # Ensure variant codon GGC is included
    ]
)
def test_generate_hdr_library_with_variant_specific(recodonised_hdr_template, spec_aa_codon_start, upstream_homology_arm, downstream_homology_arm, variant_specific, codon_to_check):
    result = generate_hdr_library(recodonised_hdr_template, spec_aa_codon_start, upstream_homology_arm, downstream_homology_arm, variant_specific=variant_specific)

    # Check if the variant-specific codon is in any of the mutated sequences
    assert any(codon_to_check in seq for seq in result)


# Parameterized test 3: Test with out-of-frame codon start
@pytest.mark.parametrize(
    "recodonised_hdr_template, spec_aa_codon_start, upstream_homology_arm, downstream_homology_arm, expected_shift",
    [
        ('ATGACCAGC', 0, 'AAGCTG', 'CTGAAG', 0),  # No shift needed (in frame)
        ('ATGACCAGC', 1, 'AAGCTG', 'CTGAAG', 1),  # Shift by 1
        ('ATGACCAGC', 2, 'AAGCTG', 'CTGAAG', 2),  # Shift by 2
    ]
)
def test_generate_hdr_library_out_of_frame(recodonised_hdr_template, spec_aa_codon_start, upstream_homology_arm, downstream_homology_arm, expected_shift):
    result = generate_hdr_library(recodonised_hdr_template, spec_aa_codon_start, upstream_homology_arm, downstream_homology_arm)
    
    # Ensure the mutated sequence length is correct after shifting
    mutated_sequence = result[0]
    assert len(mutated_sequence) > 0  # Check non-empty result
    assert upstream_homology_arm in mutated_sequence  # Ensure homology arms are included
    assert downstream_homology_arm in mutated_sequence  # Ensure homology arms are included
