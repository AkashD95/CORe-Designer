import pytest
from app_uniprot_id_to_genomic_loci_to_guide_generation import get_codon_index

# Dummy translate function (replace with actual if available)
def translate(dna_seq):
    """
    Generates a protein sequence from a DNA sequence inputted as a string. It will only convert if the exon is in upper
    case

    Arguments:
        dna_seq(str) - input DNA sequence
        
    Returns:
        protein_sequence(str) - output protein sequence 
    """
    
    #Generate dictionary of codons for dna to protein translation
    
    codon_table = {
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',  # Alanine
    'TGT': 'C', 'TGC': 'C',                          # Cysteine
    'GAT': 'D', 'GAC': 'D',                          # Aspartic acid
    'GAA': 'E', 'GAG': 'E',                          # Glutamic acid
    'TTT': 'F', 'TTC': 'F',                          # Phenylalanine
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',  # Glycine
    'CAT': 'H', 'CAC': 'H',                          # Histidine
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',              # Isoleucine
    'AAA': 'K', 'AAG': 'K',                          # Lysine
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',  # Leucine
    'ATG': 'M',                                      # Methionine (Start codon)
    'AAT': 'N', 'AAC': 'N',                          # Asparagine
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',  # Proline
    'CAA': 'Q', 'CAG': 'Q',                          # Glutamine
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',  # Arginine
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',  # Serine
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',  # Threonine
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',  # Valine
    'TGG': 'W',                                      # Tryptophan
    'TAT': 'Y', 'TAC': 'Y',                          # Tyrosine
    'TAA': '*', 'TAG': '*', 'TGA': '*'               # Stop codons
    }

# Unified test
@pytest.mark.parametrize(
    "loci, loci_exon, expected_aas, test_type",
    [
        # Case 1: valid cases
        ("aaaATGgggTTTcccGGTaaa", "ATGTTTGGT", ["M", "F", "G"], "valid"),

        # Case 2: invalid exon length
        ("aaaATGgggTTTcccGGaaa", "ATGTTTGG", None, "invalid_length"),
    
        # Case 3: no exon
        ("aaaatgccctttgggggt", "", None, "empty"),
        
    ]
)
def test_get_codon_index(loci, loci_exon, expected_aas, test_type):
    if test_type == "valid":
        result = get_codon_index(loci, loci_exon)
        assert list(result.keys()) == list(range(len(expected_aas)))
        assert [result[i]["Amino Acid"] for i in range(len(expected_aas))] == expected_aas
        for entry in result.values():
            assert all(base[1].isupper() for base in [entry["base_1"], entry["base_2"], entry["base_3"]])

    elif test_type in {"invalid_length", "empty"}:
        with pytest.raises(AssertionError):
            get_codon_index(loci, loci_exon)