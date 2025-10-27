import pytest
from app_selected_guides_to_hdr_template import recodonise_hdr_template

@pytest.mark.parametrize(
    "hdr_template, upstream_distance, downstream_distance, expected_reconstruction_status",
    [
        # 1. Simple case: Single codon HDR, in frame, no introns
        # Codon is ATG (start codon, M) -> Should stay the same because Methionine (M) must not change
        ("ATG", 0, 0, "Pass"),

        # 2. Recodonisation case: Alanine (GCT)
        # Alanine codon GCT, commonly preferred codon is GCC
        ("GCT", 0, 0, "Pass"),

        # 3. Multiple codons, simple exon only, codons should shift to similar ones
        ("GCTGCCGCG", 0, 0, "Pass"),
        
        # 4. Exons and introns mixed (lowercase introns)
        # exon: GCT intron: aaa exon: GCC
        ("GCTaaaGCC", 0, 0, "Pass"),

        # 5. Out of frame start - needs shifting
        ("aGCTGCCGCG", 1, 6, "Pass"),  # start 1 nucleotide off (intronic)

        # 6. Edge case: stop codon (TAA) - must not recodonise
        ("TAA", 0, 0, "Pass"),

        # 7. Complex case: Mixed codons, different preference frequencies
        ("GGTGGC", 0, 3, "Pass"),  # Glycine codons GGT (lower freq), GGC (higher freq)

        # 8. Frame shifted exon needing +2 correction
        ("aaGCTGCCGCG", 2, 6, "Pass"),  # shift by 2 to reach the frame
    ]
)
def test_recodonise_hdr_template(hdr_template, upstream_distance, downstream_distance, expected_reconstruction_status):
    """
    Tests recodonise_hdr_template under various scenarios:
    - Different codon sets
    - With and without introns
    - In-frame and out-of-frame scenarios
    """
    recodonised_template, spec_codon_start, recodonisation_check = recodonise_hdr_template(
        hdr_template,
        upstream_distance,
        downstream_distance,
        codon_usage_table = {'GGG': 669768, 'GGA': 669873, 'GGT': 437126, 'GGC': 903565, 'GAG': 1609975, 'GAA': 1177632, 'GAT': 885429, 'GAC': 1020595, 'GTG': 1143534, 'GTA': 287712, 'GTT': 448607, 'GTC': 588138, 'GCG': 299495, 'GCA': 643471, 'GCT': 750096, 'GCC': 1127679, 'AGG': 486463, 'AGA': 494682, 'AGT': 493429, 'AGC': 791383, 'AAG': 1295568, 'AAA': 993621, 'AAT': 689701, 'AAC': 776603, 'ATG': 896005, 'ATA': 304565, 'ATT': 650473, 'ATC': 846466, 'ACG': 246105, 'ACA': 614523, 'ACT': 533609, 'ACC': 768147, 'TGG': 535595, 'TGA': 63237, 'TGT': 430311, 'TGC': 513028, 'TAG': 32109, 'TAA': 40285, 'TAT': 495699, 'TAC': 622407, 'TTG': 525688, 'TTA': 311881, 'TTT': 714298, 'TTC': 824692, 'TCG': 179419, 'TCA': 496448, 'TCT': 618711, 'TCC': 718892, 'CGG': 464485, 'CGA': 250760, 'CGT': 184609, 'CGC': 423516, 'CAG': 1391973, 'CAA': 501911, 'CAT': 441711, 'CAC': 613713, 'CTG': 1611801, 'CTA': 290751, 'CTT': 536515, 'CTC': 796638, 'CCG': 281570, 'CCA': 688038, 'CCT': 713233, 'CCC': 804620}
    )

    # --- Check if recodonisation preserved sequence length ---
    assert len(recodonised_template) == len(hdr_template), "Length mismatch after recodonisation."

    # --- Check that the recodonisation check status matches expected ---
    assert recodonisation_check == expected_reconstruction_status, f"Recodonisation check failed: {recodonisation_check} != {expected_reconstruction_status}"

    # --- Optional: Check that input HDR sequence was modified (if expected) ---
    if hdr_template.upper() not in ["ATG", "TAA"]:
        assert recodonised_template != hdr_template, "Recodonisation did not modify the template when expected."

