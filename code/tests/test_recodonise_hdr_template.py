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
        downstream_distance
    )

    # --- Check if recodonisation preserved sequence length ---
    assert len(recodonised_template) == len(hdr_template), "Length mismatch after recodonisation."

    # --- Check that the recodonisation check status matches expected ---
    assert recodonisation_check == expected_reconstruction_status, f"Recodonisation check failed: {recodonisation_check} != {expected_reconstruction_status}"

    # --- Optional: Check that input HDR sequence was modified (if expected) ---
    if hdr_template.upper() not in ["ATG", "TAA"]:
        assert recodonised_template != hdr_template, "Recodonisation did not modify the template when expected."

