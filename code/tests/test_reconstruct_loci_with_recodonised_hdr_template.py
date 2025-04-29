import pytest

from app_selected_guides_to_hdr_template import reconstruct_loci_with_recodonised_hdr_template

@pytest.mark.parametrize(
    "loci, hdr_template_recodonised, upstream_end, downstream_start, expected_success",
    [
        # --- Passing Case ---
        ("AAAACCCCGGGGTTTT", "GGGG", 4, 8, True),  # Correct size HDR, matches region replaced

        # --- Failing Case ---
        ("AAAACCCCGGGGTTTT", "GGGGA", 4, 8, False),  # Extra nucleotide in HDR template → length mismatch
    ]
)
def test_reconstruct_loci_with_recodonised_hdr_template(loci, hdr_template_recodonised, upstream_end, downstream_start, expected_success):
    """
    Test reconstructing loci with recodonised HDR template,
    checking for both successful and failed reconstructions based on sequence lengths.
    """
    result = reconstruct_loci_with_recodonised_hdr_template(
        loci,
        hdr_template_recodonised,
        upstream_end,
        downstream_start
    )

    if expected_success:
        assert result is not None, "Expected reconstruction to succeed, but it failed."
        loci_plus_hdr, section = result
        assert len(loci_plus_hdr) == len(loci), "Reconstructed loci length mismatch."
    else:
        assert result is None, "Expected reconstruction to fail, but it succeeded."