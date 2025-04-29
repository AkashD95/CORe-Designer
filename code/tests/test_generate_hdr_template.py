import pytest

from app_selected_guides_to_hdr_template import generate_hdr_template

@pytest.mark.parametrize(
    "loci, upstream_guide_position, downstream_guide_position, homology_overlap, expected_upstream_start, expected_downstream_end",
    [
        # --- Case 1: Normal case, guides somewhere in the middle ---
        (
            "ATGC" * 100,     # 400 bp
            100,              # upstream guide at position 100
            200,              # downstream guide at position 200
            20,               # homology overlap
            80,               # expected upstream homology start
            220              # expected downstream homology end 
        ),

        # --- Case 2: Upstream guide very close to beginning ---
        (
            "ATGC" * 50,      # 200 bp
            5,                # upstream guide at position 5 (very early)
            50,               # downstream guide at position 50
            10,               # homology overlap
            0,                # expected upstream homology start (can't go below 0)
            60                # expected downstream homology end
        ),

        # --- Case 3: Downstream guide very close to end ---
        (
            "ATGC" * 50,      # 200 bp
            100,              # upstream guide at 100
            195,              # downstream guide near the end
            20,               # homology overlap
            80,               # expected upstream homology start
            200               # expected downstream homology end (cannot exceed loci length)
        ),
    ]
)
def test_generate_hdr_template_cases(
    loci,
    upstream_guide_position,
    downstream_guide_position,
    homology_overlap,
    expected_upstream_start,
    expected_downstream_end
):
    hdr_template, upstream_homology_arm, upstream_homology_arm_start, upstream_homology_arm_end, downstream_homology_arm, downstream_homology_arm_start, downstream_homology_arm_end = generate_hdr_template(
        loci, upstream_guide_position, downstream_guide_position, homology_overlap
    )

    # --- Check that types are correct ---
    assert isinstance(hdr_template, str)
    assert isinstance(upstream_homology_arm, str)
    assert isinstance(downstream_homology_arm, str)
    assert isinstance(upstream_homology_arm_start, int)
    assert isinstance(upstream_homology_arm_end, int)
    assert isinstance(downstream_homology_arm_start, int)
    assert isinstance(downstream_homology_arm_end, int)

    # --- Check region extraction ---
    assert hdr_template == loci[upstream_guide_position:downstream_guide_position]

    # --- Check upstream homology arm extraction ---
    assert upstream_homology_arm == loci[expected_upstream_start:upstream_guide_position]

    # --- Check downstream homology arm extraction ---
    assert downstream_homology_arm == loci[downstream_guide_position:expected_downstream_end+1]


    # --- Check boundaries ---
    assert upstream_homology_arm_start == expected_upstream_start
    assert downstream_homology_arm_end == expected_downstream_end
