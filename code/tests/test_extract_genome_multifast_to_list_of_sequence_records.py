import pytest
from Bio.SeqRecord import SeqRecord
from app_utils import extract_genome_multifast_to_list_of_sequence_records  # Replace with your actual module

@pytest.mark.parametrize(
    "fasta_content, expected_ids, expected_seqs",
    [   #Case 1: Mock chromosomes
        (
            """>chr1
ATGCGTACGTAGCTAG
>chr2
GCTAGCTAGCTAACGT""",
            ["chr1", "chr2"],
            ["ATGCGTACGTAGCTAG", "GCTAGCTAGCTAACGT"]
        ),
    ]
)
def test_extract_genome_multifast_to_list_of_sequence_records(tmp_path, fasta_content, expected_ids, expected_seqs):
    # Create a temporary FASTA file
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    # Call the function with the temp file
    records = extract_genome_multifast_to_list_of_sequence_records(str(fasta_file))

    # Assertions
    assert isinstance(records, list)
    assert all(isinstance(record, SeqRecord) for record in records)
    assert [record.id for record in records] == expected_ids
    assert [str(record.seq) for record in records] == expected_seqs
