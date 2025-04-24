import pytest
import pandas as pd
from unittest.mock import patch, Mock

from app_uniprot_id_to_genomic_loci_to_guide_generation import extract_genomic_information_from_uniprot_id  # Replace with actual module name

# Mock response data
mock_response_data = [{
    "accession": "P12345",
    "name": "MockProtein",
    "taxid": 9606,
    "sequence": "MADEUPSEQUENCE",
    "gnCoordinate": {
        "genomicLocation": {
            "chromosome": "1",
            "start": 123456,
            "end": 123999,
            "reverseStrand": False,
            "nucleotideId": "NC_000001.11",
            "assemblyName": "GRCh38",
            "exon": [
                {
                    "id": 1,
                    "proteinLocation": {
                        "begin": {"position": 1},
                        "end": {"position": 20}
                    },
                    "genomeLocation": {
                        "begin": {"position": 123456},
                        "end": {"position": 123475}
                    }
                },
                {
                    "id": 2,
                    "proteinLocation": {
                        "begin": {"position": 21},
                        "end": {"position": 40}
                    },
                    "genomeLocation": {
                        "begin": {"position": 123800},
                        "end": {"position": 123820}
                    }
                }
            ]
        },
        "ensemblGeneId": "ENSG00000123456",
        "ensemblTranscriptId": "ENST00000123456",
        "ensemblTranslationId": "ENSP00000123456"
    }
}]


@patch("app_uniprot_id_to_genomic_loci_to_guide_generation.requests.get")  # Patch where the function is defined
def test_extract_genomic_information_from_uniprot_id(mock_get):
    mock_response = Mock()
    mock_response.status_code = 200
    mock_response.json.return_value = mock_response_data
    mock_get.return_value = mock_response

    uniprot_id = "P12345"
    df = extract_genomic_information_from_uniprot_id(uniprot_id)

    # Assertions
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "accession" in df.columns
    assert df.loc[0, "accession"] == "P12345"
    assert df.loc[0, "gnCoordinate.ensemblGeneId"] == "ENSG00000123456"
    assert df.loc[0, "exon_id"] == "1,2"
    assert df.loc[0, "exon_proteinLocation.begin.position"] == "1,21"
    assert df.loc[0, "exon_genomeLocation.begin.position"] == "123456,123800"
