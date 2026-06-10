# CORe Designer

CORe Designer offers a protein-centric approach for the generation of guide RNAs (gRNAs) for CRISPR mutagenesis studies. It also supports the generation of homology-directed repair (HDR) templates and integration-specific primer designs, following the **CORe methodology** developed by the **Child Lab at Imperial College London**.

This repository contains all the code required to run CORe Designer from the terminal and package it for distribution.

---

## 🚀 Quick Start

For most users, CORe Designer is packaged into an easy-to-use **Graphical User Interface (GUI)**. You can access the application directly via the releases sidebar on this Github.

### GUI Features & Navigation
The interface is split into two main operational windows depending on your input data:

* **Uniprot Locus Tab:** Use this tab if you are working with a protein that already has genomic coordinates associated with it on UniProt.
* **Custom Locus Tab:** Use this tab if your target protein does not have pre-associated genomic coordinates on UniProt.

---

## 🧪 Testing

All tests for this project are written using the `pytest` framework. To run the test suite locally, navigate to the test directory and run the console command:

```bash
cd code/tests/
pytest
