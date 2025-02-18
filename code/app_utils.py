#Utilty functions
from Bio import SeqIO

def clean_loci(loci):
    """
    Cleans and processes loci of interest by removing introns, appending upstream and downstream sequences, and ensuring uppercase format.

    Arguments:
        loci (str) - input loci DNA sequence
        
    Returns:
        loci - cleaned loci sequence 
        loci_exon - exon sequence extracted from loci (assuming introns are coded as lower case)     
    """
    #Remove any leading or trailing blank spaces from loci
    loci = loci.strip()
    
    #Extract exons 
    loci_exon = ""
    for base in loci:
        if base.isupper():
            loci_exon += base

    return loci, loci_exon

def reverse_complement(dna_seq):
    """ 
    Return the reverse complement of a DNA sequence. 
    Arguments:
        dna_seq - DNA input sequence as a string
        
    Returns:
        reverse_complement - Reverse complement of DNA input sequence 
    
    """
    complement = {'a': 't', 'A':'T', 't': 'a', 'T':'A','c': 'g','C':'G', 'g': 'c', 'G':'C'}
    
    reverse_complement_seq = ''.join(complement.get(base, base) for base in reversed(dna_seq))
        
    return reverse_complement_seq

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
    
    # Initialize the protein sequence
    protein = ""

    # Iterate over the DNA sequence in steps of 3 (codon length)
    for i in range(0, len(dna_seq) - len(dna_seq) % 3, 3):
        
        # Define a codon as 3 bases
        codon = dna_seq[i:i + 3]
        # Translate the codon into an amino acid based on the dictionary
        if codon.isupper():
            protein += codon_table.get(codon, '?')  # '?' is a placeholder for unknown codons
        else:
            protein += 'x'

    return protein

def extract_genome_multifast_to_list_of_sequence_records(genome_multi_fasta):
    '''
    Takes in genomic multi fasta input and list of sequence records out 

    Args:
    genome_multi_fasta(fasta): multi fasta containing sequence records for a genome

    Returns:
    genome_sequence_records(list): list of sequences in biopython seqRecord format
    '''
    print('Starting genome extraction...')
    genome_sequence_records = []
    for record in SeqIO.parse(genome_multi_fasta, "fasta"):
        genome_sequence_records.append(record)
    print('Converted to list of sequence records')
    
    return genome_sequence_records