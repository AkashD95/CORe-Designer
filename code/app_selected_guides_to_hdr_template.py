from app_utils import translate, clean_loci

def generate_hdr_template(loci, upstream_guide_position, downstream_guide_position, homology_overlap):
    """
    Finds the sequence between two guide RNA sequences in a genomic loci and adds an homology overlap.
    
    Arguments:
    loci (str): Genomic loci that the guides were generated from.
    upstream_guide_position: position of the start of the upstream guide PAM
    downstream_guide_position: position of the start of the downstream guide PAM
    homology_overlap (int): Number of nucleotides to pad around the sequence found between the guides. Default = 100

    Returns:
    hdr_template (str): Sequence found between the start of the PAM sites for the gRNA cut sites.
    hdr_template_plus_homology (str): Sequence found between the gRNAs, padded by homology overlap.
    selected_guides_distance_from_aa: distance from target amino acid for selected guides. Upstream guide is distance to start of the codon. Downstream guide is distance to end of the codon.
    """
     
    # Get the start and end indices in the loci of the region between guides
    start_hdr_template = int(upstream_guide_position)
    print(start_hdr_template)
    end_hdr_template = int(downstream_guide_position)
    print(end_hdr_template)

        
    # Calculate the indices for the start and end of homology arms
    upstream_homology_arm_start = max(0, start_hdr_template - homology_overlap)
    upstream_homology_arm_end = upstream_guide_position
    downstream_homology_arm_start = downstream_guide_position
    downstream_homology_arm_end = min(len(loci), end_hdr_template + homology_overlap)
        
    # Extract the homology arms
    hdr_template = loci[start_hdr_template:end_hdr_template] 
    print(hdr_template)
    upstream_homology_arm = loci[upstream_homology_arm_start:upstream_homology_arm_end]
    print(upstream_homology_arm)
    downstream_homology_arm = loci[downstream_homology_arm_start:downstream_homology_arm_end+1]
    print(downstream_homology_arm)
    
    return hdr_template, upstream_homology_arm,upstream_homology_arm_start, upstream_homology_arm_end, downstream_homology_arm, downstream_homology_arm_start, downstream_homology_arm_end

def recodonise_hdr_template(hdr_template, upstream_guide_distance, downstream_guide_distance, codon_usage_table):
    """
    Recodonise hdr template to codon with nearest amino acid frequency.
    Arguments: 
    hdr_template(str): hdr template sequence to be used used as the template for recodonisation
    upstream_guide_distance (int): distance from start of upstream guide RNA to the start of the codon for the amino acid of interest
    downstream_guide_distance (int): distance from the start of the downstream guide RNA to the end of the codon for the amino acid of interest
    codon_usage_table (dict): a python dictionary of codon as key and frequency as value
    
    Returns: 
    recodonised_hdr_template(str): recodonised template sequence to be used base template for hdr library generation
    spec_aa_codon_start (int): 
    recodonisation_check (str)

    """

    # 2 dictionaries to be used by the function
    # codon_dict pairs each codon to the corresponding amino acid
    codon_dict = {
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 
    'AAT': 'N', 'AAC': 'N', 
    'GAT': 'D', 'GAC': 'D',
    'TGT': 'C', 'TGC': 'C',
    'GAA': 'E', 'GAG': 'E',
    'CAA': 'Q', 'CAG': 'Q',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'CAT': 'H', 'CAC': 'H',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'AAA': 'K', 'AAG': 'K',
    'ATG': 'M',  
    'TTT': 'F', 'TTC': 'F',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'TGG': 'W',
    'TAT': 'Y', 'TAC': 'Y',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TAA': '*', '*': '*', 'TGA': '*'
    }

    #aa_to_codon pairs each amino acid to a codon
    aa_to_codon = {
    'I': ['ATA', 'ATC', 'ATT'], # Isoleucine
    'M': ['ATG'], # Methionine
    'T': ['ACA', 'ACC', 'ACG', 'ACT'], #Threonine
    'N': ['AAC', 'AAT'], # Asparginine 
    'K': ['AAA', 'AAG'], #Lysine
    'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],  # Serine
    'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],  # Arginine
    'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], # Leucine
    'P': ['CCA', 'CCC', 'CCG', 'CCT'], # Proline
    'H': ['CAC', 'CAT'], # Histidine 
    'Q': ['CAA', 'CAG'], # Gluatmine
    'V': ['GTA', 'GTC', 'GTG', 'GTT'], # Valine
    'A': ['GCA', 'GCC', 'GCG', 'GCT'], # Alanine
    'D': ['GAC', 'GAT'], # Aspartic Acid
    'E': ['GAA', 'GAG'], # Glutamic Acid 
    'G': ['GGA', 'GGC', 'GGG', 'GGT'], # Glycine
    'F': ['TTC', 'TTT'], #Phenylalanine
    'Y': ['TAC', 'TAT'], #Tyrosine
    'C': ['TGC', 'TGT'], #Cysteine
    'W': ['TGG'], #Tryptophan
    '*': ['TAA', 'TAG', 'TGA']  # Stop codons
}
    
    #First check that the hdr_template_exon is in frame and shift start point of recodonisation to be in frame
    #this is done by identifying the target codon and frameshifting to that reference
    hdr_template, hdr_template_exon = clean_loci(hdr_template)
    spec_aa_codon_start = 0 + abs(upstream_guide_distance)
    spec_aa_codon_end = len(hdr_template) - downstream_guide_distance
    
    spec_aa_codon = hdr_template[spec_aa_codon_start:spec_aa_codon_end]
    print(f'hdr template(length:{len(hdr_template)}):')
    print(hdr_template)
    
    print(f'hdr template exon (length:{len(hdr_template_exon)}):')
    print(hdr_template_exon)
    
    print(f'Target codon is at position: {spec_aa_codon_start} - {spec_aa_codon_end - 1}')
    print(spec_aa_codon)
    
    print('Translation:')
    print(f'{translate(spec_aa_codon)}')

    #Correct spec_aa_codon_start and end for possible introns and the start and the end
    exon_start_index = next((i for i, char in enumerate(hdr_template) if char.isupper()), -1)
    spec_aa_codon_start_exon_corrected = spec_aa_codon_start - exon_start_index

    
    if spec_aa_codon_start_exon_corrected % 3 == 0:
        print('Exons are in frame')
        n = 0
    elif spec_aa_codon_start_exon_corrected % 3 == 1:
        print('Exons are out of frame by 1 position. Shifting start of recodonisation by 1')
        n = 1
    elif spec_aa_codon_start_exon_corrected % 3 == 2:
        print('Exons are out of frame by 2 positions. Shifting start of recodonisation by 2')
        n = 2
        
    #changing codons to aa
    aa_sequence = ""
    hdr_template_exon_dict = {}
    exon_start_index = next((i for i, char in enumerate(hdr_template) if char.isupper()), -1)
    for i in range(n, len(hdr_template_exon)-2, 3):
        codon_start = i
        codon_end = i+3
        codon = hdr_template_exon[codon_start:codon_end]
        aa = codon_dict[codon]
        aa_sequence += aa
        hdr_template_exon_dict[i // 3] = codon
   
    print(f'sequence to be recodonised(length:{len(hdr_template[n::])}):')
    print(hdr_template[n::])
    print('aa sequence')
    print(aa_sequence)
    
    recodonised_sequence = ""
    for m, aa in enumerate(aa_sequence):
        current_codon = hdr_template_exon_dict[m]
        print('Current codon:')
        print(current_codon)
        print('Translation:')
        print(translate(current_codon))
        

        if aa in ['M', "W", "*"]:
            print('M, W  or * will not be recodonised')
            recodonised_sequence += current_codon
            print('recodonised_sequence:')
            print(recodonised_sequence)
            
        else:
            aa_dict = {codon : codon_usage_table[codon] for codon in aa_to_codon[aa]}
            codon_freq_list = sorted(aa_dict.items(), key=lambda x: x[1], reverse=True)
            print('Codon frequency list:')
            print(codon_freq_list)
            for i, codon_freq in enumerate(codon_freq_list):
                new_codon = ''
                if (codon_freq[0] == current_codon) and (i == 0):
                    new_codon = codon_freq_list[1][0] 
                    print('Recodonising down a rank. New codon:')
                    print(new_codon)
                    recodonised_sequence += new_codon # lower down a 
                    print('recodonised_sequence:')
                    print(recodonised_sequence)

                elif codon_freq[0] == current_codon:
                    new_codon = codon_freq_list[i-1][0]
                    print('Recodonising up a rank. New codon:')
                    print(new_codon)
                    recodonised_sequence += new_codon  # up a rank
                    print('recodonised_sequence:')
                    print(recodonised_sequence)
                    
    print(f'recodonised sequence (length:{len(recodonised_sequence)})')
    print(recodonised_sequence)
    
    print('translating recodonised sequence')
    recodonised_sequence_translation = translate(recodonised_sequence)
    print(recodonised_sequence_translation)
    
    assert aa_sequence == recodonised_sequence_translation, "Recodonisation failed"
    print("Recodonised translation matches original translation")
    
    # Now reattach the recodonised hdr region back to the hdr template
    # find region upstream to reattach to recodonised sequence. +n + 1 correct for the earlier frameshift correction
    exon_start_index = next((i for i, char in enumerate(hdr_template) if char.isupper()), -1)
    upstream_hdr_template = hdr_template[0:exon_start_index+n]
    print('Upstream hdr template')
    print(upstream_hdr_template)
    downstream_hdr_template = hdr_template[exon_start_index+n+len(recodonised_sequence)::]
    print('Downstream hdr template')
    print(downstream_hdr_template)
    
    hdr_template_recodonised = upstream_hdr_template + recodonised_sequence + downstream_hdr_template
    print(f'recodonised hdr template(length:{len(hdr_template_recodonised)})')
    print(hdr_template_recodonised)
    
    if len(hdr_template) == len(hdr_template_recodonised):
        print('Length of hdr template = length of recodonised hdr template')
        recodonisation_check = 'Pass'
    else:
        print('hdr template reconstruction failed')
        recodonisation_check = 'Fail'
        
    # return the recodonized sequence
    return hdr_template_recodonised, spec_aa_codon_start, recodonisation_check

def reconstruct_loci_with_recodonised_hdr_template(loci, hdr_template_recodonised, upstream_homology_arm_end, downstream_homology_arm_start):
    '''
    This function integrates the recodonised hdr template back into the loci. 

    Args:
    loci (str): Genomic loci that the guides were generated from. 
    hdr_template_recodonised (str): hdr template recodonised to be synonymouse to the wild type 
    upstream_homology_arm_end (int): position in the loci of the end of the upstream homology arm
    downstream_homology_arm_start (int): position in the loci of the start of the downstream homology arm

    Returns:
    section_loci_plus_recodonised_hdr_template (str): Section of the loci with the recodonised hdr template integrated. Used for visual inspection and automated primer design
    '''
    upstream_loci = loci[0:int(upstream_homology_arm_end)]
    downstream_loci = loci[int(downstream_homology_arm_start)::]
    
    loci_plus_recodonised_hdr_template = upstream_loci + hdr_template_recodonised + downstream_loci
    print('Full loci with recodonised HDR template integrated:')
    print(loci_plus_recodonised_hdr_template)

    if len(loci) == len(loci_plus_recodonised_hdr_template):
        print(f'Reconstruction successful -> loci ({len(loci)}) equals loci with hdr template integrated ({len(loci_plus_recodonised_hdr_template)})')
        
        section_loci_plus_recodonised_hdr_template = loci_plus_recodonised_hdr_template[int(upstream_homology_arm_end) - 1000:int(downstream_homology_arm_start) + 1000]
        print('Section of loci with recodonised HDR template integrated:')
        print(section_loci_plus_recodonised_hdr_template)
        
        return loci_plus_recodonised_hdr_template, section_loci_plus_recodonised_hdr_template
    else:
        print(f'Reconstruction failed -> loci ({len(loci)}) not equal to loci with hdr template integrated ({len(loci_plus_recodonised_hdr_template)})')
        return
    
def generate_hdr_library(recodonised_hdr_template, spec_aa_codon_start, upstream_homology_arm, downstream_homology_arm, variant_specific = None, codon_list = [('I', 'ATC'), ('M', 'ATG'), ('T', 'ACC'), ('N', 'AAC'), ('K', 'AAG'), ('S', 'AGC'), ('R', 'CGG'), ('L', 'CTG'), ('P', 'CCC'), ('H', 'CAC'), ('Q', 'CAG'), ('V', 'GTG'), ('A', 'GCC'), ('D', 'GAC'), ('E', 'GAG'), ('G', 'GGC'), ('F', 'TTC'), ('Y', 'TAC'), ('C', 'TGC'), ('W', 'TGG'), ('*', 'TAA')]):
        '''
        Replaces a codon at a given index with codons of the highest frequency in a given organism.
        Arguments:
        hdr_template (str) - homologous directed repair template 
        spec_aa_codon_start (int) - start index of codon to be replaced 
        list_of_codons for library (str) - list of codons for library. Default is max frequencies for all 20 amino acids + 1 for stop codon
        variant_specific (str) = nucleotides for variant specific substitution 
        upstream_homology_arm  
        downstream_homology_arm
        
        Returns:
        mutated_sequence (str) - sequence with new codon mutated in
        '''
        print(f'Recodonised template (length:{len(recodonised_hdr_template)}')
        print(recodonised_hdr_template)
        
        #mutate sequence on index position
        mutated_sequences = []
        current_codon = recodonised_hdr_template[spec_aa_codon_start:spec_aa_codon_start+3]
        print(f'Codon to be mutated {current_codon} at position {spec_aa_codon_start}-{spec_aa_codon_start+2}')
        current_amino_acid = translate(current_codon)
        print(current_amino_acid)
        
        for aa_codon in codon_list:
            print(aa_codon)
            new_codon = str(aa_codon[1])
            print(new_codon)
            mutated_sequence = recodonised_hdr_template[:spec_aa_codon_start] + new_codon + recodonised_hdr_template[spec_aa_codon_start+3:]
            print(f'mutated_sequence (length:{len(mutated_sequence)})')
            print(mutated_sequence)
            if len(mutated_sequence) == len(recodonised_hdr_template):
                mutated_sequences.append(mutated_sequence)
        
        # Add variant specific codon if it's suppled and not already in the codon list  
        if variant_specific is not None and variant_specific not in codon_list:
            mutated_sequence = recodonised_hdr_template[:spec_aa_codon_start] + variant_specific + recodonised_hdr_template[spec_aa_codon_start+3:]
            mutated_sequences.append(mutated_sequence)

        print(f'hdr template library (size:{len(mutated_sequences)})')
        print(mutated_sequences)
        
        if spec_aa_codon_start % 3 == 0:
            print('Exons are in frame')
            n = 0
        elif spec_aa_codon_start % 3 == 1:
            print('Exons are out of frame by 1 position. Shifting start of recodonisation by 1')
            n = 1
        elif spec_aa_codon_start % 3 == 2:
            print('Exons are out of frame by 2 positions. Shifting start of recodonisation by 2')
            n = 2

        print('hdr template library (amino_acids)')
        
        aa_mutated_sequences = [] 
        for i in mutated_sequences:
            aa_mutated_sequence = translate(i[n::])
            aa_mutated_sequences.append(aa_mutated_sequence)

        print('Translation of mutated sequences:')      
        print(aa_mutated_sequences)

        #add homology arms 
        mutated_sequences_plus_homology_arms = [f"{upstream_homology_arm}{sequence}{downstream_homology_arm}" for sequence in mutated_sequences]


        return mutated_sequences_plus_homology_arms