import primer3

def design_integration_specific_primers(hdr_template_recodonised, loci_plus_recodonised_hdr_template, upstream_homology_arm_end, upstream_guide_distance, downstream_homology_arm_start, downstream_guide_distance):
    '''
    This function will design integration specific primers to target redocodonised HDR templates. These primer pairs will be used to generated amplicons for illumina sequencing

    Args:
    hdr_template_recodonised (str): hdr template recodonised to be synonymous to the wild type
    loci_plus_recodonised_hdr_template (str): Full genomic loci with recodonised hdr template included
    upstream_homology_arm_end (int): numerical position of end of upstream homology arm in the loci 
    upstream_guide_distance (int): numerical distance between end of the upstream homology arm and the start of the codon for the amino acid of interest
    downstream_homology_arm_start (int): numerical position of start of the downstream homology arm in the loci 
    downstream_guide_distance (int): numerical distance between start of the upstream homology arm and the end of the codon for the amino acid of interest
    
    Returns:
    integration_specific_primers (DataFrame): dataframe containing forward and reverse primer sequences and Tm's
    '''
    
    #Make sure integration-specific primer falls within the exonic portion of the recodonised region
    exon_start_index = upstream_homology_arm_end
    upstream_primer_range = upstream_guide_distance
    if loci_plus_recodonised_hdr_template[exon_start_index].islower():
        first_upper = next((i for i, char in enumerate(hdr_template_recodonised) if char.isupper()), -1)
        exon_start_index = first_upper + upstream_homology_arm_end
        upstream_primer_range = upstream_guide_distance + first_upper

    exon_end_index = downstream_homology_arm_start
    downstream_primer_range = downstream_guide_distance
    if loci_plus_recodonised_hdr_template[exon_end_index].islower():
        last_upper = next((i for i, char in enumerate(hdr_template_recodonised[::-1]) if char.isupper()), -1)
        exon_end_index = downstream_homology_arm_start - last_upper
        downstream_primer_range = downstream_guide_distance - last_upper

    print("exon start index: " + str(exon_start_index))
    print("exon end index: " + str(exon_end_index))

    #Add logic to make sure that all recodonised region specific primers at at most 50 bp from the codon of interest
    upstream_recodonised_primer_offset = 0
    if abs(upstream_primer_range) > 50:
        upstream_recodonised_primer_offset = abs(upstream_primer_range) - 50

    downstream_recodonised_primer_offset = 0
    if downstream_primer_range > 50:
        downstream_recodonised_primer_offset = downstream_primer_range - 50

    
    try: 
        if(len(loci_plus_recodonised_hdr_template) - exon_end_index < 1000):
            downstream_primer_search_length = len(loci_plus_recodonised_hdr_template) - exon_end_index
        else:
            downstream_primer_search_length = 1000

        print("fwd primer start: " + str(exon_start_index + upstream_recodonised_primer_offset))
        print("fwd primer search distance: " + str(abs(upstream_primer_range) - 5 - upstream_recodonised_primer_offset))
        print("rev primer start: " + str(exon_end_index))
        print("rev primer search distance: " + str(downstream_primer_search_length))

        primers = primer3.bindings.design_primers(
        seq_args= {
            'SEQUENCE_TEMPLATE': str(loci_plus_recodonised_hdr_template),
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [
                exon_start_index + upstream_recodonised_primer_offset,
                abs(upstream_primer_range) - 5 - upstream_recodonised_primer_offset,
                exon_end_index,
                downstream_primer_search_length
        ]
        },
        global_args={
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 15,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 58.0,
            'PRIMER_MAX_TM': 62.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 500.0,
            'PRIMER_PRODUCT_SIZE_RANGE': [[450,550]]
        })

        #print(primers)
        # If no primers are returned, raise an exception
        if primers['PRIMER_LEFT_NUM_RETURNED'] == 0:
            raise ValueError("No primers found")
        
        left_primer_sequence = primers['PRIMER_LEFT'][0]['SEQUENCE']
        print(left_primer_sequence)
        left_primer_tm = primers['PRIMER_LEFT'][0]['TM']
        print(left_primer_tm)
        left_primer_start = primers['PRIMER_LEFT'][0]['COORDS'][0]
        left_primer_len = primers['PRIMER_LEFT'][0]['COORDS'][1]
        
        right_primer_sequence = primers['PRIMER_RIGHT'][0]['SEQUENCE']
        print(right_primer_sequence)
        right_primer_tm = primers['PRIMER_RIGHT'][0]['TM']
        print(right_primer_tm)
        right_primer_start = primers['PRIMER_RIGHT'][0]['COORDS'][0]
        right_primer_len = primers['PRIMER_RIGHT'][0]['COORDS'][1]
    
        return left_primer_sequence, left_primer_tm, right_primer_sequence, right_primer_tm
    
    except ValueError as e:
        print(e)
        # Handle the case where no primers are found and continue by making overlap primer on the right side
        if(exon_start_index < 1000):
            upstream_search_position = 0
            upstream_primer_search_length = exon_start_index
        else:
            upstream_search_position = exon_start_index - 1000
            upstream_primer_search_length = 1000

        print("fwd primer start: " + str(upstream_search_position))
        print("fwd primer search distance: " + str(upstream_primer_search_length))
        print("rev primer start: " + str(exon_end_index - downstream_primer_range + 5))
        print("rev primer search distance: " + str(downstream_primer_range - 5 - downstream_recodonised_primer_offset))

        primers = primer3.bindings.design_primers(
        seq_args= {
            'SEQUENCE_TEMPLATE': str(loci_plus_recodonised_hdr_template),
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [
                upstream_search_position,
                upstream_primer_search_length,
                exon_end_index - downstream_primer_range + 5,
                downstream_primer_range - 5 - downstream_recodonised_primer_offset
        ]
        },
        global_args={
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 15,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 58.0,
            'PRIMER_MAX_TM': 62.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 500.0,
            'PRIMER_PRODUCT_SIZE_RANGE': [[450,550]]
        })

        #print(primers)
        # If no primers are returned, raise an exception
        if primers['PRIMER_LEFT_NUM_RETURNED'] == 0:
            raise ValueError("No integration-specific primer found after target codon")

        left_primer_sequence = primers['PRIMER_LEFT'][0]['SEQUENCE']
        print(left_primer_sequence)
        left_primer_tm = primers['PRIMER_LEFT'][0]['TM']
        print(left_primer_tm)
        left_primer_start = primers['PRIMER_LEFT'][0]['COORDS'][0]
        left_primer_len = primers['PRIMER_LEFT'][0]['COORDS'][1]
        
        right_primer_sequence = primers['PRIMER_RIGHT'][0]['SEQUENCE']
        print(right_primer_sequence)
        right_primer_tm = primers['PRIMER_RIGHT'][0]['TM']
        print(right_primer_tm)
        right_primer_start = primers['PRIMER_RIGHT'][0]['COORDS'][0]
        right_primer_len = primers['PRIMER_RIGHT'][0]['COORDS'][1]

        return left_primer_sequence, left_primer_tm, right_primer_sequence, right_primer_tm