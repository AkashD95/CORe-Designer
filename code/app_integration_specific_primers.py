import primer3

def generate_integration_specific_primers(loci, hdr_template_recodonised, upstream_homology_arm_end_position, downstream_homology_arm_start_position):
    # Construct the recodonised loci
    upstream_loci = loci[0:upstream_homology_arm_end_position]
    downstream_loci = loci[downstream_homology_arm_start_position::]
    loci_plus_recodonised_hdr_template = upstream_loci + hdr_template_recodonised + downstream_loci
    if len(loci) == len(loci_plus_recodonised_hdr_template):
        print('Reconstruction successful')
        recodonised_loci_reconstruction_check = 'Pass'
    else:
        print('Reconstruction failed')
        recodonised_loci_reconstruction_check = 'Fail'

    # Run primer3 workflow
    seq_args = {} 
    seq_args['SEQUENCE_TEMPLATE'] = loci
    seq_args['SEQUENCE_FORCE_LEFT_END'] = upstream_homology_arm_end_position + 20

    try:    
        primers = primer3.bindings.design_primers(
        seq_args= seq_args,
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
            'PRIMER_PRODUCT_SIZE_RANGE': [
                [400,500]
            ],
        })
        print(primers)
        # If no primers are returned, raise an exception
        if primers['PRIMER_LEFT_NUM_RETURNED'] == 0:
            raise ValueError("No primers found")
        
        left_primer_sequence = primers['PRIMER_LEFT'][0]['SEQUENCE']
        left_primer_tm = primers['PRIMER_LEFT'][0]['TM']
        right_primer_sequence = primers['PRIMER_RIGHT'][0]['SEQUENCE']
        right_primer_tm = primers['PRIMER_RIGHT'][0]['TM']
        return(loci_plus_recodonised_hdr_template, recodonised_loci_reconstruction_check, left_primer_sequence, left_primer_tm, right_primer_sequence, right_primer_tm)
    
    except ValueError as e:
        print(e)
        # Handle the case where no primers are found and continue by making overlap primer on the right side
        seq_args = {} 
        seq_args['SEQUENCE_TEMPLATE'] = loci
        seq_args['SEQUENCE_FORCE_RIGHT_END'] = downstream_homology_arm_start_position - 20
        primers = primer3.bindings.design_primers(
        seq_args= seq_args,
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
            'PRIMER_PRODUCT_SIZE_RANGE': [
                [400,500]
            ],
        })
        print(primers)
        left_primer_sequence = primers['PRIMER_LEFT'][0]['SEQUENCE']
        left_primer_tm = primers['PRIMER_LEFT'][0]['TM']
        right_primer_sequence = primers['PRIMER_RIGHT'][0]['SEQUENCE']
        right_primer_tm = primers['PRIMER_RIGHT'][0]['TM']
        return(loci_plus_recodonised_hdr_template, recodonised_loci_reconstruction_check, left_primer_sequence, left_primer_tm, right_primer_sequence, right_primer_tm)