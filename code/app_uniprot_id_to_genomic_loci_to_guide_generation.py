import pandas as pd
import requests
from app_utils import translate, reverse_complement, clean_loci

'''
This script will take in uniprot IDs and extract the genomic loci in format (introns in lowercase and exons in uppercase). This format is used by downstream processes to design guide RNAs,
homology directed repair templates and amplicon specific primers.
'''


def extract_genomic_information_from_uniprot_id(uniprot_id):
    '''
    Takes in uniprot ID as a string input and pings the Uniprot API to extract genomic coordinates of the protein and exons.
    Metadata such as name, taxID, protein sequence, genome assembly name,  ENSEMBL GeneID, ENSEMBL Transcript ID and ENSEMBL Translations IDs is included alongside the extracted coordinates.  

    Args:
    uniprot_id (str): uniprot ID

    Returns:
    genomic_information (pd.DataFrame): DataFrame containing genomic coordinates of the protein of interest alongside exon positions and metadata 
    '''
    genomic_information = pd.DataFrame()
    try:
        print(f'Searching for UniProt ID: {uniprot_id}')
        requestURL_protein = f"https://www.ebi.ac.uk/proteins/api/coordinates/{uniprot_id}"
        response_protein = requests.get(requestURL_protein, headers={"Accept": "application/json"})
        
        # Check if the request was successful
        response_protein.raise_for_status()
        
        # Load JSON response
        response_protein = response_protein.json()
        
        # Check if response is not empty
        if response_protein:
            response_protein_normalise = pd.json_normalize(
                response_protein, 
                record_path=['gnCoordinate', 'genomicLocation', 'exon'], 
                meta=['accession', 'name', 'taxid', 'sequence', 
                      ['gnCoordinate', 'genomicLocation', 'chromosome'], 
                      ['gnCoordinate', 'genomicLocation', 'start'], 
                      ['gnCoordinate', 'genomicLocation', 'end'], 
                      ['gnCoordinate', 'genomicLocation', 'reverseStrand'], 
                      ['gnCoordinate', 'genomicLocation', 'nucleotideId'], 
                      ['gnCoordinate', 'genomicLocation', 'assemblyName'], 
                      ['gnCoordinate', 'ensemblGeneId'], 
                      ['gnCoordinate', 'ensemblTranscriptId'], 
                      ['gnCoordinate', 'ensemblTranslationId']],
                record_prefix='exon_'
            )

            # Group and aggregate exon information
            response_protein_normalise = response_protein_normalise.groupby([
                'accession', 'name', 'taxid', 'sequence', 
                'gnCoordinate.genomicLocation.chromosome', 
                'gnCoordinate.genomicLocation.start', 
                'gnCoordinate.genomicLocation.end', 
                'gnCoordinate.genomicLocation.reverseStrand', 
                'gnCoordinate.genomicLocation.nucleotideId', 
                'gnCoordinate.genomicLocation.assemblyName', 
                'gnCoordinate.ensemblGeneId', 
                'gnCoordinate.ensemblTranscriptId', 
                'gnCoordinate.ensemblTranslationId'
            ]).agg({
                'exon_id': lambda x: ','.join(map(str, x)),
                'exon_proteinLocation.begin.position': lambda x: ','.join(map(str, x)),                    
                'exon_proteinLocation.end.position': lambda x: ','.join(map(str, x)),
                'exon_genomeLocation.begin.position': lambda x: ','.join(map(str, x)),                    
                'exon_genomeLocation.end.position': lambda x: ','.join(map(str, x))
            }).reset_index()

            # Concatenate to the main DataFrame
            genomic_information = pd.concat([genomic_information, response_protein_normalise], ignore_index=True)
        else:
            print(f"No data found for UniProt ID: {uniprot_id}")
            
    except Exception as e:
        print(f"An error occurred: {e}")

    return genomic_information

def recursive_translate(genomic_loci_exons, genomic_information_loci_extracted):
    # Initial check if "ATG" exists
    index = genomic_loci_exons.find("ATG")
    
    # If "ATG" is not found, return a fail
    if index == -1:
        print("ATG not found in the sequence")
        genomic_information_loci_extracted['translation_check'] = 'Fail'
        return genomic_information_loci_extracted

    # List to keep track of all ATG positions
    atg_positions = []
    while index != -1:
        atg_positions.append(index)
        # Search for the next "ATG" after the current index
        index = genomic_loci_exons.find("ATG", index + 1)

    # Now, iterate through each ATG position and try translating
    for atg_index in atg_positions:
        atg_start_sequence = genomic_loci_exons[atg_index:]
        print(f"Translating from ATG position {atg_index} /coding DNA length({len(atg_start_sequence)})...")
        print(atg_start_sequence)
        translation = translate(atg_start_sequence)
        print(f"Translation length ({len(translation)}): {translation}")
        # Check if translation matches the expected sequence
        if translation == genomic_information_loci_extracted['sequence'].iloc[0]:
            genomic_information_loci_extracted['translation_check'] = 'Pass'
            print("Translation matched. Marking as 'Pass'.")
            genomic_information_loci_extracted['loci_exons'] = atg_start_sequence
            return genomic_information_loci_extracted

        else:
            print('Translation match failed. Checking next index position.')
    # If no translation matches after checking all positions, mark as 'Fail'
    genomic_information_loci_extracted['translation_check'] = 'Fail'
    print("No matching translation found. Marking as 'Fail'.")
    return genomic_information_loci_extracted

def extract_genomic_loci_from_genomic_information(genomic_information, genome_sequence_records, UTR_length=1000):
    '''                  
    Takes the output from the Uniprot API from extract_genomic_information_from_uniprot_id function and extracts the genomic loci using a given reference genome. It appends the genomic loci, 
    genomic loci exons and performs a translation check between the protein sequence meta data outputted from extract_genomic_information_from_uniprot_id and the translation of the genomic loci exons
    outputted from this function 

    Args:
    genomic_information (pd.DataFrame): dataframe input from the output of extract_genomic_information_from_uniprot_id

    Returns:
    genomic_information_loci_extracted (pd.DataFrame): dataframe output with extracted loci, loci_exons and translation checks
    '''
    
    print(genomic_information)
    chromosome_id = str(genomic_information['gnCoordinate.genomicLocation.chromosome'].iloc[0])
    records_dict = {seq.id: seq for seq in genome_sequence_records}
    match = records_dict.get(chromosome_id)
    if match:
        is_reverse_strand = str(genomic_information['gnCoordinate.genomicLocation.reverseStrand'].iloc[0]) == 'True'
        uniprot_start = int(genomic_information['gnCoordinate.genomicLocation.start'].iloc[0])
        uniprot_end = int(genomic_information['gnCoordinate.genomicLocation.end'].iloc[0])
        
        # Adjust start and end for reverse strand
        if is_reverse_strand:
            uniprot_start, uniprot_end = uniprot_end, uniprot_start
           
        # Preparing exons DataFrame by splitting on commas
        df_exons_id = genomic_information[['exon_id']].copy()
        df_exons_id['exon_id'] = df_exons_id['exon_id'].str.split(',')
        
        df_exons_start = genomic_information[['exon_genomeLocation.begin.position']].copy()
        df_exons_start['exon_genomeLocation.begin.position'] = df_exons_start['exon_genomeLocation.begin.position'].str.split(',')
        
        df_exons_end = genomic_information[['exon_genomeLocation.end.position']].copy()
        df_exons_end['exon_genomeLocation.end.position'] = df_exons_end['exon_genomeLocation.end.position'].str.split(',')
        
        # Exploding columns
        df_exons_id = df_exons_id.explode('exon_id')
        df_exons_start = df_exons_start.explode('exon_genomeLocation.begin.position')
        df_exons_end = df_exons_end.explode('exon_genomeLocation.end.position')
        
        df_exons = pd.concat([df_exons_id, df_exons_start, df_exons_end], axis=1)
        df_exons.columns = ['exon_id', 'exon_start', 'exon_end']
        
        if is_reverse_strand:
            df_exons['exon_start'], df_exons['exon_end'] = df_exons['exon_end'], df_exons['exon_start']
        
        extracted_gene = match.seq[(uniprot_start - 1) - UTR_length:(uniprot_end - 1) + UTR_length].lower()
        extracted_gene = str(extracted_gene)
        print(f'Extracted gene (length:{len(extracted_gene)}):')
        print(extracted_gene)

        # Convert each exon sequence to uppercase in the extracted gene sequence
        for _, row in df_exons.iterrows():
            if row['exon_start'] != 'nan' and row['exon_end']!= 'nan':
                exon_start = int(float(row['exon_start'])) - uniprot_start + UTR_length
                print("Exon start:")
                print(exon_start)
                exon_end = int(float(row['exon_end'])) - uniprot_start + UTR_length
                print("Exon end:")
                print(exon_end)
            else:
                exon_start = exon_end = 0
                
            extracted_exon = extracted_gene[exon_start:exon_end+1]
            print('Extracted exon:')
            print(extracted_exon)
            
            extracted_gene = (
                    extracted_gene[:exon_start]
                    + extracted_exon.upper()
                    + extracted_gene[exon_end+1:]
                )
            
            print('Extracted gene:')
            print(extracted_gene)
        
        if is_reverse_strand:
            annotated_gene = reverse_complement(extracted_gene)
        else:
            annotated_gene = extracted_gene
        
        genomic_loci, genomic_loci_exons = clean_loci(annotated_gene)
        genomic_information_loci_extracted = genomic_information.assign(genomic_loci=genomic_loci, genomic_loci_exon=genomic_loci_exons)
        
        #Print reference loci and amino acid sequence 
        print(f'Genomic loci exons (length = {len(genomic_loci_exons)})')
        print(genomic_loci_exons)
        print(f"Reference amino acid sequence (length = {len(genomic_information_loci_extracted['sequence'].iloc[0])}):")
        print(genomic_information_loci_extracted['sequence'].iloc[0])

        # Translate the extracted loci exons and match them against the uniprot sequence 
        if genomic_loci_exons.startswith("ATG"):
            translation = translate(genomic_loci_exons)
            print('Translating exons...')
            print(translation)
            if translation == genomic_information_loci_extracted['sequence'].iloc[0]:
                genomic_information_loci_extracted = genomic_information_loci_extracted.assign(translation_check='Pass')
                print('Translation matches protein sequence:')
                return genomic_information_loci_extracted
            else:
                genomic_information_loci_extracted = genomic_information_loci_extracted.assign(translation_check='Fail')   
                print('Translation does not match protein sequence (running string search as possible fix):')
                return genomic_information_loci_extracted
        else:
            genomic_information_loci_extracted = recursive_translate(genomic_loci_exons, genomic_information_loci_extracted)
            return genomic_information_loci_extracted

# guide generation for specific amino acid
def find_guides(loci, pam):
    """ 
    Finds guide RNA sequences on the positive and negative strands of the loci by 
    searching for PAM sites and extracting 20 bases upstream from the PAM position. 
    
    Arguments:
        loci - DNA input sequence as a string
        pam - protospacer adjacent motif for Cas protein
        
    Returns:
        crRNA_DNA_start_position - position of crRNA_DNA start site (adjusted for python 0 indexing)
        crRNA_DNA_sequence - sequence 20bp upstream of PAM
        strand - which strand PAM was found on (positive or negative)
    """
    
    # Initialize lists to store gRNA information
    pam_found = []
    crRNA_DNA_sequence_start_position = []
    crRNA_DNA_sequence = []
    strand = []
    
    # Define PAM motifs
    if pam == "NGG":
        print(f'Searching forward strand for {pam} sites...')
        print(loci)
        for n in range(len(loci) - 2):  # Ensure n + 2 in the following line is within bounds
            if (loci[n] in ['A','a','T','t','C','c','G','g'] and loci[n+1] in ['G','g'] and loci[n+2] in ['G','g']):
                if n - 21 >= 0:  # Ensure there are enough bases before the PAM
                    pam_found.append(loci[n:n+3])
                    crRNA_DNA_sequence_start_position.append(n) # crRNA DNA sequence start site
                    crRNA_DNA_sequence.append(loci[n-20:n])  # Extract sequence (adjusted for python 0 indexing and python slicing being non inclusive of the final value)
                    strand.append("forward")
        
        print(f'Searching reverse strand for {pam} sites...')
        loci_reverse_complement = reverse_complement(loci)
        print(loci_reverse_complement)
        for n in range(len(loci_reverse_complement) - 2):  # Ensure n + 2 in the following line is within bounds
            if (loci_reverse_complement[n] in ['A','a','T','t','C','c','G','g'] and loci_reverse_complement[n+1] in ['G','g'] and loci_reverse_complement[n+2] in ['G','g']):
                if n - 21 >= 0:  # Ensure there are enough bases before the PAM
                    pam_found.append(loci_reverse_complement[n:n+3])
                    crRNA_DNA_sequence_start_position.append(len(loci) - (n))  #crRNA DNA sequence start site
                    crRNA_DNA_sequence.append(loci_reverse_complement[n-20:n])  # Extract sequence
                    strand.append("reverse")
                    
    else:
        raise ValueError("PAM not recognised")

    return pam_found, crRNA_DNA_sequence_start_position, crRNA_DNA_sequence, strand

#get codon index function
def get_codon_index(loci, loci_exon):
    """ 
    Maps the positions of the protein coding DNA sequences (exons) along the loci. It gives information 
    regarding the amino acid residue in each position.

    Args:
    loci (str): extracted DNA loci sequence 
    loci_exon (str): extracted DNA loci exon sequence

    Returns:
    protein_dict (dict): a dictionary where each base in the loci exon is mapped to a protein and its location in the loci
    """
    assert len(loci_exon) % 3 == 0, "The length of the coding sequence is not a multiple of 3" #Check that exon is a multiple of 3
    Base_df = pd.DataFrame(list(loci), columns = ["Base"]) #List of all the bases in the genomic loci
    Base_df["Position"] = Base_df.index
    Base_df = Base_df[Base_df['Base'].str.isupper()].reset_index(drop=True) # Only
    pos1 = Base_df[Base_df.index % 3 == 0].reset_index(drop=True) #Generate a list of every third base from the first base
    pos2 = Base_df[Base_df.index % 3 == 1].reset_index(drop=True) #Generate a list of every third base from the second base
    pos3 = Base_df[Base_df.index % 3 == 2].reset_index(drop=True) #Generate a list of every third base from the third base
    aas = list(translate(loci_exon))
    protein_dict = {}
    for aa in range(len(aas)):
        protein_dict[aa] = {"Amino Acid": aas[aa],
                        "base_1": (pos1["Position"][aa],
                                pos1["Base"][aa]),
                        "base_2": (pos2["Position"][aa],
                                pos2["Base"][aa]),
                        "base_3": (pos3["Position"][aa],
                                pos3["Base"][aa])}

    return protein_dict

def get_GC_content(grna_dna_sequence):
    
    """ 
    Return the G/C content of the guide RNA sequence as a percentage.

    args: 
    grna_dna_sequence (str) : dna complement of gRNA variable region

    returns:
    gc_percentage (float): %GC of the dna complement of the gRNA variable region

    """
    gc_count = 0
    base_count = 0
    for char in grna_dna_sequence:
        if char.isalpha(): #Count the length of the guide RNA
            gc_count += 1
        if char in ['C','c','G','g']:
            base_count += 1 #Count the number of Gs or Cs in the guide RNA
    gc_percentage = float(base_count)/float(gc_count) * 100 #Calculate a percentage of Gs and Cs in the length of the guide RNA
    gc_percentage = round(gc_percentage, 2) #Round to two decimal places
    return gc_percentage

def guide_RNA_notes(grna_dna_sequence, gc_percentage):
    """ 
    Return additional information on the genereated guide RNA sequences.
    args:
    grna_dna_sequence (str) : dna complement of gRNA variable region
    gc_percentage (float): %GC of the dna complement of the gRNA variable region
    
    returns:
    notes (list): list of additional information for the guide RNA sequence 

    """

    # Check if guides have been identified
    if grna_dna_sequence == "No 5' guide could be identified" or grna_dna_sequence == "No 3' guide could be identified":
        return ""
    
    # Initialise the notes  
    notes = []
    for n in range(len(grna_dna_sequence)-3):
        if grna_dna_sequence[n:n + 4] == 'TTTT': #Check for four thymines in a row
            notes.append('PolyT present.')
    if grna_dna_sequence[0] != 'G': #Check the guide RNA starts with a 'G' at the most 5' position
        notes.append('No leading G.')
    if gc_percentage >= 75:
        notes.append('G/C content over 75%.') #Check if the G/C content of the guide is more than or equal to 75%
    return notes

def specific_function_guide_RNA_generation(spec_aa_pos, loci, reference_genome, minimum_distance = 5, maximum_distance = 86, pam = 'NGG', ):
    """ 
    Identifies and produces metrics for guides upstream and downstream from a specific amino acid residue.
    
    Arguments:
        spec_aa_pos - specific amino acid residue selected by position
        loci - DNA loci for protein of interest 
        reference_genome - reference genome used
        minimum_distance - minimum search distance to find guide RNA (in bases pairs from the codon for the specific amino acid)
        minimum_distance - maximum search distance to find guide RNA (in bases pairs from the codon for the specific amino acid)
        pam - pam sequence e.g. NGG for spCas9
        
            
    Returns:
        guides (pd.Dataframe): pandas dataframe containing information about all guide pairs
    """
  
    loci, loci_exon = clean_loci(loci)
    print(f'loci (length:{len(loci)})')
    print(loci)
    print(f'loci exon (length:{len(loci_exon)})')
    print(loci_exon)
    
    # convert the CDS to a protein dictionary
    protein_dict = get_codon_index(loci, loci_exon)
    
    # display information about the selected amino acid
    spec_aa_pos = int(spec_aa_pos)
    maximum_distance = int(maximum_distance)
    aa_selected = protein_dict[spec_aa_pos - 1] #-1 corrects for python 0 indexing
    # Confirm the selection is the correct amino acid
    print(f'The amino acid you have selected is {aa_selected} at aa position {spec_aa_pos} and loci postion {protein_dict[spec_aa_pos - 1]["base_1"][0]} - {protein_dict[spec_aa_pos - 1]["base_3"][0]}')  
    print(f'The PAM you have selected is {pam}')
    
    # generate lists of all the potential guides within a subsection genomic loci and a list of all guides on the reverse complement.
    loci_region_start = protein_dict[spec_aa_pos - 1]["base_1"][0] - maximum_distance - 25
    loci_region_end = protein_dict[spec_aa_pos - 1]["base_3"][0] + maximum_distance + 25
    loci_search_region = loci[loci_region_start:loci_region_end] # Make a smaller region of the loci for loci off target checking + easier visualisation
    print(f'Searching for gRNAs in loci region {loci_region_start}-{loci_region_end}...')
    
    #loci check region for appending onto output dataframe for easy checks and to be used in the 
    loci_check_region = loci[protein_dict[spec_aa_pos - 1]["base_1"][0] - 1000: protein_dict[spec_aa_pos - 1]["base_3"][0] + 1000]

    # Identify PAMs in the inputted genomic loci
    pam_found , guide_positions, gRNA_list, guide_strands = find_guides(loci_search_region, pam)  
    
    # convert the lists into dataframes
    guide_df = pd.DataFrame()
    
    # Add a column for the pam
    guide_df["PAM"] = pam_found
    
    # Add a column of guide RNA sequences to the dataframe
    guide_df["gRNA variable region sequence"] = gRNA_list  
    
    # Add a column of guide RNA cut site positions to the dataframe
    guide_df["Position on + strand"] = guide_positions + loci_region_start
                        
    # Add a column to strand relative to the inputted genomic loci
    guide_df["Strand"] = guide_strands
    
    # Sort the guide RNAs by their position in the genomic loci
    guide_df = guide_df.sort_values(by=["Position on + strand"]).reset_index(drop=True)
    
    # Stops code if guide dataframe is empty
    assert not guide_df.empty, "No guide RNAs identified in sequence"
    
    print(f'{(len(guide_df))} guides found')
    
   #get the distance of the gRNAs from the start and end base of the amino acid
    start_base = protein_dict[spec_aa_pos - 1]["base_1"][0] # -1 corrects for python 0 index
    print(start_base)
    end_base = protein_dict[spec_aa_pos - 1]["base_3"][0] # -1 corrects for python 0 index
    print(end_base)
    guide_df["Distance from start of Amino Acid (bp)"] = guide_df["Position on + strand"] - start_base
    guide_df["Distance from end of Amino Acid (bp)"] = guide_df["Position on + strand"] - end_base - 1 # -1 ensures that distances are not inclusive of the end base
  
  

    # define a dataframe for if there are no guides within the distance range
    noguides = pd.DataFrame(columns=["loci_region","PAM","gRNA variable region sequence","Position on + strand","Strand","Distance from Amino Acid (bp)"])  
       
    # get the upstream guide RNA dataframe
    upstream_guides = []
    upstream_guides = guide_df.query(f"`Distance from start of Amino Acid (bp)` <= {-minimum_distance}").query(f"`Distance from start of Amino Acid (bp)` >= {-maximum_distance}")
    upstream_guides = upstream_guides.rename(columns={"Distance from start of Amino Acid (bp)": "Distance from Amino Acid (bp)"})
    upstream_guides = upstream_guides.drop("Distance from end of Amino Acid (bp)", axis=1)
    print('upstream_guides')
    print(upstream_guides)
    
    if len(upstream_guides) == 0:
        upstream_guides = noguides
    
    # get the downstream guide RNA dataframe
    downstream_guides = []
    downstream_guides = guide_df.query(f"`Distance from end of Amino Acid (bp)` >= {minimum_distance}").query(f"`Distance from end of Amino Acid (bp)` <= {maximum_distance}")
    downstream_guides = downstream_guides.rename(columns={"Distance from end of Amino Acid (bp)": "Distance from Amino Acid (bp)"})
    downstream_guides = downstream_guides.drop("Distance from start of Amino Acid (bp)", axis=1)
    print('downstream_guides')
    print(downstream_guides)
    
    if len(downstream_guides) == 0:
        downstream_guides = noguides
        guides = pd.concat([upstream_guides,downstream_guides])
        return guides

    # join the dataframes together
    guides = pd.concat([upstream_guides,downstream_guides])
    guides = guides.reset_index(drop = True)
    
    # get the G/C content
    guides["G/C Content (%)"] = guides.apply(
        lambda row: (
            get_GC_content(row["gRNA variable region sequence"])
            if row["gRNA variable region sequence"] != ""
            else ""
        ),
        axis=1
    )
    # get the notes
    guides["Notes"] = guides.apply(
        lambda row: (
            guide_RNA_notes(row["gRNA variable region sequence"], row["G/C Content (%)"])
            if row["gRNA variable region sequence"] != ""
            else ""
        ),
        axis=1
    )  
    
    # Perform off target searches
    guides["off target count (loci)"] = guides.apply(
    lambda row: loci.count(row['gRNA variable region sequence'] + row['PAM']) + 
                 loci.count(reverse_complement(row['gRNA variable region sequence'] + row['PAM'])) - 1,  
    axis=1
)
    guides["off target count (genome)"] = guides.apply(
        lambda row: sum(
            record.seq.count((row['gRNA variable region sequence'] + row['PAM']).upper()) +
            record.seq.count(reverse_complement((row['gRNA variable region sequence'] + row['PAM'] ).upper()))
            for record in reference_genome
        ) - 1,  # Adjusting for the self-match if necessary
        axis=1
    )
    
    # Add a column for loci_region
    guides["loci_region"] = str(loci_check_region)
    # Reorganise the guide RNA dataframe
    guides = guides[
        [
            "PAM",
            "Distance from Amino Acid (bp)",
            "gRNA variable region sequence",
            "Position on + strand",
            "Strand",
            "G/C Content (%)",
            "off target count (loci)",
            "off target count (genome)",
            "Notes",
            "loci_region"
        ]
    ]
    guides = guides.reset_index(drop = True)
    
    return guides