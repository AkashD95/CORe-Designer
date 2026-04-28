import sys
import pandas as pd 
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem, QPushButton
from ui.CORe_designer import Ui_MainWindow  # Import the generated UI class
from app_uniprot_id_to_genomic_loci_to_guide_generation import extract_genomic_information_from_uniprot_id, extract_genomic_loci_from_genomic_information,find_guides, get_codon_index, get_GC_content, guide_RNA_notes, specific_function_guide_RNA_generation
from app_selected_guides_to_hdr_template import generate_hdr_template, recodonise_hdr_template, reconstruct_loci_with_recodonised_hdr_template, generate_hdr_library
from app_integration_specific_primers import design_integration_specific_primers
from app_utils import extract_genome_multifast_to_list_of_sequence_records
import time
import logging

# Configure logging
logging.basicConfig(
    filename='CORe_designer.log',
    level=logging.ERROR,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

#Class for pop up table needed for guide RNA selection
class PopUpTableGuides(QDialog):
    def __init__(self, df, parent=None):
        '''
        Pop-up table for selecting guide RNAs.

        Args:
            df (pd.DataFrame): DataFrame containing guide RNA information.
            parent (QWidget, optional): Parent widget for the dialog.
        '''

        super().__init__(parent)
        self.setWindowTitle("Select two gRNAs. One upstream and one downstream. If generating homology directed repair templates we recommend a reverse strand upstream guide and a forward strand downstream guide.")
        self.resize(400, 300)
        
        # Layout for the dialog
        layout = QVBoxLayout(self)
        
        # Create the table
        self.table = QTableWidget(self)
        self.table.setRowCount(df.shape[0])  # Set number of rows
        self.table.setColumnCount(df.shape[1])  # Set number of columns
        self.table.setHorizontalHeaderLabels(df.columns.tolist())  # Set column headers
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.MultiSelection)  # Allow multiple selection
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        
        # Fill the table with data from the DataFrame
        for row_idx, row_data in df.iterrows():
            for col_idx, value in enumerate(row_data):
                self.table.setItem(row_idx, col_idx, QTableWidgetItem(str(value)))
        
        layout.addWidget(self.table)
        
        # Add a select button
        self.select_button = QPushButton("Select")
        self.select_button.clicked.connect(self.get_selected_rows)
        layout.addWidget(self.select_button)
        
        #Initialise variable
        self.selected_guides = None
    
    def get_selected_rows(self):
        selected_items = self.table.selectedItems()
        
        # Organize selected items by rows
        selected_rows = set(item.row() for item in selected_items)
        if len(selected_rows) != 2:
            QMessageBox.warning(self, "Selection Error", "Please select two guide RNAs.")
            return
        
        # Retrieve data for each selected row
        self.selected_guides = [
            [self.table.item(row, col).text() for col in range(self.table.columnCount())]
            for row in selected_rows
        ]
        self.accept()  # Close the dialog

class PopUpTableLoci(QDialog):
    def __init__(self, df, parent=None):
        '''
        Pop-up table for selecting genomic loci

        Args:
            df (pd.DataFrame): DataFrame containing guide RNA information.
            parent (QWidget, optional): Parent widget for the dialog.
        '''

        super().__init__(parent)
        self.setWindowTitle("Select desired transcript ID.")
        self.resize(400, 300)
        
        # Layout for the dialog
        layout = QVBoxLayout(self)
        
        # Create the table
        self.table = QTableWidget(self)
        self.table.setRowCount(df.shape[0])  # Set number of rows
        self.table.setColumnCount(df.shape[1])  # Set number of columns
        self.table.setHorizontalHeaderLabels(df.columns.tolist())  # Set column headers
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.MultiSelection)  # Allow multiple selection
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        
        # Fill the table with data from the DataFrame
        for row_idx, row_data in df.iterrows():
            for col_idx, value in enumerate(row_data):
                self.table.setItem(row_idx, col_idx, QTableWidgetItem(str(value)))
        
        layout.addWidget(self.table)
        
        # Add a select button
        self.select_button = QPushButton("Select")
        self.select_button.clicked.connect(self.get_selected_rows)
        layout.addWidget(self.select_button)
        
        #Initialise variable
        self.selected_genomic_information = None
    
    def get_selected_rows(self):
        selected_items = self.table.selectedItems()
        
        # Organize selected items by rows
        selected_rows = set(item.row() for item in selected_items)
        if len(selected_rows) != 1:
            QMessageBox.warning(self, "Selection Error", "Please select two guide RNAs.")
            return
        
        # Retrieve data for each selected row
        self.selected_genomic_information = [
            [self.table.item(row, col).text() for col in range(self.table.columnCount())]
            for row in selected_rows
        ]
        self.accept()  # Close the dialog

class MainApp(QtWidgets.QMainWindow):
    def __init__(self):
        '''
        Initialize the main application window.
        '''
        super().__init__()
        self.initialize_ui()

    def initialize_ui(self):
        '''
        Set up the user interface, variables and  connect buttons to functions.
        '''

        self.ui = Ui_MainWindow()  # Create an instance of the UI class
        self.ui.setupUi(self)  # Set up the UI on this QMainWindow instance

        #For uniprot locus tab
        # Attributes to store variables
        self.selected_reference_genome = None  
        self.selected_codon_usage_table = None
        self.saved_text_uniprot_id = None
        self.saved_text_amino_acid_number = None
        self.saved_text_custom_locus = None
        self.saved_text_gRNA_distance = None
        self.loci = None
        self.genomic_information_loci_extracted_plus_guides = None
        self.homology_arm_length = None
        self.hdr_template_recodonised = None
        self.upstream_homology_arm_end_position = None
        self.downstream_homology_arm_start_position = None
        self.genomic_information_loci_extracted_plus_guides_plus_hdr = None
        self.loci_plus_recodonised_hdr_template = None

        #Call on text inputs from the QLineEdit boxes
        self.ui.uniprot_id_line_edit.editingFinished.connect(self.save_text_uniprot_id)
        self.ui.amino_acid_position_line_edit.editingFinished.connect(self.save_text_amino_acid_position)
        self.ui.gRNA_distance_line_edit.editingFinished.connect(self.save_text_gRNA_distance)
        self.ui.homology_arm_length_lineEdit.editingFinished.connect(self.save_text_homology_arm_length)
        
        # Connect buttons to functions
        self.ui.upload_reference_genome_button.clicked.connect(self.upload_reference_genome_action)
        self.ui.upload_codon_usage_table_button.clicked.connect(self.upload_codon_usage_table_action)
        self.ui.get_sgRNAs_button.clicked.connect(self.get_gRNAs_action)
        self.ui.generate_hdr_template_button.clicked.connect(self.get_hdr_template_action)
        self.ui.generate_integration_specifc_primer_button.clicked.connect(self.get_integration_specific_primers_action)
        self.ui.save_button.clicked.connect(self.save_action)
        self.ui.reset_button.clicked.connect(self.reset_action)  

        #For custom locus tab
        #Call on text inputs from the QLineEdit boxes
        self.ui.upload_codon_usage_table_button_custom.clicked.connect(self.upload_codon_usage_table_action)
        self.ui.custom_loci_line_edit.editingFinished.connect(self.save_text_custom_locus)
        self.ui.amino_acid_position_line_edit_custom.editingFinished.connect(self.save_text_amino_acid_position_custom)
        self.ui.gRNA_distance_line_edit_custom.editingFinished.connect(self.save_text_gRNA_distance_custom)
        self.ui.homology_arm_length_lineEdit_custom.editingFinished.connect(self.save_text_homology_arm_length_custom)
        
        # Connect buttons to functions
        self.ui.get_sgRNAs_button_custom.clicked.connect(self.get_gRNAs_action_custom)
        self.ui.generate_hdr_template_button_custom.clicked.connect(self.get_hdr_template_action_custom)
        self.ui.generate_integration_specifc_primer_button_custom.clicked.connect(self.get_integration_specific_primers_action_custom)
        self.ui.save_button_custom.clicked.connect(self.save_action_custom)
        self.ui.reset_button_custom.clicked.connect(self.reset_action)

            
    #Functions to save text
    def save_text_uniprot_id(self):
        '''
        Save the UniProt ID from the QLineEdit.
        '''
        print('Saving uniprot id:')
        # Get the text from the QLineEdit
        text = self.ui.uniprot_id_line_edit.text()
        print(text)

        # Check if the text is not empty
        if text.strip():  
            # Save the text to a variable
            self.saved_text_uniprot_id = text
        else:
            QMessageBox.warning(self, "No Input", "Please enter text before leaving the field.")
    
    def save_text_amino_acid_position(self):
        '''
        Save the amino acid position from the QLineEdit.
        '''
        print('Saving amino acid position:')
        # Get the text from the QLineEdit
        text = self.ui.amino_acid_position_line_edit.text()
        print(text)

        # Check if the text is not empty
        if text.strip():  
            self.saved_text_amino_acid_position = text
        else:
            QMessageBox.warning(self, "No Input", "Please enter text before leaving the field.")

    def save_text_amino_acid_position_custom(self):
        '''
        Save the amino acid position from the QLineEdit.
        '''
        print('Saving amino acid position:')
        # Get the text from the QLineEdit
        text = self.ui.amino_acid_position_line_edit_custom.text()
        print(text)


        if text.strip():  # Check if the text is not empty
            self.saved_text_amino_acid_position = text
        else:
            QMessageBox.warning(self, "No Input", "Please enter text before leaving the field.")

    def save_text_custom_locus(self):
        '''
        Save the custom locus from the QLineEdit.
        '''
        print('Saving custom locus:')
        # Get the text from the QLineEdit
        text = self.ui.custom_loci_line_edit.text()
        print(text)

        # Check if the text is not empty
        if text.strip():  
            self.saved_text_custom_locus = text
        else:
            QMessageBox.warning(self, "No Input", "Please enter text before leaving the field.")

    def save_text_gRNA_distance(self):
        '''
        Save the gRNA distance from the QLineEdit.
        '''
        print('Saving gRNA distance:')
        # Get the text from the QLineEdit
        text = self.ui.gRNA_distance_line_edit.text()
        print(text)

        # Check if the text is not empty
        if text.strip():  
            self.saved_gRNA_distance = text
        else:
            QMessageBox.warning(self, "No Input", "Please enter text before leaving the field.")

    def save_text_gRNA_distance_custom(self):
        '''
        Save the gRNA distance from the QLineEdit.
        '''
        print('Saving gRNA distance:')
        # Get the text from the QLineEdit
        text = self.ui.gRNA_distance_line_edit_custom.text()
        print(text)

        # Check if the text is not empty
        if text.strip():  
            self.saved_gRNA_distance = text
        else:
            QMessageBox.warning(self, "No Input", "Please enter text before leaving the field.")
    
    def save_text_homology_arm_length(self):
        '''
        Save the homology arm length from the QLineEdit.
        '''
        # Get the text from the QLineEdit
        print('Saving homology arm length:')
        text = self.ui.homology_arm_length_lineEdit.text()
        print(text)

        # Check if the text is not empty
        if text.strip():  
            self.homology_arm_length = text
        else:
            QMessageBox.warning(self, "No Input", "Please enter text before leaving the field.")

    def save_text_homology_arm_length_custom(self):
        '''
        Save the homology arm length from the QLineEdit.
        '''
        # Get the text from the QLineEdit
        print('Saving homology arm length:')
        text = self.ui.homology_arm_length_lineEdit_custom.text()
        print(text)

        # Check if the text is not empty
        if text.strip():    
            self.homology_arm_length = text
        else:
            QMessageBox.warning(self, "No Input", "Please enter text before leaving the field.")
    
    #Functions to populate results tables 
    def populate_results_table(self, dataframe):
        """
        Populate QTableWidget with the contents of a Pandas DataFrame.
        """
        self.ui.results_table.setRowCount(len(dataframe))
        self.ui.results_table.setColumnCount(len(dataframe.columns))

        # Set the column headers
        self.ui.results_table.setHorizontalHeaderLabels(dataframe.columns)

        # Fill the table
        for row in range(len(dataframe)):
            for column in range(len(dataframe.columns)):
                value = str(dataframe.iloc[row, column])
                item = QTableWidgetItem(value)
                self.ui.results_table.setItem(row, column, item)
    
    def populate_results_table_custom(self, dataframe):
        """
        Populate QTableWidget with the contents of a Pandas DataFrame.
        """
        self.ui.results_table_custom.setRowCount(len(dataframe))
        self.ui.results_table_custom.setColumnCount(len(dataframe.columns))

        # Set the column headers
        self.ui.results_table_custom.setHorizontalHeaderLabels(dataframe.columns)

        # Fill the table
        for row in range(len(dataframe)):
            for column in range(len(dataframe.columns)):
                value = str(dataframe.iloc[row, column])
                item = QTableWidgetItem(value)
                self.ui.results_table_custom.setItem(row, column, item)
    
     
    def populate_hdr_library_table(self, dataframe):
        """
        Populate QTableWidget with the contents of a Pandas DataFrame.
        """
        self.ui.hdr_library_table.setRowCount(len(dataframe))
        self.ui.hdr_library_table.setColumnCount(len(dataframe.columns))

        # Set the column headers
        self.ui.hdr_library_table.setHorizontalHeaderLabels(dataframe.columns)

        # Fill the table
        for row in range(len(dataframe)):
            for column in range(len(dataframe.columns)):
                value = str(dataframe.iloc[row, column])
                item = QTableWidgetItem(value)
                self.ui.hdr_library_table.setItem(row, column, item)
    
    def populate_hdr_library_table_custom(self, dataframe):
        """
        Populate QTableWidget with the contents of a Pandas DataFrame.
        """
        self.ui.hdr_library_table_custom.setRowCount(len(dataframe))
        self.ui.hdr_library_table_custom.setColumnCount(len(dataframe.columns))
        # Set the column headers
        self.ui.hdr_library_table_custom.setHorizontalHeaderLabels(dataframe.columns)

        # Fill the table
        for row in range(len(dataframe)):
            for column in range(len(dataframe.columns)):
                value = str(dataframe.iloc[row, column])
                item = QTableWidgetItem(value)
                self.ui.hdr_library_table_custom.setItem(row, column, item)

    #Functions to run after buttons are pressed
    def upload_reference_genome_action(self):
        '''
        Function to upload a reference genome file.
        '''
        # Open a file dialog to select a file
        file_name, _ = QFileDialog.getOpenFileName(None, "Select File", "", "All Files (*.*);;Text Files (*.txt);;FASTA Files (*.fasta)")

        if file_name:  # If a file was selected
            QMessageBox.information(None, "File Selected", f"You selected: {file_name}")
            self.selected_reference_genome = file_name  # Store the file path
            print(f"Selected file: {file_name}")
        else:
            QMessageBox.warning(None, "No File", "No file was selected.")

    def upload_codon_usage_table_action(self):
        '''
        Function to upload a codon usage table file.
        '''
        # Open a file dialog to select a file
        file_name, _ = QFileDialog.getOpenFileName(None, "Select File", "", "All Files (*.*);;Text Files (*.txt);;FASTA Files (*.fasta)")

        if file_name:  # If a file was selected
            QMessageBox.information(None, "File Selected", f"You selected: {file_name}")
            self.selected_codon_usage_table = file_name  # Store the file path
            print(f"Selected file: {file_name}")
        else:
            QMessageBox.warning(None, "No File", "No file was selected.")
    
    def get_gRNAs_action(self):
        '''
        Function for workflow if get_'sgRNAs' button is pressed
        '''

        # First get the uniprot ID and process the file against the uploaded reference genome
        # Define the process when get_'sgRNAs' button is clicked
        start = time.perf_counter()
        # Error handling if file has not been selected 
        if not self.selected_reference_genome:  # Check if a file has been selected
            QMessageBox.information(None, "No reference genome file found ", "Please upload a file first.")
            return
        
        # Convert the file into a list of genomic records 
        QMessageBox.information(None, "Starting the guide generation process", "Click OK to continue (this may take a couple of minutes), a table with generated guides for selection will pop up after this process is complete.")
        try:
            genome_sequence_records = extract_genome_multifast_to_list_of_sequence_records(self.selected_reference_genome)
            print('genome sequence records:')
            print(genome_sequence_records)
            uniprot_id = self.saved_text_uniprot_id
            print('uniprot id:')
            print(uniprot_id)
            amino_acid_position = self.saved_text_amino_acid_position
            print('amino_acid_position')
            print(amino_acid_position)
            genomic_information = extract_genomic_information_from_uniprot_id(uniprot_id)
            print(genomic_information)
            QMessageBox.information(None, "Uniprot record found", "Select desired transcript ID from the pop up table. Press OK to continue.")
            
            genomic_information_pop_up = PopUpTableLoci(genomic_information, self)
            if genomic_information_pop_up.exec_() == QDialog.Accepted:
                selected_genomic_information = genomic_information_pop_up.selected_genomic_information
                selected_genomic_information_df = pd.DataFrame(selected_genomic_information, columns=genomic_information.columns)

            genomic_information = selected_genomic_information_df.reset_index(drop=True)

            genomic_information_loci_extracted = extract_genomic_loci_from_genomic_information(genomic_information,genome_sequence_records)
            print(genomic_information_loci_extracted)
            
            #Add amino acid position to this dataframe
            genomic_information_loci_extracted['Amino_acid_position'] = amino_acid_position

            #Extract the loci to run through the guide RNA screen
            loci = genomic_information_loci_extracted['genomic_loci'].iloc[0]
            self.loci = loci #create class variable to be used in other functions
            print('loci')
            print(loci)
            grna_distance = self.saved_gRNA_distance
            guides = specific_function_guide_RNA_generation(amino_acid_position, loci, genome_sequence_records, minimum_distance=5, maximum_distance = grna_distance)
            
            print('All guides')
            print(guides)

            end = time.perf_counter()
            time_taken = end - start 
            print(time_taken)
        
        except Exception as e:  
            print("An error occurred during guide RNA generation:", str(e))
            logging.error("An error occurred during guide RNA generation", exc_info=True)
            QMessageBox.critical(None, "Error", f"An error occurred and information has been logged.")
            return

        QMessageBox.information(None, "Guides are generated", "Select two guides from the pop up table. Press OK to continue.")
        
        # Display the pop-up table at the end of this function
        guide_pop_up = PopUpTableGuides(guides, self)
        if guide_pop_up.exec_() == QDialog.Accepted:
            selected_guides = guide_pop_up.selected_guides
            selected_guide_df = pd.DataFrame(selected_guides, columns=guides.columns)
            self.selected_guide_df = selected_guide_df.sort_values(by='Distance from Amino Acid (bp)').reset_index(drop=True) 

             
        selected_guides = self.selected_guide_df
        print('Selected guides are:')
        print(selected_guides)
        
        #add guide information to the genomic_information_loci_extracted DataFrame
        upstream_guide = selected_guides.loc[[0]].reset_index(drop = True)
        print(upstream_guide)
        upstream_guide.columns = ['Upstream guide RNA ' + col for col in upstream_guide.columns]
        print(upstream_guide)
        genomic_information_loci_extracted_plus_guides = pd.concat([genomic_information_loci_extracted, upstream_guide], axis = 1)
        downstream_guide = selected_guides.loc[[1]].reset_index(drop = True)
        print(downstream_guide)
        downstream_guide.columns = ['Downstream guide RNA ' + col for col in downstream_guide]
        print(downstream_guide)
        genomic_information_loci_extracted_plus_guides = pd.concat([genomic_information_loci_extracted_plus_guides, downstream_guide], axis = 1)
        print(genomic_information_loci_extracted_plus_guides)
        self.genomic_information_loci_extracted_plus_guides = genomic_information_loci_extracted_plus_guides
        # Update the results table
        self.populate_results_table(genomic_information_loci_extracted_plus_guides)
    
    def get_gRNAs_action_custom(self):
        # First get the uniprot ID and process the file against the uploaded reference genome
        # Define the process when get_'sgRNAs' button is clicked
        start = time.perf_counter() 
        # Convert the file into a list of genomic records 
        QMessageBox.information(None, "Starting the guide generation process", "Click OK to continue (this may take a couple of minutes), a table with generated guides for selection will pop up after this process is complete.")
        try:
            amino_acid_position = self.saved_text_amino_acid_position
            #Extract the loci to run through the guide RNA screen
            loci = self.saved_text_custom_locus
            self.loci = loci #create class variable to be used in other functions
            grna_distance = self.saved_gRNA_distance
            genome_sequence_records = None 
            guides = specific_function_guide_RNA_generation(amino_acid_position, loci, genome_sequence_records, minimum_distance=5, maximum_distance = grna_distance)
            print('All guides')
            print(guides)
            end = time.perf_counter()
            time_taken = end - start 
            print(time_taken)

        except Exception as e:  
            print("An error occurred during guide RNA generation:", str(e))
            logging.error("An error occurred during guide RNA generation", exc_info=True)
            QMessageBox.critical(None, "Error", f"An error occurred and information has been logged.")
            return

        QMessageBox.information(None, "Guides are generated", "Select two guides from the pop up table. Press OK to continue.")
        # Display the pop-up table at the end of this function
        guide_pop_up = PopUpTable(guides, self)
        if guide_pop_up.exec_() == QDialog.Accepted:
            selected_guides = guide_pop_up.selected_guides
            selected_guide_df = pd.DataFrame(selected_guides, columns=guides.columns)
            self.selected_guide_df = selected_guide_df.sort_values(by='Distance from Amino Acid (bp)').reset_index(drop=True) 
        
        selected_guides = self.selected_guide_df
        print('Selected guides are:')
        print(selected_guides)
        
        #add guide information to the genomic_information_loci_extracted DataFrame
        upstream_guide = selected_guides.loc[[0]].reset_index(drop = True)
        print(upstream_guide)
        upstream_guide.columns = ['Upstream guide RNA ' + col for col in upstream_guide.columns]
        print(upstream_guide)
        downstream_guide = selected_guides.loc[[1]].reset_index(drop = True)
        print(downstream_guide)
        downstream_guide.columns = ['Downstream guide RNA ' + col for col in downstream_guide]
        print(downstream_guide)
        genomic_information_loci_extracted_plus_guides = pd.concat([upstream_guide, downstream_guide], axis = 1)
        print(genomic_information_loci_extracted_plus_guides)
        self.genomic_information_loci_extracted_plus_guides = genomic_information_loci_extracted_plus_guides
        
        #Update the results table
        self.populate_results_table_custom(genomic_information_loci_extracted_plus_guides)

    def get_hdr_template_action(self):
        '''function for downstream workflow if generate HDR template is pressed''' 
        start = time.perf_counter()
        # message to start function
        QMessageBox.information(None, "Generating HDR template", "Click OK to begin HDR template generation.")

        #Extract factors from get_gRNAs_action
        selected_guides = self.selected_guide_df
        print('Selected guides are:')
        print(selected_guides)

        #add guide information to the genomic_information_loci_extracted DataFrame
        upstream_guide = selected_guides.loc[[0]]
        upstream_guide.columns = ['Upstream guide RNA ' + col for col in upstream_guide.columns]
        downstream_guide = selected_guides.loc[[1]]
        downstream_guide.columns = ['Downstream guide RNA ' + col for col in downstream_guide]
        loci = str(self.loci)
        upstream_guide_position = int(upstream_guide['Upstream guide RNA Position on + strand'])
        downstream_guide_position = int(downstream_guide['Downstream guide RNA Position on + strand'])

        #Get homology arm length
        homology_arm_length = int(self.homology_arm_length)
        print('Inputted homology arm length:')
        print(homology_arm_length)

        hdr_template, upstream_homology_arm, upstream_homology_arm_start_position, upstream_homology_arm_end_position, downstream_homology_arm, downstream_homology_arm_start_position, downstream_homology_arm_end_postion = generate_hdr_template(loci, upstream_guide_position, downstream_guide_position, homology_overlap = homology_arm_length)
        print(hdr_template)

        upstream_guide_distance = int(upstream_guide['Upstream guide RNA Distance from Amino Acid (bp)'])
        downstream_guide_distance = int(downstream_guide['Downstream guide RNA Distance from Amino Acid (bp)'])
        
        try:
            #Get codon usage table and process into dictionary
            codon_usage_table_file = self.selected_codon_usage_table
            codon_usage_table = pd.read_excel(codon_usage_table_file)
            codon_usage_table_dict = dict(zip(codon_usage_table['Codon'], codon_usage_table['Frequency']))
            hdr_template_recodonised, spec_aa_codon_start, recodonisation_check = recodonise_hdr_template(hdr_template, upstream_guide_distance, downstream_guide_distance, codon_usage_table = codon_usage_table_dict)
            print('Recodonised HDR template:')
            print(hdr_template_recodonised)
        
            #Save to class space so that variables can be used in other functions 
            self.upstream_guide_distance = upstream_guide_distance
            self.downstream_guide_distance = downstream_guide_distance
            self.upstream_homology_arm_end_position = hdr_template_recodonised
            self.upstream_homology_arm_end_position = upstream_homology_arm_end_position
            self.downstream_homology_arm_start_position = downstream_homology_arm_start_position

            hdr_list = [hdr_template, hdr_template_recodonised, upstream_homology_arm, downstream_homology_arm]
            print(hdr_list)
            hdr_df = pd.DataFrame(data = [hdr_list], columns= ['HDR template', 'recodonised HDR template', 'Upstream homology arm', 'Downstream homology arm' ])
            print(hdr_df)

            genomic_information_loci_extracted_plus_guides = self.genomic_information_loci_extracted_plus_guides
            genomic_information_loci_extracted_plus_guides_plus_hdr = pd.concat([genomic_information_loci_extracted_plus_guides, hdr_df], axis = 1)
            self.genomic_information_loci_extracted_plus_guides_plus_hdr = genomic_information_loci_extracted_plus_guides_plus_hdr
            print(genomic_information_loci_extracted_plus_guides_plus_hdr)

            #Reconstruct the loci with the recodonised HDR template   
            loci_plus_recodonised_hdr_template, section_loci_plus_recodonised_hdr_template = reconstruct_loci_with_recodonised_hdr_template(loci, hdr_template_recodonised, upstream_homology_arm_end_position, downstream_homology_arm_start_position)
            self.loci_plus_recodonised_hdr_template = loci_plus_recodonised_hdr_template

            #Generate hdr_library for target of interest
            #Take the top frequency codon for each amino acid subsitution
            codon_usage_table_sorted = codon_usage_table.sort_values(['Amino_Acid', 'Frequency'], ascending=[True, False])
            
            # Keep the first row for each AminoAcid (highest frequency)
            top_codon = codon_usage_table_sorted.groupby('Amino_Acid', as_index=False).first()
            top_codon_list = list(zip(top_codon['Amino_Acid'], top_codon['Codon']))
            hdr_library = generate_hdr_library(hdr_template_recodonised, spec_aa_codon_start, upstream_homology_arm, downstream_homology_arm, codon_list = top_codon_list)
            print('Generated HDR library:')
            end = time.perf_counter()
            time_taken = end - start 
            print(time_taken)
            print(hdr_library)
            hdr_library = pd.DataFrame(hdr_library, columns=["HDR_library"])
            print('HDR templates generated')
        except Exception as e:
            print("An error occurred during HDR template recodonisation:", str(e))
            logging.error("An error occurred during HDR template recodonisation", exc_info=True)
            QMessageBox.critical(None, "Error", f"An error occurred and information has been logged.")
            return

        # Update the results table
        QMessageBox.information(None, "HDR template generated", "Click OK to add the information to the results table.")
        self.populate_results_table(genomic_information_loci_extracted_plus_guides_plus_hdr)
        self.populate_hdr_library_table(hdr_library)
    
    def get_hdr_template_action_custom(self):
        '''function for downstream workflow if generate HDR template is pressed''' 
        start = time.perf_counter()
        #message to start function
        QMessageBox.information(None, "Generating HDR template", "Click OK to begin HDR template generation.")

        #Extract factors from get_gRNAs_action
        selected_guides = self.selected_guide_df
        print('Selected guides are:')
        print(selected_guides)

        #add guide information to the genomic_information_loci_extracted DataFrame
        upstream_guide = selected_guides.loc[[0]]
        upstream_guide.columns = ['Upstream guide RNA ' + col for col in upstream_guide.columns]
        downstream_guide = selected_guides.loc[[1]]
        downstream_guide.columns = ['Downstream guide RNA ' + col for col in downstream_guide]
        loci = str(self.loci)
        print(loci)
        upstream_guide_position = int(upstream_guide['Upstream guide RNA Position on + strand'])
        print(upstream_guide_position)
        downstream_guide_position = int(downstream_guide['Downstream guide RNA Position on + strand'])
        print(downstream_guide_position)
        #Get homology arm length
        homology_arm_length = int(self.homology_arm_length)
        print('Inputted homology arm length:')
        print(homology_arm_length)
        hdr_template, upstream_homology_arm, upstream_homology_arm_start_position, upstream_homology_arm_end_position, downstream_homology_arm, downstream_homology_arm_start_position, downstream_homology_arm_end_postion = generate_hdr_template(loci, upstream_guide_position, downstream_guide_position, homology_overlap = homology_arm_length)
        print(hdr_template)

        upstream_guide_distance = int(upstream_guide['Upstream guide RNA Distance from Amino Acid (bp)'])
        downstream_guide_distance = int(downstream_guide['Downstream guide RNA Distance from Amino Acid (bp)'])
        
        try:
            #Get codon usage table and process into dictionary
            codon_usage_table_file = self.selected_codon_usage_table
            codon_usage_table = pd.read_excel(codon_usage_table_file)
            codon_usage_table_dict = dict(zip(codon_usage_table['Codon'], codon_usage_table['Frequency']))
            
            
            hdr_template_recodonised, spec_aa_codon_start, recodonisation_check = recodonise_hdr_template(hdr_template, upstream_guide_distance, downstream_guide_distance, codon_usage_table=codon_usage_table_dict)
            print(hdr_template_recodonised)
            
            #Save to class space so that variables can be used in other functions 
            self.upstream_guide_distance = upstream_guide_distance
            self.downstream_guide_distance = downstream_guide_distance
            self.upstream_homology_arm_end_position = hdr_template_recodonised
            self.upstream_homology_arm_end_position = upstream_homology_arm_end_position
            self.downstream_homology_arm_start_position = downstream_homology_arm_start_position

            hdr_list = [hdr_template, hdr_template_recodonised, upstream_homology_arm, downstream_homology_arm]
            print(hdr_list)
            hdr_df = pd.DataFrame(data = [hdr_list], columns= ['HDR template', 'recodonised HDR template', 'Upstream homology arm', 'Downstream homology arm' ])
            print(hdr_df)

            genomic_information_loci_extracted_plus_guides = self.genomic_information_loci_extracted_plus_guides
            genomic_information_loci_extracted_plus_guides_plus_hdr = pd.concat([genomic_information_loci_extracted_plus_guides, hdr_df], axis = 1)
            self.genomic_information_loci_extracted_plus_guides_plus_hdr = genomic_information_loci_extracted_plus_guides_plus_hdr
            #print(genomic_information_loci_extracted_plus_guides_plus_hdr)

            #Reconstruct the loci with the recodonised HDR template   
            loci_plus_recodonised_hdr_template, section_loci_plus_recodonised_hdr_template = reconstruct_loci_with_recodonised_hdr_template(loci, hdr_template_recodonised, upstream_homology_arm_end_position, downstream_homology_arm_start_position)
            self.loci_plus_recodonised_hdr_template = loci_plus_recodonised_hdr_template

            #Generate hdr_library
            #Get codon usage table and process into a list of tuples of most popular codons for each amino acid
            codon_usage_table_sorted = codon_usage_table.sort_values(['Amino_Acid', 'Frequency'], ascending=[True, False])
            # Keep the first row for each AminoAcid (highest frequency)
            top_codon = codon_usage_table_sorted.groupby('Amino_Acid', as_index=False).first()
            top_codon_list = list(zip(top_codon['Amino_Acid'], top_codon['Codon']))
            hdr_library = generate_hdr_library(hdr_template_recodonised, spec_aa_codon_start, upstream_homology_arm, downstream_homology_arm, codon_list = top_codon_list)
            print('Generated HDR library')
            end = time.perf_counter()
            time_taken = end - start 
            print(time_taken)
            hdr_library = pd.DataFrame(hdr_library, columns=["HDR_library"])
            print('HDR templates generated')

        except Exception as e:
            print("An error occurred during HDR template recodonisation:", str(e))
            logging.error("An error occurred during HDR template recodonisation", exc_info=True)
            QMessageBox.critical(None, "Error", f"An error occurred and information has been logged.")
            return

       #Update the results table
        QMessageBox.information(None, "HDR template generated", "Click OK to add the information to the results table.")
        self.populate_results_table_custom(genomic_information_loci_extracted_plus_guides_plus_hdr)
        self.populate_hdr_library_table_custom(hdr_library)
       
    
    def get_integration_specific_primers_action(self):
        '''function for downstream workflow if generate integration specific primers button is pressed'''
        #message to start function
        start = time.perf_counter()
        QMessageBox.information(None, "Generating integration specific primers", "Click OK to begin primer generation.")
        
        #Load up variables from the class space
        loci = str(self.loci)
        hdr_template_recodonised = str(self.hdr_template_recodonised)
        upstream_homology_arm_end_position = int(self.upstream_homology_arm_end_position)
        upstream_guide_distance = int(self.upstream_guide_distance)
        downstream_homology_arm_start_position = int(self.downstream_homology_arm_start_position)
        downstream_guide_distance = int(self.downstream_guide_distance)
        loci_plus_recodonised_hdr_template = str(self.loci_plus_recodonised_hdr_template)

        #Run primer 3 function
        try:
            primer_list = design_integration_specific_primers(hdr_template_recodonised, loci_plus_recodonised_hdr_template, upstream_homology_arm_end_position, upstream_guide_distance, downstream_homology_arm_start_position, downstream_guide_distance)
            print(primer_list)
            primer_df = pd.DataFrame(data = [primer_list], columns= ['left_primer_sequence', 'left_primer_tm','right_primer_sequence','right_primer_tm'])
            end = time.perf_counter()
            time_taken = end - start 
            print(time_taken)
            QMessageBox.information(None, "Primers generated", "Click OK to add the information to the results table.")
            genomic_information_loci_extracted_plus_guides_plus_hdr = self.genomic_information_loci_extracted_plus_guides_plus_hdr
            genomic_information_loci_extracted_plus_guides_plus_hdr_plus_primers = pd.concat([genomic_information_loci_extracted_plus_guides_plus_hdr, primer_df], axis = 1)
        
        except Exception as e:
            print("An error occurred during primer design:", str(e))
            logging.error("An error occurred during primer design", exc_info=True)
            QMessageBox.critical(None, "Error", f"An error occurred and information has been logged.")
            return
        
        self.populate_results_table(genomic_information_loci_extracted_plus_guides_plus_hdr_plus_primers)
    
    def get_integration_specific_primers_action_custom(self):
        '''function for downstream workflow if generate integration specific primers button is pressed'''
        #message to start function
        start = time.perf_counter()
        QMessageBox.information(None, "Generating integration specific primers", "Click OK to begin primer generation.")
        
        #Load up variables from the class space
        loci = str(self.loci)
        hdr_template_recodonised = str(self.hdr_template_recodonised)
        print('HDR template recodonised')
        print(hdr_template_recodonised)

        upstream_homology_arm_end_position = int(self.upstream_homology_arm_end_position)
        print('upstream homology arm end position')
        print(upstream_homology_arm_end_position)

        upstream_guide_distance = int(self.upstream_guide_distance)
        print('upstream guide distance')
        print(upstream_guide_distance)

        downstream_homology_arm_start_position = int(self.downstream_homology_arm_start_position)
        print('downstream homology arm start position')
        print(downstream_homology_arm_start_position)

        downstream_guide_distance = int(self.downstream_guide_distance)
        print('downstream guide distance')
        print(upstream_guide_distance)

        loci_plus_recodonised_hdr_template = str(self.loci_plus_recodonised_hdr_template)
        print('locus plus recodonised hdr template')
        print(loci_plus_recodonised_hdr_template)
        
        #Run primer 3 function
        try:
            primer_list = design_integration_specific_primers(hdr_template_recodonised, loci_plus_recodonised_hdr_template, upstream_homology_arm_end_position, upstream_guide_distance, downstream_homology_arm_start_position, downstream_guide_distance)
            #print(primer_list)
            primer_df = pd.DataFrame(data = [primer_list], columns= ['left_primer_sequence', 'left_primer_tm','right_primer_sequence','right_primer_tm'])
            end = time.perf_counter()
            time_taken = end - start 
            print(time_taken)
            QMessageBox.information(None, "Primers generated", "Click OK to add the information to the results table.")
            genomic_information_loci_extracted_plus_guides_plus_hdr = self.genomic_information_loci_extracted_plus_guides_plus_hdr
            genomic_information_loci_extracted_plus_guides_plus_hdr_plus_primers = pd.concat([genomic_information_loci_extracted_plus_guides_plus_hdr, primer_df], axis = 1)
        except Exception as e:
            print("An error occurred during primer design:", str(e))
            logging.error("An error occurred during primer design", exc_info=True)
            QMessageBox.critical(None, "Error", f"An error occurred and information has been logged.")
            return
        self.populate_results_table_custom(genomic_information_loci_extracted_plus_guides_plus_hdr_plus_primers)
        

    def save_action(self):
        """Save the contents of the results_table to a CSV file."""
        QMessageBox.information(None, "Saving main results table", "Click OK to select save location.")
        # Extract the data from the results_table
        data = []
        rows = self.ui.results_table.rowCount()
        columns = self.ui.results_table.columnCount()

        # Collect data from the table
        for row in range(rows):
            row_data = []
            for column in range(columns):
                item = self.ui.results_table.item(row, column)
                row_data.append(item.text() if item else "")
            data.append(row_data)

        # Create a DataFrame from the table data
        column_headers = [self.ui.results_table.horizontalHeaderItem(i).text() for i in range(columns)]
        df = pd.DataFrame(data, columns=column_headers)

        # Open a file dialog to choose the file location
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getSaveFileName(
            self,
            "Save Excel File",
            "",
            "Excel Files (*.xlsx);;All Files (*)",
            options=options
        )

        if file_name:
            # Ensure .xlsx extension
            if not file_name.lower().endswith(".xlsx"):
                file_name += ".xlsx"

            try:
                # Save the DataFrame to an Excel file
                df.to_excel(file_name, index=False, engine="openpyxl")
                print(f"Data saved to {file_name}")
            except Exception as e:
                print(f"Error saving Excel file: {e}")
        
        QMessageBox.information(None, "Saving hdr library table", "Click OK to select save location.")

        # Extract the data from the results_table
        data = []
        rows = self.ui.hdr_library_table.rowCount()
        columns = self.ui.hdr_library_table.columnCount()

        # Collect data from the table
        for row in range(rows):
            row_data = []
            for column in range(columns):
                item = self.ui.hdr_library_table.item(row, column)
                row_data.append(item.text() if item else "")
            data.append(row_data)

        # Create a DataFrame from the table data
        column_headers = [self.ui.hdr_library_table.horizontalHeaderItem(i).text() for i in range(columns)]
        df = pd.DataFrame(data, columns=column_headers)

        # Open a file dialog to choose the file location
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getSaveFileName(
            self,
            "Save Excel File",
            "",
            "Excel Files (*.xlsx);;All Files (*)",
            options=options
        )

        if file_name:
            # Ensure .xlsx extension
            if not file_name.lower().endswith(".xlsx"):
                file_name += ".xlsx"

            try:
                # Save the DataFrame to an Excel file
                df.to_excel(file_name, index=False, engine="openpyxl")
                print(f"Data saved to {file_name}")
            except Exception as e:
                print(f"Error saving Excel file: {e}")
    
    def save_action_custom(self):
        """Save the contents of the results_table to a CSV file."""

        QMessageBox.information(None, "Saving main results table", "Click OK to select save location.")

        # Extract the data from the results_table
        data = []
        rows = self.ui.results_table_custom.rowCount()
        columns = self.ui.results_table_custom.columnCount()

        # Collect data from the table
        for row in range(rows):
            row_data = []
            for column in range(columns):
                item = self.ui.results_table_custom.item(row, column)
                row_data.append(item.text() if item else "")
            data.append(row_data)

        # Create a DataFrame from the table data
        column_headers = [self.ui.results_table_custom.horizontalHeaderItem(i).text() for i in range(columns)]
        df = pd.DataFrame(data, columns=column_headers)

        # Open a file dialog to choose the file location
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getSaveFileName(
            self,
            "Save Excel File",
            "",
            "Excel Files (*.xlsx);;All Files (*)",
            options=options
        )

        if file_name:
            # Ensure .xlsx extension
            if not file_name.lower().endswith(".xlsx"):
                file_name += ".xlsx"

            try:
                # Save the DataFrame to an Excel file
                df.to_excel(file_name, index=False, engine="openpyxl")
                print(f"Data saved to {file_name}")
            except Exception as e:
                print(f"Error saving Excel file: {e}")
        
        """Save the contents of the results_table to a CSV file."""
        QMessageBox.information(None, "Saving hdr library table", "Click OK to select save location.")

        # Extract the data from the results_table
        data = []
        rows = self.ui.hdr_library_table_custom.rowCount()
        columns = self.ui.hdr_library_table_custom.columnCount()

        # Collect data from the table
        for row in range(rows):
            row_data = []
            for column in range(columns):
                item = self.ui.hdr_library_table_custom.item(row, column)
                row_data.append(item.text() if item else "")
            data.append(row_data)

        # Create a DataFrame from the table data
        column_headers = [self.ui.hdr_library_table_custom.horizontalHeaderItem(i).text() for i in range(columns)]
        df = pd.DataFrame(data, columns=column_headers)

        # Open a file dialog to choose the file location
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getSaveFileName(
            self,
            "Save Excel File",
            "",
            "Excel Files (*.xlsx);;All Files (*)",
            options=options
        )

        if file_name:
            # Ensure .xlsx extension
            if not file_name.lower().endswith(".xlsx"):
                file_name += ".xlsx"

            try:
                # Save the DataFrame to an Excel file
                df.to_excel(file_name, index=False, engine="openpyxl")
                print(f"Data saved to {file_name}")
            except Exception as e:
                print(f"Error saving Excel file: {e}")

    def reset_action(self):
        """Reset the entire application to its initial state."""
        # Remove the current UI
        self.centralWidget().deleteLater()
        # Reinitialize everything
        self.initialize_ui()


if __name__ == "__main__":
    print('Running app.py')
    app = QtWidgets.QApplication(sys.argv)
    main_window = MainApp()
    main_window.show()
    sys.exit(app.exec_())