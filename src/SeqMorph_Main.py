# SeqMorph Main.py (Updated)
import re
import os
import logging
import json
from typing import Dict
from input_module import InputHandler, SequenceFetcher
from mutation_module import MutationEngine, CODON_TABLE
from analysis_module import SequenceAnalysisReport
from sequence_structure import SequenceStructure
from gene_library import SequenceValidation

# Configure logging
logging.basicConfig(
    filename='seqmorph.log',
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def main():
    logging.info("Program started")
    print("   Welcome to SeqMorph!  \n==========================")
    print(r"""
   _____            __  __                  _       __        ___  
  / ____|          |  \/  |                | |     /_ |      / _ \ 
 | (___   ___  __ _| \  / | ___  _ __ _ __ | |__    | |     | | | |
  \___ \ / _ \/ _` | |\/| |/ _ \| '__| '_ \| '_ \   | |     | | | |
  ____) |  __/ (_| | |  | | (_) | |  | |_) | | | |  | |  _  | |_| |
 |_____/ \___|\__, |_|  |_|\___/|_|  | .__/|_| |_|  |_| (_)  \___/ 
                 | |                 | |       - BY PRANAVA UPPARLAPALLI                     
                 |_|                 |_|                          
    """)

    handler = InputHandler()
    fetcher = SequenceFetcher()
    sequence_data = {}
    seq_structure = SequenceStructure()
    

    try:
        # Simplified input handling
        accession_id = input("Enter accession ID (NCBI/UniProt): ").strip()
        logging.info(f"Accession ID provided: {accession_id}")
        
        # Fetch and validate sequence
        raw_sequence = fetcher.fetch_sequence(accession_id)
        seq_type = SequenceValidation.detect_seq(raw_sequence)
        
        if not SequenceValidation.validate(raw_sequence, seq_type):
            raise ValueError(f"Invalid {seq_type} sequence")

        # Add to sequence structure
        seq_structure.add_sequence(
            accession_id=accession_id,
            sequence=raw_sequence,
            seq_type=seq_type,
            source="NCBI" if is_ncbi_accession(accession_id) else "UniProt"
        )

        # Initialize mutation engine
        mutation_engine = MutationEngine(seq_structure)
        
        # Get mutation parameters
        def get_valid_mutation_type():
            while True:
                mutation_type = input("Mutation type (point/insertion/deletion): ").lower()
                if mutation_type in {"point", "insertion", "deletion"}:
                    return mutation_type
                else:
                    print("Invalid mutation type. Must be point, insertion, or deletion.")
        mutation_type = get_valid_mutation_type()
        mutation_rate = float(input("Mutation rate (%): "))
        start = int(input("Start position: "))
        end = int(input("End position: "))

        # Apply mutations
        original = seq_structure.get_sequence(accession_id)
        mutated = mutation_engine.apply_mutations(
            accession_id, mutation_rate, mutation_type, start, end
        )
        # Save outputs
        output_dir = "output"
        os.makedirs(output_dir, exist_ok=True)
        # Save analysis report
        report = SequenceAnalysisReport(original, mutated, seq_type).compare_sequences()

        # Save mutated sequence
        with open(f"{output_dir}/{accession_id}_mutated.fasta", "w") as f:
            f.write(f">Mutated_{accession_id}\n{mutated}")

        # Save analysis report
        report_path = f"{output_dir}/{accession_id}_report.json"
        with open(report_path, "w") as f:
            json.dump(report, f, indent=4)

        print(f"Analysis report saved at: {report_path}")
        analysis_report = SequenceAnalysisReport(original, mutated, seq_type)
        report = analysis_report.compare_sequences()
        analysis_report.plot_comparison(report, save_path=f"{output_dir}/{accession_id}_plots")

        print("Operation completed successfully!")
        logging.info("Pipeline completed")

    except Exception as e:
        logging.error(f"Critical error: {str(e)}")
        print(f"Error: {str(e)}")

def is_ncbi_accession(accession_id: str) -> bool:
    return bool(re.match(r"^[A-Z]{2}_?\d+", accession_id))

if __name__ == "__main__":
    main()