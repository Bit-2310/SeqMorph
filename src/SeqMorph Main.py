import re
import os
import logging
from input_module import InputHandler, SequenceFetcher
from mutation_module import MutationUtils as mUtils, SequenceAnalyzer, CODON_TABLE
from analysis_module import SequenceAnalysisReport as SeqR
from sequence_structure import Sequence_Structure as SS


# Configure logging
logging.basicConfig(
    filename='seqmorph.log',
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def main():
    # Print welcome message with ASCII art
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

    # Step 1: Get the input sequences
    handler = InputHandler()
    fetcher = SequenceFetcher()
    
    input_method = "1"  # Going to include other methods soon (ex: manual input and loading file)
    sequence_data = {}

    try:
        if input_method == "1":
            # Automatically detect whether it is an NCBI or UniProt ID
            accession_id = input("Enter the accession ID (NCBI or UniProt): ").strip()
            logging.info(f"Accession ID provided: {accession_id}")

            # Check if the accession ID matches NCBI or UniProt pattern
            if re.match(r"^(NC_|NM_|NW_)", accession_id):
                # NCBI ID detected
                seq = fetcher.fetch_sequence(accession_id)
                seq_type = handler.sequence_type(seq)
                sequence_data[accession_id] = {'sequence': seq, 'type': seq_type}  # Store the sequence data
                print(f"Sequence fetched from NCBI: {accession_id}")
                logging.info(f"Sequence fetched from NCBI for accession ID: {accession_id}")
            elif re.match(r"^[A-Za-z0-9]{6,}", accession_id):
                # UniProt ID detected (simple alphanumeric check)
                seq = fetcher.fetch_sequence(accession_id)
                seq_type = handler.sequence_type(seq)
                sequence_data[accession_id] = {'sequence': seq, 'type': seq_type}  # Store the sequence data
                print(f"Sequence fetched from UniProt: {accession_id}")
                logging.info(f"Sequence fetched from UniProt for accession ID: {accession_id}")
            else:
                print("Invalid accession ID format. Please provide a valid NCBI or UniProt ID.")
                logging.warning("Invalid accession ID format provided.")
                return  # Exit if the ID format doesn't match expected patterns

        # Step 2: Initialize Sequence_Structure with sequence data
        if sequence_data:
            seqs = SS(sequence_data)

            # Add sequences using `add_sequence` method
            for accession, data in sequence_data.items():
                seqs.add_sequence(accession, data['sequence'], data['type'])
            logging.info("Sequences added to Sequence_Structure.")

            # Step 3: Display Metadata
            if len(seqs.Sequences) != 0:
                print("\nMetadata summary: \n====================== \n", seqs.MetaData)
                print("\nSequence type summary: \n====================== \n", seqs.get_summary())
                logging.info("Displayed metadata summary.")

                # Step 4: Mutation inputs & applying mutation (only one call)
                mutation_utils = mUtils(seqs)
                mutation_rate, range_values, mutation_type = mutation_utils.mutation_inputs(seqs, accession_id)

                if mutation_rate and range_values and mutation_type:
                    original_seq = seqs.get_sequence(accession_id)  # Fetch original sequence for later reference
                    mutated_sequence = mutation_utils.apply_mutation(
                        accession_id, mutation_rate, range_values[0], range_values[1], mutation_type
                    )
                    print("Mutation function completed.")
                    logging.info(f"Mutation applied successfully to sequence {accession_id}.")

                    # Step 5: Generate mutation report using SequenceAnalyzer
                    analyzer = SequenceAnalyzer(CODON_TABLE)
                    synonymous_count, non_synonymous_count = analyzer.analyze_synonymous(original_seq, mutated_sequence)
                    mutation_report = analyzer.generate_mutation_report(synonymous_count, non_synonymous_count)
                    mutation_report['total_mutations'] = len(mutation_utils.get_mutation_log())
                    print("\nMutation Report: \n======================")
                    print(f"Total Mutations: {mutation_report['total_mutations']}")
                    print(f"Synonymous Mutations: {mutation_report['synonymous_mutations']}")
                    print(f"Non-synonymous Mutations: {mutation_report['non_synonymous_mutations']}")
                    logging.info("Mutation report generated.")

                    # Step 6: Create output directory if it doesn't exist
                    output_dir = "output"
                    if not os.path.exists(output_dir):
                        os.makedirs(output_dir)
                        logging.info(f"Output directory created: {output_dir}")

                    # Step 7: Export original and mutated sequence
                    mutated_file_name = os.path.join(output_dir, f"mutated_sequence_{accession_id}.fasta")
                    analysis_file_name = os.path.join(output_dir, f"mutation_analysis_{accession_id}.txt")

                    # Export mutated sequence
                    with open(mutated_file_name, "w") as mutated_file:
                        mutated_file.write(f">Mutated Sequence {accession_id}\n")
                        mutated_file.write(mutated_sequence)
                    print(f"Mutated sequence exported to {mutated_file_name}")
                    logging.info(f"Mutated sequence exported to {mutated_file_name}")

                    # Export mutation analysis report
                    with open(analysis_file_name, "w") as analysis_file:
                        analysis_file.write("Mutation Analysis Report\n======================\n")
                        analysis_file.write(f"Total Mutations: {mutation_report['total_mutations']}\n")
                        analysis_file.write(f"Synonymous Mutations: {mutation_report['synonymous_mutations']}\n")
                        analysis_file.write(f"Non-synonymous Mutations: {mutation_report['non_synonymous_mutations']}\n")
                    print(f"Mutation analysis report exported to {analysis_file_name}")
                    logging.info(f"Mutation analysis report exported to {analysis_file_name}")

                    # Step 8: Generate detailed analysis report and plot comparisons
                    analysis_report = SeqR(original_seq, mutated_sequence, seqs.MetaData[accession_id]['Type'])
                    report = analysis_report.compare_sequences()
                    analysis_report.plot_comparison(report)
                    print(f"Detailed mutation analysis report and visual comparison have been generated.")
                    logging.info("Detailed analysis report and visual comparison generated.")
                else:
                    print("Mutation was not applied.")
                    logging.warning("Mutation inputs were incomplete. Mutation not applied.")
            else:
                print("No valid sequences loaded.")
                logging.warning("No valid sequences loaded to Sequence_Structure.")
        else:
            print("No sequence data available.")
            logging.warning("No sequence data available for processing.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")

# Run the main function
if __name__ == "__main__":
    main()
