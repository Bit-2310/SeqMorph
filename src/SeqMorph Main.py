from input_module import InputHandler
from mutation_engine import MutationEngine
from export_module import ExportHandler  # Ensure this is imported

def main():
    print("Welcome to SeqMorph: A Sequence Mutation Simulator")

    # Input Handling
    handler = InputHandler()

    # Get sequence input
    print("How would you like to input sequences?")
    print("1. Load from files (batch processing)")
    print("2. Enter a single sequence manually")
    choice = input("Enter your choice (1 or 2): ").strip()

    sequences = []
    if choice == "1":
        print("Select input files with sequences.")
        file_paths = handler.select_multiple_files()
        if not file_paths:
            print("No files selected. Exiting.")
            return
        try:
            for file_path in file_paths:
                sequences.append(handler.load_sequence_file(file_path))
        except Exception as e:
            print(f"Error loading files: {e}")
            return
    elif choice == "2":
        sequence = input("Enter your sequence: ").strip()
        if not sequence:
            print("No sequence entered. Exiting.")
            return
        sequences.append(sequence)
    else:
        print("Invalid choice. Exiting.")
        return

    for idx, sequence in enumerate(sequences, 1):
        print(f"\nProcessing sequence {idx}: {sequence}")
        try:
            sequence_type = handler._detect_sequence_type(sequence)
            print(f"Detected sequence type: {sequence_type}")
        except Exception as e:
            print(f"Error detecting sequence type: {e}")
            continue

        # Ask how many times to run the program for this sequence
        repeat_times = input("How many times do you want to run mutations on this sequence? (default: 1): ").strip()
        repeat_times = int(repeat_times) if repeat_times.isdigit() else 1

        for run in range(1, repeat_times + 1):
            print(f"\nRun {run}/{repeat_times} for sequence {idx}")

            # Mutation Options
            print("Select mutation type:")
            print("1. Substitute")
            print("2. Insert")
            print("3. Delete")
            print("4. Frameshift")
            mutation_choice = input("Enter the number corresponding to the mutation type: ").strip()
            mutation_types = {"1": "substitute", "2": "insert", "3": "delete", "4": "frameshift"}
            mutation_type = mutation_types.get(mutation_choice)
            if not mutation_type:
                print("Invalid choice. Skipping this run.")
                continue

            mutation_rate = input("Enter mutation rate as a percentage (default 10%): ").strip()
            mutation_rate = float(mutation_rate) if mutation_rate else 10.0

            start = input("Enter start position for mutation range (optional): ").strip()
            start = int(start) if start else None

            end = input("Enter end position for mutation range (optional): ").strip()
            end = int(end) if end else None

            try:
                mutated_sequence = MutationEngine.mutate(
                    sequence, mutation_type, mutation_rate, start, end
                )
                print(f"Mutated sequence: {mutated_sequence}")
            except Exception as e:
                print(f"Error during mutation: {e}")
                continue

            # Highlight Mutations
            highlight_choice = input("Do you want to highlight the mutations? (y/n): ").strip().lower()
            if highlight_choice == "y":
                highlighted_sequence = MutationEngine.highlight_mutations(sequence, mutated_sequence)
                print(f"Highlighted Mutations: {highlighted_sequence}")

            # Export Results
            export_choice = input("Do you want to save the mutated sequence to a file? (y/n): ").strip().lower()
            if export_choice == "y":
                output_file = input(f"Enter output file name for sequence {idx}, run {run} (default: mutated_output_{idx}_run{run}.txt): ").strip()
                output_file = output_file if output_file else f"mutated_output_{idx}_run{run}.txt"
                try:
                    ExportHandler.export_sequence(sequence, mutated_sequence, output_file, f"Sequence_{idx}_Run_{run}")
                    print(f"Mutated sequence saved to {output_file}")
                except Exception as e:
                    print(f"Error saving file: {e}")

if __name__ == "__main__":
    main()
