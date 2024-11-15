
# SeqMorph

SeqMorph is a powerful tool for simulating mutations in DNA, RNA, and protein sequences. Designed for bioinformatics enthusiasts and researchers, SeqMorph helps you explore how mutations impact sequences, their biological functions, and downstream analyses.

---

## **Key Features**

- **Sequence Type Detection**: Automatically detects whether a sequence is DNA, RNA, or protein.
- **Custom Mutations**: Simulate random or region-specific mutations in your sequences.
- **Sense/Antisense Detection**: Identify whether RNA strands are sense or antisense.
- **Sequence Analysis**:
  - Calculate GC content for DNA and RNA.
  - Reverse transcription for RNA.
  - Complement and reverse complement for DNA.
  - Hydrophobicity scoring for proteins.
- **Interactive File Management**: Use a file manager to select sequence files dynamically.
- **Extensibility**: Modular design allows easy addition of new features and tools.

---

## **How It Works**

1. **Input Handling**:
   - Upload `.fasta` or `.txt` sequence files.
   - Automatically validate and parse the sequences.

2. **Mutation Simulation**:
   - Introduce substitutions, insertions, deletions, or frameshifts.
   - Specify mutation rates or apply mutations to specific regions.

3. **Analysis**:
   - Explore sequence-level features (e.g., GC content, complement).
   - Visualize mutation density or sequence changes.

4. **Export**:
   - Save mutated sequences and logs for further analysis.

---

## **Installation**

SeqMorph requires Python 3.8 or above. To get started:

1. Clone the repository:
   ```bash
   git clone https://github.com/your-repo/seqmorph.git
   cd seqmorph
   ```
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

---

## **Usage**

### **Run the Tool**
```bash
python src/seqmorph_main.py
```

### **Example Input Files**
1. `test.fasta`:
   ```
   >Example DNA Sequence
   ATGCGTACGTTAGC
   ```
2. `test.txt`:
   ```
   AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA
   ```
3. `protein_test.txt`:
   ```
   MAMAPRTEINSTRING
   ```

### **Features in Action**
- Select a file using the interactive file manager.
- Detect the sequence type (DNA, RNA, or Protein).
- Perform sequence-specific operations, such as:
  - DNA: Complement, reverse complement, and GC content.
  - RNA: Reverse transcription, GC content, sense/antisense strand detection.
  - Protein: Hydrophobicity scoring.

---

## **Planned Features**

- Advanced mutation engine for precise simulations.
- Visualizations for mutation density and sequence alignments.
- CLI support for batch processing.
- Integration with databases like NCBI GenBank and UniProt.

---

## **Contributing**

Contributions are welcome! Feel free to fork the repository, make changes, and submit a pull request. Please ensure your code follows the established coding guidelines and includes relevant tests.

---

## **License**

This project is licensed under the MIT License. See the `LICENSE` file for details.
