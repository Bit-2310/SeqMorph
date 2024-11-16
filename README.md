# SeqMorph
SeqMorph is a bioinformatics tool designed for simulating mutations in DNA, RNA, and protein sequences. It allows users to explore the impact of mutations on sequences and their biological functions. With an interactive interface, batch processing support, and detailed analysis features, SeqMorph is ideal for researchers and bioinformatics enthusiasts.

---

## **Key Features**

- **Sequence Type Detection**: Automatically detects whether the input is DNA, RNA, or protein.
- **Custom Mutations**: Perform substitutions, insertions, deletions, or frameshift mutations.
- **Unified Gene Library**: Handles both DNA and RNA sequence operations, reducing duplicate code.
- **Custom Mutations**: Perform substitutions, insertions, deletions, frameshift mutations, and point mutations with transition/transversion probability weights.
- **Batch Processing**: Apply mutations to multiple sequences in a single run.
- **Sense/Antisense Detection**: Identify whether RNA strands are sense or antisense.
- **Sequence Analysis**:
  - Calculate GC content for DNA and RNA.
  - Reverse transcription for RNA.
  - Complement and reverse complement for DNA.
  - Hydrophobicity scoring for proteins.
- **Dynamic File Input**:
  - Select `.fasta` or `.txt` files through an interactive file manager.
  - Retrieve sequences from NCBI or UniProt using accession numbers.
  - Enter sequences manually.
- **Mutation Attributes**: Track mutation logs with details such as positions, types, and original vs mutated bases.
- **Export**:
  - Save detailed outputs, including sequence attributes and mutation logs, in a structured file format.

---

## **How It Works**

1. **Input Handling**:
   - Choose from multiple input methods: file upload, accession number, or manual entry.
   - Validate and parse sequences automatically.

2. **Mutation Simulation**:
   - Specify mutation type (`substitute`, `insert`, `delete`, `frameshift`, `silent`, `missense`, `nonsense`).
   - Define mutation rate (default: 10%).
   - Apply mutations to full sequences or specific regions (start index and end index).

3. **Analysis**:
   - Explore sequence-level features:
     - GC content.
     - Nucleotide frequency.
     - Hydrophobicity.
     - Transition/transversion ratio.
   - Highlight mutations compared to the original sequence.

4. **Export**:
   - Save original and mutated sequences, along with their attributes, in a detailed file.

---

## **Installation**

SeqMorph requires Python 3.8 or higher. To install:

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
python src/SeqMorph\ Main.py
```

### **Example Input Files**
1. `example_dna.fasta`:
   ```
   >Example DNA Sequence
   ATGCGTACGTTAGC
   ```
2. `example_rna.txt`:
   ```
   AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA
   ```
3. `example_protein.txt`:
   ```
   MAMAPRTEINSTRING
   ```

### **Features in Action**

1. Select a sequence source:
   - File upload.
   - Manual sequence entry.
2. Detect sequence type (DNA, RNA, or protein).
3. Choose mutation type, rate, and target regions.
4. Analyze sequence attributes:
   - DNA/RNA: GC content, nucleotide frequency, complement, reverse complement, transition/transversion counts.
   - Protein: Hydrophobicity score.
5. Export results:
   - Save the original and mutated sequences, along with their attributes and mutation logs, to a structured output file.

---

## **Planned Features**

- Visualizations for mutation density and sequence changes.
- Advanced alignment tools for comparing original and mutated sequences.
- Integration with online databases for automated data retrieval (*Currently W.I.P*).
- Enhanced reporting with customizable export formats (e.g., CSV, JSON).

---

## **Contributing**

Contributions are welcome! Fork the repository, create a new branch, make your changes, and submit a pull request. Ensure that your code adheres to the project's coding standards and includes relevant tests.

---

## **License**

This project is licensed under the MIT License. See the `LICENSE` file for details.
