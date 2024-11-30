
# **SeqMorph**
SeqMorph is an advanced bioinformatics tool tailored for simulating mutations in DNA, RNA, and protein sequences. 
It empowers researchers to investigate the effects of genetic variations on sequences and their biological functions through highly customizable mutation simulations and comprehensive analyses.
---

## **Key Features**

- **Sequence Type Detection**: Automatically identifies whether the input is DNA, RNA, or protein.
- **Custom Mutations**: Simulate a variety of mutations including:
  - Substitutions
  - Insertions
  - Deletions
  - Frameshifts
  - Point mutations with transition/transversion probability weights
  - Translocations
  - Inversions
  - Duplications
- **Batch Processing**: Apply mutations to multiple sequences in one run.
- **Unified Gene Library**: Simplifies sequence operations for DNA and RNA.
- **Sense/Antisense Detection**: Detect whether RNA strands are sense or antisense.
- **Sequence Analysis**:
  - Calculate GC content for DNA/RNA.
  - Reverse transcription for RNA.
  - Complement and reverse complement for DNA.
  - Hydrophobicity scoring for proteins.
  - Synonymous vs. non-synonymous mutation analysis for protein-coding regions.
- **Dynamic Input Options**:
  - Retrieve sequences using NCBI or UniProt accession numbers.
  - Select `.fasta` or `.txt` files through an interactive file manager.
- **Mutation Attributes**: Logs details such as mutation positions, types, and original vs. mutated bases.
- **Export**:
  - Save mutated sequences and detailed mutation logs in FASTA and text formats.
  - Generate visualizations for sequence comparisons and mutation density.

---

## **How It Works**

### **Input Handling**
- Fetch sequences using NCBI or UniProt accession numbers.
- Upload sequences via `.fasta` or `.txt` files.
- Detect sequence types automatically (DNA, RNA, or Protein).

### **Mutation Simulation**
1. Define mutation rate (default: 5%).
2. Select mutation range (start and end positions).
3. Choose mutation type:
   - Point Mutation
   - Substitution
   - Insertion
   - Deletion
   - Translocation
   - Inversion
   - Duplication

### **Analysis**
- Generate a mutation report with synonymous and non-synonymous counts for protein sequences.
- Calculate GC content, nucleotide frequency, and hydrophobicity scores.
- Highlight mutation density across sequences.

### **Export**
- Save original and mutated sequences with metadata.
- Export mutation analysis reports and visual comparisons.

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

### **Example Workflow**
1. **Input Handling**:
   - Enter an NCBI or UniProt accession ID (e.g., `NC_000913`).
   - Sequence is automatically fetched and validated.

2. **Mutation Simulation**:
   - Define the mutation rate (e.g., 5%).
   - Specify mutation range (e.g., positions 1-100).
   - Choose the mutation type (e.g., substitution).

3. **Analysis and Export**:
   - View mutation logs, including positions and types.
   - Save mutated sequences in `.fasta` format.
   - Export detailed mutation analysis reports.

---

## **Planned Features**
- Manual sequence input via the terminal.
- Visualizations for mutation density and sequence comparisons.
- Advanced sequence alignment tools.
- Customizable export formats (e.g., CSV, JSON).

---

## **Contributing**

We welcome contributions! Fork the repository, create a new branch, make your changes, and submit a pull request. Ensure your code adheres to the project's standards and includes relevant tests.

---

## **License**

This project is licensed under the MIT License. See the `LICENSE` file for details.

---
