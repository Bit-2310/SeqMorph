# **SeqMorph**

SeqMorph is a bioinformatics tool for simulating mutations in DNA, RNA, and protein sequences. Researchers can explore genetic variations and their impact through advanced mutation simulations, sequence analysis, and statistical testing.

---

## **Key Features**

### 🧬 **Sequence Analysis & Mutation Simulation**

- **Automatic Sequence Type Detection** (DNA, RNA, Protein)
- **Customizable Mutation Types**:
  - Substitutions (Point Mutations with transition/transversion probabilities)
  - Insertions, Deletions, and Frameshifts
  - Translocations, Inversions, and Duplications
- **Batch Processing** for multiple sequences
- **Biological Significance Tests**:
  - Synonymous vs. Non-synonymous mutations
  - GC-content comparison using **t-tests**
  - Nucleotide/K-mer frequency differences via **Chi-square tests**
- **Performance Improvements**:
  - **Optimized algorithms** making SeqMorph **twice as fast** as previous versions
- **Dynamic Input Options**:
  - Fetch sequences via **NCBI** or **UniProt** Accession IDs
  - Load sequences from `.fasta` or `.txt` files

### 📊 **Data Visualization & Export**

- Graphical analysis of:
  - **GC-content comparison**
  - **Mutation density plots**
  - **K-mer frequency distributions**
- **Export Options**:
  - Save mutated sequences and reports in **FASTA, JSON, CSV**
  - Log **mutation positions, types, and original vs. mutated sequences**

---

## **Installation**

SeqMorph requires **Python 3.8+**. To install dependencies:

```bash
pip install biopython numpy scipy matplotlib
```

Alternatively, use:

```bash
pip install -r requirements.txt
```

---

## **Usage**

### **1️⃣ Running SeqMorph**

```bash
python src/SeqMorph_Main.py
```

### **2️⃣ Example Workflow**

#### **Step 1: Load a Sequence**

- Enter an NCBI or UniProt **Accession ID** (e.g., `NC_000913`)
- Or **upload a FASTA file** via interactive file selection

#### **Step 2: Simulate Mutations**

1. Select **mutation type** (`point`, `insertion`, `deletion`, etc.)
2. Define **mutation rate** (e.g., 5%)
3. Set the **mutation range** (e.g., `1-100`)
4. SeqMorph applies **context-aware** Markov-based mutations

#### **Step 3: Analyze & Export**

- View **mutation logs** (positions & types)
- Export **mutated sequences in FASTA**
- Generate **mutation density & sequence comparison plots**
- Export detailed **mutation analysis reports** (JSON, CSV)

---

## **Debugging & Testing**

### ✅ **Check Dependencies**

If you get a `ModuleNotFoundError`, install missing packages:

```bash
pip install biopython numpy scipy matplotlib
```

### 🔬 **Run Unit Tests**

To verify SeqMorph’s functionality, execute:

```bash
python -m unittest discover tests/
```

### 🛠 **Troubleshooting**

- **Missing Biopython?** Run `pip install biopython`
- **Mutation Engine Issues?** Ensure **input sequences** are valid DNA/RNA
- **Incorrect ORF Detection?** Update `gene_library.py`

---

## **Planned Features**

- **Improved Sequence Alignment Tools**
- **More Advanced Mutation Modeling**
- **CLI implementation**

---

## **Contributing**

Contributions are welcome! To contribute:

1. **Fork** the repository
2. **Create a new branch** (`feature-xyz`)
3. **Submit a pull request** after testing

---

## **License**

This project is licensed under the **MIT License**. See the `LICENSE` file for details.
