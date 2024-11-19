class ProteinUtils:
    @staticmethod
    def hydrophobicity(protein_sequence):
        """
        Calculate the average hydrophobicity of a protein sequence using the Kyte-Doolittle scale.
        """
        hydrophobicity_scale = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }

        score = sum(hydrophobicity_scale.get(aa, 0) for aa in protein_sequence)
        return score / len(protein_sequence) if protein_sequence else 0

    @staticmethod
    def molecular_weight(protein_sequence):
        """
        Calculate the approximate molecular weight of a protein sequence in Daltons (Da).
        Uses average weights of amino acids.
        """
        amino_acid_weights = {
            'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
            'Q': 146.2, 'E': 147.1, 'G': 75.1, 'H': 155.2, 'I': 131.2,
            'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
            'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
        }

        weight = sum(amino_acid_weights.get(aa, 0) for aa in protein_sequence)
        # Adding the weight of a water molecule for peptide bonds
        weight += (len(protein_sequence) - 1) * 18.015 if protein_sequence else 0
        return weight

    @staticmethod
    def charge(protein_sequence, pH=7.4):
        """
        Calculate the net charge of a protein sequence at a given pH.
        Uses pKa values for ionizable amino acids.
        """
        pKa_values = {
            'D': 3.9,  # Aspartic acid
            'E': 4.1,  # Glutamic acid
            'H': 6.0,  # Histidine
            'C': 8.3,  # Cysteine
            'Y': 10.1,  # Tyrosine
            'K': 10.5,  # Lysine
            'R': 12.5   # Arginine
        }

        # Terminal groups
        n_term_pKa = 9.6
        c_term_pKa = 2.1

        # Calculate the contributions to the charge
        charge = 0
        for aa in protein_sequence:
            if aa in pKa_values:
                if aa in 'DE':  # Negatively charged
                    charge -= 1 / (1 + 10**(pH - pKa_values[aa]))
                elif aa in 'HKRYC':  # Positively charged
                    charge += 1 / (1 + 10**(pKa_values[aa] - pH))

        # Add terminal contributions
        charge += 1 / (1 + 10**(pH - n_term_pKa))  # N-terminal
        charge -= 1 / (1 + 10**(c_term_pKa - pH))  # C-terminal

        return charge

    @staticmethod
    def isoelectric_point(protein_sequence):
        """
        Estimate the isoelectric point (pI) of a protein sequence.
        Uses a bisection approach to approximate the pH where net charge is 0.
        """
        lower_pH = 0.0
        upper_pH = 14.0
        tolerance = 0.01

        while upper_pH - lower_pH > tolerance:
            mid_pH = (lower_pH + upper_pH) / 2.0
            net_charge = ProteinUtils.charge(protein_sequence, pH=mid_pH)
            if net_charge > 0:
                lower_pH = mid_pH
            else:
                upper_pH = mid_pH

        return (lower_pH + upper_pH) / 2.0

    @staticmethod
    def validate(sequence):
        """
        Validate the protein sequence.
        """
        valid_bases = set("ACDEFGHIKLMNPQRSTVWYX")
        return set(sequence.upper()).issubset(valid_bases)
