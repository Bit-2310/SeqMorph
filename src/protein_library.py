class ProteinUtils:
    # Class-level constants for better reuse
    HYDROPHOBICITY_SCALE = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }

    AMINO_ACID_WEIGHTS = {
        'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
        'Q': 146.2, 'E': 147.1, 'G': 75.1, 'H': 155.2, 'I': 131.2,
        'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
        'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
    }

    PKA_VALUES = {
        'D': 3.9, 'E': 4.1, 'H': 6.0, 'C': 8.3,
        'Y': 10.1, 'K': 10.5, 'R': 12.5
    }

    N_TERM_PKA = 9.6
    C_TERM_PKA = 2.1

    @staticmethod
    def hydrophobicity(protein_sequence: str) -> float:
        """Calculate average hydrophobicity (Kyte-Doolittle scale)."""
        if not protein_sequence:
            raise ValueError("Protein sequence cannot be empty")
        scores = [ProteinUtils.HYDROPHOBICITY_SCALE.get(aa, 0) 
                  for aa in protein_sequence]
        return sum(scores) / len(protein_sequence)

    @staticmethod
    def molecular_weight(protein_sequence: str) -> float:
        """Calculate molecular weight (Da), subtracting water for peptide bonds."""
        if not protein_sequence:
            raise ValueError("Protein sequence cannot be empty")
        
        # Sum individual amino acid weights
        total = sum(ProteinUtils.AMINO_ACID_WEIGHTS.get(aa, 0) 
                    for aa in protein_sequence)
        
        # Subtract water for each peptide bond: (n-1) bonds × 18.015 Da
        water_loss = (len(protein_sequence) - 1) * 18.015
        return total - water_loss

    @staticmethod
    def charge(protein_sequence: str, pH: float = 7.4) -> float:
        """Calculate net charge using pre-mapped contributions."""
        charge_dict = {
            'D': (-1, ProteinUtils.PKA_VALUES['D']),
            'E': (-1, ProteinUtils.PKA_VALUES['E']),
            'H': (+1, ProteinUtils.PKA_VALUES['H']),
            'C': (+1, ProteinUtils.PKA_VALUES['C']),
            'Y': (+1, ProteinUtils.PKA_VALUES['Y']),
            'K': (+1, ProteinUtils.PKA_VALUES['K']),
            'R': (+1, ProteinUtils.PKA_VALUES['R'])
        }

        charge = 0
        for aa in protein_sequence:
            if aa in charge_dict:
                sign, pKa = charge_dict[aa]
                charge += sign / (1 + 10**(sign * (pH - pKa)))

        # Add N-terminal (amine) and C-terminal (carboxyl) charges
        charge += 1 / (1 + 10**(pH - ProteinUtils.N_TERM_PKA))  # N-term
        charge -= 1 / (1 + 10**(ProteinUtils.C_TERM_PKA - pH))  # C-term
        return charge

    @staticmethod
    def isoelectric_point(protein_sequence: str, 
                         tolerance: float = 0.01, 
                         max_iter: int = 100) -> float:
        """Optimized pI estimation with convergence checks."""
        lower, upper = 2.0, 12.0  # Narrowed initial range
        for _ in range(max_iter):
            mid = (lower + upper) / 2
            current_charge = ProteinUtils.charge(protein_sequence, mid)
            if abs(current_charge) < tolerance:
                return mid
            if current_charge > 0:
                lower = mid
            else:
                upper = mid
        return (lower + upper) / 2

    @staticmethod
    def validate(sequence: str) -> bool:
        """Validate protein sequence (allows 'X' for unknown residues)."""
        valid = set("ACDEFGHIKLMNPQRSTVWYX")
        return all(aa.upper() in valid for aa in sequence)