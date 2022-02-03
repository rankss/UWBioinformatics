# Constants/Defaults

class Block:
    NUCLEOTIDES = "ACGT"
    AMINO_ACIDS = "CSTPAGNDEQHRKMILVFYW"

class Matrix:
    EXAMPLE_MATRIX = {
        'A': {
            'A': 1,
            'G': 0,
            'C': 0,
            'T': 0
        },
        'G': {
            'A': 0,
            'G': 1,
            'C': 0,
            'T': 0
        },
        'C': {
            'A': 0,
            'G': 0,
            'C': 1,
            'T': 0
        },
        'T': {
            'A': 0,
            'G': 0,
            'C': 0,
            'T': 1
        },
    }

    BLOSUM62 = {
        # TODO
    }
    