# Constants/Defaults

class BlockConstants:

    NUCLEOTIDES = "ACGT"
    AMINO_ACIDS = "CSTPAGNDEQHRKMILVFYW"

class MatrixConstants:

    EXAMPLE_MATRIX = {
        'A': {
            'A': 1,
            'G': -1,
            'C': -1,
            'T': -1
        },
        'G': {
            'A': -1,
            'G': 1,
            'C': -1,
            'T': -1
        },
        'C': {
            'A': -1,
            'G': -1,
            'C': 1,
            'T': -1
        },
        'T': {
            'A': -1,
            'G': -1,
            'C': -1,
            'T': 1
        },
    }

    BLOSUM62 = {
        # TODO
    }
    