"""
Molecule constants
"""
ALLOWED_ATOMS = [1, 6, 7, 8, 9, 15, 16, 17, 35, 53]
SALT_METALS = [3, 11, 12, 19, 20]

"""
Protein constants
"""
THREE_LETTER_AMINO_ACID_CODES = (
    'ALA',
    'ARG',
    'ASN',
    'ASP',
    'ASX',  # asparagine or aspartic acid
    'CYS',
    'GLU',
    'GLN',
    'GLX',  # glutamine or glutamic acid
    'GLY',
    'HIS',
    'ILE',
    'LEU',
    'LYS',
    'MET',
    'PHE',
    'PRO',
    'SER',
    'THR',
    'TRP',
    'TYR',
    'VAL'
)


SINGLE_LETTER_AMINO_ACID_CODES = (  # same order as three letter codes
    'A',  # alanine
    'R',  # arginine
    'N',  # asparagine
    'D',  # aspartic acid
    'B',  # aspartic acid or asparagine
    'C',  # cysteine
    'E',  # glutamic acid
    'Q',  # glutamine
    'Z',  # glutamic acid or glutamine
    'G',  # glycine
    'H',  # histidin
    'I',  # isoleucine
    'L',  # leucine
    'K',  # lysine
    'M',  # methionine
    'F',  # phenylalanine
    'P',  # proline
    'S',  # serine
    'T',  # threonine
    'W',  # tryptophan
    'Y',  # tyrosine
    'V',  # valine
)
