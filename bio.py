def read_fasta(filepath):
    sequences = {}
    current_name = None
    current_seq = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line == "":
                continue

            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:]
                current_seq = []
            else:
                current_seq.append(line.upper())

        if current_name is not None:
            sequences[current_name] = "".join(current_seq)

    return sequences


def gc_content(sequence):
    """
    Calcule le pourcentage de bases G et C dans une séquence ADN.
    """
    sequence = sequence.upper()
    total = len(sequence)
    if total == 0:
        return 0.0

    g_count = sequence.count("G")
    c_count = sequence.count("C")

    gc_percent = (g_count + c_count) / total * 100
    return round(gc_percent, 2)


def transcribe(sequence):
    """
    Transcrit une séquence ADN en ARN (remplace T par U).
    """
    sequence = sequence.upper()
    return sequence.replace("T", "U")


CODON_TABLE = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(rna_sequence):
    """
    Traduit une séquence ARN en protéine, codon par codon.
    S'arrête au premier codon stop (*).
    """
    protein = ""
    for i in range(0, len(rna_sequence) - 2, 3):
        codon = rna_sequence[i:i+3]
        amino_acid = CODON_TABLE.get(codon, "?")
        if amino_acid == "*":
            break
        protein += amino_acid
    return protein


def reverse_complement(sequence):
    """
    Retourne le brin complémentaire inverse d'une séquence ADN.
    """
    sequence = sequence.upper()
    complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}

    complement = "".join(complement_map.get(base, "N") for base in sequence)
    reversed_complement = complement[::-1]

    return reversed_complement


def find_motif(sequence, motif):
    """
    Retourne la liste des positions (index 0-based) où le motif
    apparaît dans la séquence. Gère les occurrences qui se chevauchent.
    """
    sequence = sequence.upper()
    motif = motif.upper()
    positions = []

    for i in range(len(sequence) - len(motif) + 1):
        if sequence[i:i + len(motif)] == motif:
            positions.append(i)

    return positions


def describe(values):
    """
    Calcule moyenne, médiane, variance et écart-type d'une liste de nombres.
    """
    n = len(values)
    if n == 0:
        return None

    mean = sum(values) / n

    sorted_values = sorted(values)
    mid = n // 2
    if n % 2 == 0:
        median = (sorted_values[mid - 1] + sorted_values[mid]) / 2
    else:
        median = sorted_values[mid]

    variance = sum((x - mean) ** 2 for x in values) / n
    std_dev = variance ** 0.5

    return {
        "mean": round(mean, 2),
        "median": round(median, 2),
        "variance": round(variance, 2),
        "std_dev": round(std_dev, 2),
    }


seqs = read_fasta("test.fasta")
for name, seq in seqs.items():
    print(name, "->", seq)
    print("GC content:", gc_content(seq), "%")
    rna = transcribe(seq)
    print("ARN:", rna)
    print("Protéine:", translate(rna))
    print("Brin complémentaire inverse:", reverse_complement(seq))
    print("Positions du motif 'ATG':", find_motif(seq, "ATG"))

gc_values = [gc_content(seq) for seq in seqs.values()]
print("\nStatistiques sur le GC content de toutes les séquences:")
print(describe(gc_values))