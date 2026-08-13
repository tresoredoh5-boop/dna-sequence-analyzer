# DNA Sequence Analyzer

A simple Python tool to parse FASTA files and perform basic DNA sequence analysis.

## Scientific question

Can a lightweight, dependency-free Python script correctly reproduce fundamental
molecular biology operations (GC content, transcription, translation) on DNA sequences?

## Data

Small DNA sequences in FASTA format (`test.fasta`), used to validate the script's
outputs against manual biological calculations.

## Method

1. Parse a FASTA file into sequence records
2. Calculate GC content (%)
3. Transcribe DNA into RNA (T → U)
4. Translate RNA into a protein sequence using the standard genetic code,
   stopping at the first stop codon

## How to run

\`\`\`
python bio.py
\`\`\`

## Example output

\`\`\`
sequence_2 -> ATGGGCTATAGCTA
GC content: 42.86 %
ARN: AUGGGCUAUAGCUA
Protéine: MGYS
\`\`\`

## Limitations

- Single reading frame only (no frameshift handling)
- No handling of ambiguous bases (N, etc.)
- No support for reverse complement (planned)

## Next steps

- Add reverse complement function
- Add motif search
- Compare outputs with Biopython for validation

## Author

Built as part of an ongoing bioinformatics learning path (biology background,
self-taught programming), documented step by step.