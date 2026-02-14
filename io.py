# Parsing FASTA file for the alignment execution.
def parse_fasta(path):
    """Parse a single-sequence FASTA file and return the sequence string."""
    seq = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                continue
            seq.append(line)
    return ''.join(seq)

