#Tool Aligner
#Author: Manojit Mazumder
#Date: 18-02-2026

import argparse
import os
import sys
from aligner.algorithms import needleman_wunsch_affine, smith_waterman_affine
from aligner.heatmap import save_plot
from aligner.io import parse_fasta


def load_sequence(file_path):
    """Load sequence from a file path or stdin ('-').

    Auto-detects FASTA by content (lines starting with '>') and returns
    the first sequence found. For plain text it collapses whitespace.
    """
    # Support reading from stdin when path is '-'
    if file_path == '-':
        content = sys.stdin.read().strip()
    else:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        with open(file_path, 'r') as f:
            content = f.read().strip()

    if not content:
        return ''

    # FASTA content -> parse first record
    if content.startswith('>'):
        # If we read from a file path, prefer the dedicated parser in `aligner.io`
        # (this also ensures the `parse_fasta` import is used so linters/Pylance
        # can resolve the symbol). If reading from stdin, parse the content.
        if file_path != '-':
            return parse_fasta(file_path)

        seqs = []
        current = []
        for line in content.splitlines():
            if line.startswith('>'):
                if current:
                    seqs.append(''.join(current))
                    current = []
            else:
                current.append(line.strip())
        if current:
            seqs.append(''.join(current))

        if not seqs:
            return ''
        return seqs[0]

    # Plain text sequence - collapse whitespace/newlines
    return ''.join(content.split())


def write_alignment(outpath, result, seq1, seq2, algorithm):
    with open(outpath, 'w') as f:
        f.write(f"Algorithm: {algorithm}\n")
        f.write(f"Score: {result['score']}\n")
        f.write(f"Matches: {result['matches']} / {result['alignment_length']}\n")
        f.write(f"Identity: {result['identity_percent']:.2f}%\n\n")
        f.write(result['alignment1'] + '\n')
        f.write(result['alignment2'] + '\n')


def main():
    # Parse command-line arguments
    p = argparse.ArgumentParser(description='Simple DNA/RNA aligner', 
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog="""
Examples:
  python run_aligner.py --file1 seq1.fasta --file2 seq2.fasta
  python run_aligner.py --file1 seq1.txt --file2 seq2.txt
  python run_aligner.py --seq1 ATCG --seq2 ATCCG (direct input, for testing)
                                """)
    
    # File input arguments (preferred method)
    p.add_argument('--file1', help='Input file for sequence 1 (FASTA or plain text)')
    p.add_argument('--file2', help='Input file for sequence 2 (FASTA or plain text)')
    
    # Legacy direct sequence input arguments
    p.add_argument('--seq1', help='Direct sequence input for sequence 1 (for testing)')
    p.add_argument('--seq2', help='Direct sequence input for sequence 2 (for testing)')
    
    # Legacy FASTA arguments (still supported)
    p.add_argument('--fasta1', help='FASTA file for sequence 1 (deprecated, use --file1)')
    p.add_argument('--fasta2', help='FASTA file for sequence 2 (deprecated, use --file1)')
    
    p.add_argument('--algorithm', choices=['needleman-wunsch', 'smith-waterman'], default='needleman-wunsch')
    p.add_argument('--outdir', default='output')
    args = p.parse_args()

    # Load sequences with priority: file > fasta (legacy) > direct input
    seq1 = None
    seq2 = None
    
    # Try file input first
    if args.file1:
        try:
            seq1 = load_sequence(args.file1)
        except Exception as e:
            p.error(f"Error reading file1: {e}")
    
    if args.file2:
        try:
            seq2 = load_sequence(args.file2)
        except Exception as e:
            p.error(f"Error reading file2: {e}")
    
    # Fall back to legacy FASTA arguments
    if not seq1 and args.fasta1:
        try:
            seq1 = load_sequence(args.fasta1)
        except Exception as e:
            p.error(f"Error reading fasta1: {e}")
    
    if not seq2 and args.fasta2:
        try:
            seq2 = load_sequence(args.fasta2)
        except Exception as e:
            p.error(f"Error reading fasta2: {e}")
    
    # Fall back to direct input
    if not seq1:
        seq1 = args.seq1 or ''
    if not seq2:
        seq2 = args.seq2 or ''

    # Validate sequences
    if not seq1 or not seq2:
        p.error('Provide sequences via --file1/--file2 or --fasta1/--fasta2 or --seq1/--seq2')


    os.makedirs(args.outdir, exist_ok=True)

    # Run alignment algorithm
    if args.algorithm == 'needleman-wunsch':
        res = needleman_wunsch_affine(seq1, seq2)
    else:
        res = smith_waterman_affine(seq1, seq2)

    # Write results
    out_align = os.path.join(args.outdir, 'alignment.txt')
    write_alignment(out_align, res, seq1, seq2, args.algorithm)

    out_dot = os.path.join(args.outdir, 'heatmap.png')
    save_plot(seq1, seq2, out_dot)

    print('Wrote:', out_align)
    print('Wrote:', out_dot)


if __name__ == '__main__':
    main()
