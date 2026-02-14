import argparse
import os
from aligner.algorithms import needleman_wunsch, smith_waterman
from aligner.dotplot import save_dotplot
from aligner.io import parse_fasta


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
    p = argparse.ArgumentParser(description='Simple DNA/RNA aligner')
    p.add_argument('--seq1', help='Direct sequence input for sequence 1')
    p.add_argument('--seq2', help='Direct sequence input for sequence 2')
    p.add_argument('--fasta1', help='FASTA file for sequence 1')
    p.add_argument('--fasta2', help='FASTA file for sequence 2')
    p.add_argument('--algorithm', choices=['needleman-wunsch', 'smith-waterman'], default='needleman-wunsch')
    p.add_argument('--outdir', default='output')
    args = p.parse_args()

    # Load sequences from files or direct input
    if args.fasta1:
        seq1 = parse_fasta(args.fasta1)
    else:
        seq1 = args.seq1 or ''
    if args.fasta2:
        seq2 = parse_fasta(args.fasta2)
    else:
        seq2 = args.seq2 or ''

    # Validate sequences
    if not seq1 or not seq2:
        p.error('Provide sequences via --seq1/--seq2 or --fasta1/--fasta2')

    os.makedirs(args.outdir, exist_ok=True)

    # Run alignment algorithm
    if args.algorithm == 'needleman-wunsch':
        res = needleman_wunsch(seq1, seq2)
    else:
        res = smith_waterman(seq1, seq2)

    # Write results
    out_align = os.path.join(args.outdir, 'alignment.txt')
    write_alignment(out_align, res, seq1, seq2, args.algorithm)

    out_dot = os.path.join(args.outdir, 'dotplot.png')
    save_dotplot(seq1, seq2, out_dot)

    print('Wrote:', out_align)
    print('Wrote:', out_dot)


if __name__ == '__main__':
    main()
