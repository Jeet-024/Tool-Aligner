# Tool-Aligner
# This repository includes the source code for a CLI based tool built on Python language for aligning DNA/RNA sequences directly given as input or using fasta files. It provides a quick basic alignment result in a textual format and dotplot. 


Simple command-line tool to align DNA/RNA sequences using Needleman-Wunsch (global) and Smith-Waterman (local) algorithms, produce alignment results, identity percentage, and dot-plot images.

Usage example:

```
python run_aligner.py --seq1 ATGCTAG --seq2 ATGATAG --algorithm needleman-wunsch --outdir output
```

Outputs:
- `output/alignment.txt` — raw alignment and scores ->
```
  cat alignment.txt
```
- `output/dotplot.png` — dot-plot image ->
```
  eog dotplot.png
```

Dependencies (library required): "numpy" and "matplotlib"
Other dependencies (system requirements): WSL (Ubuntu) and Python 3.13 and above.
