import matplotlib.pyplot as plt
import numpy as np


def save_plot(s1, s2, filepath, figsize=(10, 8), cmap='RdYlGn'):
    """
    Generate a heatmap showing sequence similarity.
    
    Args:
        s1: First sequence
        s2: Second sequence
        filepath: Output file path
        figsize: Figure size (width, height)
        cmap: Colormap name ('RdYlGn', 'viridis', 'YlOrRd', etc.)
    """
    n, m = len(s1), len(s2)
    
    # Create similarity matrix
    similarity_matrix = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            # 1 for match, 0 for mismatch
            similarity_matrix[i, j] = 1 if s1[i] == s2[j] else 0
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=figsize)
    
    im = ax.imshow(similarity_matrix, cmap=cmap, aspect='auto', origin='upper', interpolation='nearest')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label='Match (1) / Mismatch (0)')
    
    # Set labels
    ax.set_xlabel('Sequence 2 Position', fontsize=12)
    ax.set_ylabel('Sequence 1 Position', fontsize=12)
    ax.set_title('Sequence Similarity Heatmap', fontsize=14, fontweight='bold')
    
    # Add grid for better readability (optional, for smaller sequences)
    if n < 100 and m < 100:
        ax.set_xticks(np.arange(-0.5, m, 10), minor=True)
        ax.set_yticks(np.arange(-0.5, n, 10), minor=True)
        ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filepath, dpi=150, bbox_inches='tight')
    plt.close()
