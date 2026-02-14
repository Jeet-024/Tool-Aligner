import matplotlib.pyplot as plt


def save_dotplot(s1, s2, filepath, dot_color='black', figsize=(6, 6)):
    n, m = len(s1), len(s2)
    xs, ys = [], []
    for i, a in enumerate(s1):
        for j, b in enumerate(s2):
            if a == b:
                xs.append(j)
                ys.append(n - 1 - i)

    plt.figure(figsize=figsize)
    plt.scatter(xs, ys, s=6, c=dot_color)
    plt.xlim(-0.5, m - 0.5)
    plt.ylim(-0.5, n - 0.5)
    plt.xlabel('Sequence 2')
    plt.ylabel('Sequence 1')
    plt.title('Dot-plot')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(filepath, dpi=150)
    plt.close()
