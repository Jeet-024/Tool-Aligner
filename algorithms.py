#Needleman-Wunsch and Smith-Waterman algorithm
def _score(a, b, match, mismatch, sub_matrix=None):
    """Return substitution score."""
    if sub_matrix:
        return sub_matrix.get((a, b), sub_matrix.get((b, a), mismatch))
    return match if a == b else mismatch


def _calculate_identity(align1, align2):
    """Identity excluding gap-only columns."""
    matches = 0
    comparable = 0

    for a, b in zip(align1, align2):
        if a != '-' and b != '-':
            comparable += 1
            if a == b:
                matches += 1

    identity = (matches / comparable * 100) if comparable > 0 else 0
    return matches, len(align1), identity


# =========================
# Needleman–Wunsch (Affine)
# =========================

def needleman_wunsch_affine(
    s1,
    s2,
    match=2,
    mismatch=-1,
    gap_open=-2,
    gap_extend=-1,
    sub_matrix=None
):
    """
    Global alignment with affine gap penalty.
    """

    m, n = len(s1), len(s2)

    # Three matrices
    M = [[0]*(n+1) for _ in range(m+1)]  # match/mismatch
    X = [[float('-inf')]*(n+1) for _ in range(m+1)]  # gap in s2
    Y = [[float('-inf')]*(n+1) for _ in range(m+1)]  # gap in s1

    trace = [[None]*(n+1) for _ in range(m+1)]

    # Initialization
    for i in range(1, m+1):
        X[i][0] = gap_open + (i-1)*gap_extend
        M[i][0] = X[i][0]

    for j in range(1, n+1):
        Y[0][j] = gap_open + (j-1)*gap_extend
        M[0][j] = Y[0][j]

    # Fill matrices
    for i in range(1, m+1):
        for j in range(1, n+1):

            score_sub = _score(s1[i-1], s2[j-1], match, mismatch, sub_matrix)

            X[i][j] = max(
                M[i-1][j] + gap_open,
                X[i-1][j] + gap_extend
            )

            Y[i][j] = max(
                M[i][j-1] + gap_open,
                Y[i][j-1] + gap_extend
            )

            M[i][j] = max(
                M[i-1][j-1] + score_sub,
                X[i][j],
                Y[i][j]
            )

            if M[i][j] == M[i-1][j-1] + score_sub:
                trace[i][j] = 'D'
            elif M[i][j] == X[i][j]:
                trace[i][j] = 'U'
            else:
                trace[i][j] = 'L'

    # Traceback
    align1, align2 = "", ""
    i, j = m, n

    while i > 0 or j > 0:
        direction = trace[i][j]
        if direction == 'D':
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            i -= 1
            j -= 1
        elif direction == 'U':
            align1 = s1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        elif direction == 'L':
            align1 = '-' + align1
            align2 = s2[j-1] + align2
            j -= 1
        else:
            break

    matches, align_len, identity = _calculate_identity(align1, align2)

    return {
        "score": M[m][n],
        "alignment1": align1,
        "alignment2": align2,
        "matches": matches,
        "alignment_length": align_len,
        "identity_percent": identity
    }


# =========================
# Smith–Waterman (Affine)
# =========================

def smith_waterman_affine(
    s1,
    s2,
    match=2,
    mismatch=-1,
    gap_open=-2,
    gap_extend=-1,
    sub_matrix=None
):
    """
    Local alignment with affine gap penalty.
    """

    m, n = len(s1), len(s2)

    M = [[0]*(n+1) for _ in range(m+1)]
    X = [[0]*(n+1) for _ in range(m+1)]
    Y = [[0]*(n+1) for _ in range(m+1)]

    trace = [[None]*(n+1) for _ in range(m+1)]

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m+1):
        for j in range(1, n+1):

            score_sub = _score(s1[i-1], s2[j-1], match, mismatch, sub_matrix)

            X[i][j] = max(
                M[i-1][j] + gap_open,
                X[i-1][j] + gap_extend
            )

            Y[i][j] = max(
                M[i][j-1] + gap_open,
                Y[i][j-1] + gap_extend
            )

            M[i][j] = max(
                0,
                M[i-1][j-1] + score_sub,
                X[i][j],
                Y[i][j]
            )

            if M[i][j] == 0:
                trace[i][j] = 'S'
            elif M[i][j] == M[i-1][j-1] + score_sub:
                trace[i][j] = 'D'
            elif M[i][j] == X[i][j]:
                trace[i][j] = 'U'
            else:
                trace[i][j] = 'L'

            if M[i][j] > max_score:
                max_score = M[i][j]
                max_pos = (i, j)

    # Proper local traceback (stop at 0)
    align1, align2 = "", ""
    i, j = max_pos

    while i > 0 and j > 0 and M[i][j] != 0:
        direction = trace[i][j]

        if direction == 'D':
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            i -= 1
            j -= 1
        elif direction == 'U':
            align1 = s1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        elif direction == 'L':
            align1 = '-' + align1
            align2 = s2[j-1] + align2
            j -= 1

    matches, align_len, identity = _calculate_identity(align1, align2)

    return {
        "score": max_score,
        "alignment1": align1,
        "alignment2": align2,
        "matches": matches,
        "alignment_length": align_len,
        "identity_percent": identity
    }
