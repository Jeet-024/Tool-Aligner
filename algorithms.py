#Needleman-Wunsch and Smith-Waterman algorithm
def _init_matrix(rows, cols, fill=0):
    return [[fill] * cols for _ in range(rows)]


def needleman_wunsch(s1, s2, match=1, mismatch=-1, gap=-1):
    n, m = len(s1), len(s2)
    dp = _init_matrix(n + 1, m + 1)
    bt = _init_matrix(n + 1, m + 1)

    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + gap
        bt[i][0] = 'U'
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + gap
        bt[0][j] = 'L'

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score_diag = dp[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            score_up = dp[i - 1][j] + gap
            score_left = dp[i][j - 1] + gap
            best = max(score_diag, score_up, score_left)
            dp[i][j] = best
            if best == score_diag:
                bt[i][j] = 'D'
            elif best == score_up:
                bt[i][j] = 'U'
            else:
                bt[i][j] = 'L'

    # Traceback
    i, j = n, m
    a1, a2 = [], []
    matches = 0
    while i > 0 or j > 0:
        move = bt[i][j]
        if move == 'D':
            a1.append(s1[i - 1])
            a2.append(s2[j - 1])
            if s1[i - 1] == s2[j - 1]:
                matches += 1
            i -= 1
            j -= 1
        elif move == 'U':
            a1.append(s1[i - 1])
            a2.append('-')
            i -= 1
        else:
            a1.append('-')
            a2.append(s2[j - 1])
            j -= 1

    aln1 = ''.join(reversed(a1))
    aln2 = ''.join(reversed(a2))
    aln_len = len(aln1)
    identity = matches / aln_len * 100 if aln_len > 0 else 0.0

    return {
        'score': dp[n][m],
        'alignment1': aln1,
        'alignment2': aln2,
        'matches': matches,
        'alignment_length': aln_len,
        'identity_percent': identity,
    }


def smith_waterman(s1, s2, match=2, mismatch=-1, gap=-1):
    n, m = len(s1), len(s2)
    dp = _init_matrix(n + 1, m + 1, 0)
    bt = _init_matrix(n + 1, m + 1, None)

    max_i = max_j = 0
    max_score = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score_diag = dp[i - 1][j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            score_up = dp[i - 1][j] + gap
            score_left = dp[i][j - 1] + gap
            best = max(0, score_diag, score_up, score_left)
            dp[i][j] = best
            if best == 0:
                bt[i][j] = None
            elif best == score_diag:
                bt[i][j] = 'D'
            elif best == score_up:
                bt[i][j] = 'U'
            else:
                bt[i][j] = 'L'

            if best > max_score:
                max_score = best
                max_i, max_j = i, j

    # Traceback from max
    i, j = max_i, max_j
    a1, a2 = [], []
    matches = 0
    while i > 0 and j > 0 and bt[i][j] is not None:
        move = bt[i][j]
        if move == 'D':
            a1.append(s1[i - 1])
            a2.append(s2[j - 1])
            if s1[i - 1] == s2[j - 1]:
                matches += 1
            i -= 1
            j -= 1
        elif move == 'U':
            a1.append(s1[i - 1])
            a2.append('-')
            i -= 1
        else:
            a1.append('-')
            a2.append(s2[j - 1])
            j -= 1

    aln1 = ''.join(reversed(a1))
    aln2 = ''.join(reversed(a2))
    aln_len = len(aln1)
    identity = matches / aln_len * 100 if aln_len > 0 else 0.0

    return {
        'score': max_score,
        'alignment1': aln1,
        'alignment2': aln2,
        'matches': matches,
        'alignment_length': aln_len,
        'identity_percent': identity,
    }

