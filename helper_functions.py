import numpy as np
from Bio.Align import substitution_matrices

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """

    m, n = len(seq1), len(seq2)

    # Initialise DP table (score matrix) and Traceback matrix
    dp = np.zeros((m + 1, n + 1))
    traceback = np.zeros((m + 1, n + 1), dtype=int)

    # Calculate gap penatlies, assuming uniform linear gap penalty from scorin function('-', 'A')
    gap_penalty_seq1 = scoring_function("-", seq2[0]) if n > 0 else -1
    gap_penalty_seq2 = scoring_function(seq1[0], "-") if m > 0 else -1


    # Fill base cases (edges) 
    for i in range(1, m + 1):
        dp[i][0] = dp[i-1][0] + gap_penalty_seq2
        traceback[i][0] = 2
    
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j - 1] + gap_penalty_seq1
        traceback[0][j] = 3 

    for i in range(1, m + 1):
        char1 = seq1[i - 1]
        for j in range(1, n + 1):
            char2 = seq2[j - 1]

            match_score = dp[i - 1][j - 1] + scoring_function(char1, char2)
            delete_score = dp[i - 1][j] + scoring_function(char1, "-")  # Up
            insert_score = dp[i][j - 1] + scoring_function("-", char2)  # Left

            # Max score selection
            max_score = max(match_score, delete_score, insert_score)
            dp[i][j] = max_score

            # Record path for traceback
            if max_score == match_score:
                traceback[i][j] = 1  # Diagonal
            elif max_score == delete_score:
                traceback[i][j] = 2  # Up
            else:
                traceback[i][j] = 3  # Left

    # 3. Traceback from bottom-right (m, n) to top-left (0, 0)
    aligned_seq1 = []
    aligned_seq2 = []

    i, j = m, n
    while i > 0 or j > 0:
        if traceback[i][j] == 1:  # Diagonal (Match/Mismatch)
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif traceback[i][j] == 2:  # Up (Gap in seq2)
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        elif traceback[i][j] == 3:  # Left (Gap in seq1)
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    # Reverse sequences to restore original order
    aligned_seq1 = "".join(reversed(aligned_seq1))
    aligned_seq2 = "".join(reversed(aligned_seq2))
    final_score = float(dp[m][n])

    return aligned_seq1, aligned_seq2, final_score



def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    m, n = len(seq1), len(seq2)
    
    # 1. Initialise scoring and tracking matrices
    H = [[0.0] * (n + 1) for _ in range(m + 1)]
    trace = [["Stop"] * (n + 1) for _ in range(m + 1)]
    
    max_score = 0.0
    max_i, max_j = 0, 0
    
    # 2. Fill matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Pass characters dynamically into the scoring_function argument
            diagonal = H[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
            up       = H[i-1][j]   + scoring_function(seq1[i-1], "-")
            left     = H[i][j-1]   + scoring_function("-", seq2[j-1])
            
            # Smith-Waterman rule: clip negative scores to 0.0
            H[i][j] = max(0.0, diagonal, up, left)
            
            # Map tracking for traceback path
            if H[i][j] > 0.0:
                if H[i][j] == diagonal:
                    trace[i][j] = "D"
                elif H[i][j] == up:
                    trace[i][j] = "U"
                elif H[i][j] == left:
                    trace[i][j] = "L"
            
            # Track max score position for matrix
            if H[i][j] > max_score:
                max_score = H[i][j]
                max_i, max_j = i, j

    # 3. Traceback
    align1, align2 = [], []
    i, j = max_i, max_j
    
    # Stop when score is 0
    while i > 0 and j > 0 and H[i][j] > 0.0 and trace[i][j] != "Stop":
        if trace[i][j] == "D":
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif trace[i][j] == "U":
            align1.append(seq1[i-1])
            align2.append("-")
            i -= 1
        elif trace[i][j] == "L":
            align1.append("-")
            align2.append(seq2[j-1])
            j -= 1

    return "".join(reversed(align1)), "".join(reversed(align2)), float(max_score)




## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)

# Load BLOSUM62 matrix
blosum62 = substitution_matrices.load("BLOSUM62")
GAP_PENALTY = -4  


def blosum62_scoring(aa1, aa2, gap_penalty=GAP_PENALTY):
    if aa1 == "-" or aa2 == "-":
        return gap_penalty
    try:
        return blosum62[aa1, aa2]
    except KeyError:
        return blosum62.get((aa1, aa2), -4)
    

# Main to test the examples given in the stubs
if __name__ == "__main__":
    import doctest
    print("Running doctest evaluation...")
    result = doctest.testmod()
    if result.failed == 0:
        print(" Success! Your code perfectly matches the assignment example requirements.")
    else:
        print(f" Fail: {result.failed} mismatches detected.")
