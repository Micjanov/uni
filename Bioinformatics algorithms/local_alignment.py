from Bio import SeqIO
from Bio.SubsMat import MatrixInfo

blosum = MatrixInfo.pam250
 
def align(seq1, seq2):
    # d is gap penalty
    d = -5
    # Initialize dynamic programming matrix
    S = [[0 for i in range(len(seq1)+1)] for j in range(len(seq2)+1)]
    
    # Initialization
    # set first row
    S[0] = [0 for i in range(len(seq1)+1)]
    # set first column
    for j,row in enumerate(S):
        row[0] = 0
    
    # Calculation
    # i is row, j is column
    best = 0
    bestpos = (0,0)

    for i in range(len(seq2)+1)[1:]:
        for j in range(len(seq1)+1)[1:]:

            f = [S[i-1][j-1] + score((seq2[i-1],seq1[j-1])), # match
                 S[i-1][j] + d,                              # gap
                 S[i][j-1] + d]                              # gap

            if max(f) < 0:
                S[i][j] = 0
            else:
                S[i][j] = max(f)
                
            if S[i][j] > best:
                best = S[i][j]
                bestpos = (i,j)
    pointer= [[0 if x < 0 else x for x in y] for y in S]

    # Traceback
    align1 = ""
    align2 = ""
    i, j = bestpos[0], bestpos[1]
    while i > 0 and j > 0 and pointer[i][j] > 0:
        score_current = S[i][j]
        score_diagonal = S[i-1][j-1]
        score_up = S[i][j-1]
        score_left = S[i-1][j]
        
        if score_current == score_diagonal + score((seq1[j-1], seq2[i-1])):
            align1 += seq1[j-1]
            align2 += seq2[i-1]
            i -= 1
            j -= 1
        elif score_current == score_up + d:
            align1 += seq1[j-1]
            align2 += '-'
            j -= 1
        elif score_current == score_left + d:
            align1 += '-'
            align2 += seq2[i-1]
            i -= 1

    align1 = align1[::-1]
    align2 = align2[::-1]

    return best, align1, align2
    
def score(pair):

    if pair in blosum:
        return blosum[pair]
    else:
        return blosum[tuple(reversed(pair))]


def read_sequences(filename):

    fa = list(SeqIO.parse(filename, "fasta"))
    seq1 = str(fa[0].seq)
    seq2 = str(fa[1].seq)

    return seq1, seq2

def local_alignment_score(file):
    
    seq1 , seq2 = read_sequences(file)
    return align(seq1, seq2)[0]

def local_alignment(file):
    
    seq1 , seq2 = read_sequences(file)
    align1 , align2 = align(seq1, seq2)[1:]
    return (align1 , align2)