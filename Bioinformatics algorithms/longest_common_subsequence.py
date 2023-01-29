#Write a function lcs that takes two arguments:
#i) the location of a FASTA file containing two DNA strings  and  and ii)
# another file location. The function must write a longest common subsequence 
#of  and  in FASTA format to the file whose location is passed as the
# second argument.
from Bio import SeqIO
def lcs(input_fasta,output_fasta):
    lijst_seq=[]
    for seq_reqord in SeqIO.parse(input_fasta, "fasta"):
        lijst_seq.append(str(seq_reqord.seq))
    
    S1=lijst_seq[0]
    S2=lijst_seq[1]
    m=len(S1)
    n=len(S2)
    
    L = [[0 for x in range(n+1)] for x in range(m+1)]

    # Building the mtrix in bottom-up way
    for i in range(m+1):
        for j in range(n+1):
            if i == 0 or j == 0:
                L[i][j] = 0
            elif S1[i-1] == S2[j-1]:
                L[i][j] = L[i-1][j-1] + 1
            else:
                L[i][j] = max(L[i-1][j], L[i][j-1])

    index = L[m][n]

    lcs_algo = [""] * (index+1)
    lcs_algo[index] = ""

    i = m
    j = n
    while i > 0 and j > 0:

        if S1[i-1] == S2[j-1]:
            lcs_algo[index-1] = S1[i-1]
            i -= 1
            j -= 1
            index -= 1

        elif L[i-1][j] > L[i][j-1]:
            i -= 1
        else:
            j -= 1
    
    begin_output='>seq01\n'
    einde_output=""
    for letter in lcs_algo:
        einde_output+=str(letter)
    string=begin_output+einde_output
    output = open(output_fasta, 'w')
    output.write(string)
    output.close()
    # Printing the sub sequences
    #return "S1 : " + S1 + "\nS2 : " + S2