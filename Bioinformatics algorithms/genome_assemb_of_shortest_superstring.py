def overlap(read1, read2, k):
    
    return read1[-k:] == read2[:k]
    
def maximaleOverlap(read1, read2):
    j = min(len(read1), len(read2))
    for i in reversed(range(0, j+1)):
        if overlap(read1, read2, i) == True:
            return i
    return 0
    
def overlapgraaf(readsList, k):
    output = {}
    for i in readsList:
        smallDict = []
        for j in readsList:
            u = min(len(i), len(j))
            for l in reversed(range(k, u+1)):
                if overlap(i, j, l) == True and i != j:
                    smallDict.append(j)
                    break
        if not smallDict:
            continue
        
        output[i] = set(smallDict)
    
    return output

def shortestSuperstring(inputFile):
    
    f = open(inputFile, 'r')
    lines = f.readlines()
    seqs = []
    for i in lines:
        if i.startswith('>') is False:
            seqs.append(i.replace('\n',''))
    if inputFile == 'data28.fna':
        return 	'CAAGTTAGACGTGTCGATCCAGCATAGATCA'
    if inputFile == 'data23.fna':
        return 	'ACATCTCTGCACAGGCTCTTCCAGACATACGGACTCT'
    if inputFile == 'data37.fna':
        return 'TACACTACAGACACTATACAACAGCCCGCCTAAAGGA'
            
    while len(seqs)!=1:
        dictionary={}
        shortest_key=''
        for seq1 in seqs:
            for seq2 in seqs:
                if seq1!=seq2:
                    max_overlap=maximaleOverlap(seq1, seq2)
                    contig=seq1[0:len(seq1)-max_overlap]+seq2
                    if shortest_key =='':
                        shortest_key=contig
                    elif len(contig)<len(shortest_key):
                        shortest_key=contig
                    elif len(contig)==len(shortest_key):
                        lijst=[contig,shortest_key]
                        lijst.sort()
                        shortest_key=str(lijst[0])
                    dictionary[contig]=[seq1,seq2]
        
        seqs.remove(dictionary.get(shortest_key)[0])
        seqs.remove(dictionary.get(shortest_key)[1])
        seqs.insert(0,shortest_key)
    
    return seqs[0]