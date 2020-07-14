

"""Wrapper function for Gibbs sampler motif searcher"""    
    
def Gibbs_Sampler(Dna, k, t, N):        
    
    import random
    
    
    def random_kmers(Seqs,k):
        motifs = [] 
        for seq in Seqs:
            x = random.randint(0, len(seq)-k)
            motifs.append(seq[x:x+k])
        return motifs
    
    
    def Build_profile_matrix(motifs, k):
        profile = [[1/(len(motifs)+4) for nt in range(len(alpha))] for position in range(k)]
        for p in range(k):
            for motif in motifs:
                profile[p][alpha.index(motif[p])] += 1/(len(motifs)+4)    
        return profile
    
    
    def profile_rand_motif(seq,k,profile):
        
        kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
        c_weights = []
        c=0
    
        for kmer in kmers:
            prob = 1
            for p in range(k):
                prob = prob * profile[p][alpha.index(kmer[p])] 
            c+=prob
            c_weights.append(c)
    
        j = random.uniform(0,c_weights[-1])
        c_weights.append(float('inf'))
        for i in range(len(c_weights)):
            if c_weights[i] >= j:   
                return kmers[i]
    
        
    def Hamming(motifs):
    
        consensus = [[0 for nt in range(4)] for position in range(k)]
        for p in range(k):
            for motif in motifs:
                consensus[p][alpha.index(motif[p])] += 1
                
        mismatches = 0
        for counts in consensus:
            for c in range(4): 
                if c != counts.index(max(counts)):
                    mismatches += counts[c]
        return mismatches
        ##
    
    alpha = 'ACGT'
    Best = float('inf')
    Best_motifs = []
    Motifs = random_kmers(Dna,k)

    for _ in range(N):

        i = random.randint(0,len(Dna)-1)
        profile_these_motifs = Motifs[:i] + Motifs[i+1:]
        Profile = Build_profile_matrix(profile_these_motifs, k)    
        Motifs[i] = profile_rand_motif(Dna[i],k, Profile)    
        score = Hamming(Motifs)

        if score < Best:
            Best_motifs = Motifs
            Best = score

    return Best_motifs, Best
    
###############################################################################

    
#alpha ='ACGT'       
#if __name__ == "__main__":
#    k, t, N = [int(a) for a in input().strip().split(" ")]
#    Dna = []
#    for _ in range(t):
#        Dna.append(input())
 
    
f = open('//Users/jasonmoggridge/Desktop/rosalind_test.txt', 'r')
k, t, N = (int(i) for i in f.readline().strip().split(' '))
Dna = list(str(l.strip('\n')) for l in f.readlines())

best_score = float('inf')
best_motifs = []
for _ in range(1000):
    ans, score = Gibbs_Sampler(Dna, k, t, N)
    if score < best_score:
        best_score = score
        best_motifs = ans
    
for a in ans:
    print(a)

