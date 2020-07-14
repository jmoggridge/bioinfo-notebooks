"""
BA11-H ->
Spectral dictionaries...
https://stepik.org/lesson/11866/step/9?unit=8305

"""




def Size_of_spectral_dictionary(spectrum, min_s, max_s):

    """ Recursive function to count all possible peptides 
    scoring between min -> max score (*note -> never let t<0 in path)"""
    
    
    def Size(i,t):
         
        # do not count if score dips below zero (silly but yeah)
        if i < 0 or t < 0:
            return 0
        
        # recursion -> get sum of peptides leading to (i,t)
        #  size(i,t) =  sum -> size(i-|aa|,t-s[i]) -> for all AAs 

        elif (i,t) not in size_memo.keys():
            size_memo[(i,t)] = sum(Size(i - aa, t - spectrum[i]) for aa in AAs)

        return size_memo[(i,t)]
    
    ## wrapper
    
    # initialize size_memo: an (i*t) matrix as dictionary
    # size: sum of peptides in spectral_dict   
    size_memo = {(0,0):1}
    sizes = []
    mass = len(spectrum)-1
    
    # compute size(m,score) for all scores t in range(acceptable scores)    
    for t in range(min_s, max_s + 1):
        sizes.append(Size(mass, t))
    return sizes

# Main --

# list of AA masses:
    
AAs = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115,\
       128, 128, 129, 131, 137, 147, 156, 163, 186]
#AAs = [4,5] # AAs for sample

# Parsing
with open('data/rosalind_ba11h.txt', 'r') as infile:
    
    spectrum = [0] + [int(i) for i in infile.readline().strip().split()] 
    t = int(infile.readline())
    max_score = int(infile.readline())

# Get array of sizes of spectral dicts for each score in range t->max
size_dict = Size_of_spectral_dictionary(spectrum, t, max_score)
print('>Size of spectral dict = ',sum(size_dict))


