#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 01:10:24 2020

@author: jasonmoggridge


DEFINITIONS

There are many assembly tools, but none of them is perfect. 
Biologists therefore need to evaluate the quality of various 
assemblers by comparing their results. In our case, once we
 have run the SPAdes assembler on a set of reads, we need to 
 test the quality of the resulting assembly.

Contig: A contiguous segment of the genome that has been
 reconstructed by an assembly algorithm

Scaffold: An ordered sequence of contigs (possibly separated 
by gaps between them) that are reconstructed by an assembly 
algorithm. The order of contigs in a correctly assembled scaffold 
corresponds to their order in the genome. Existing assemblers 
specify the approximate lengths of gaps between contigs in a 
scaffold.

                                          

N50 statistic: N50 is a statistic that is used to measure 
the quality of an assembly. N50 is defined as the maximal 
contig length for which all contigs greater than or equal
 to that length comprise at least half of the sum of the 
 lengths of all the contigs. For example, consider the five
 toy contigs with the following lengths: [10, 20, 30, 60, 70]
 . Here, the total length of contigs is 190, and contigs of 
 length 60 and 70 account for at least 50% of the total length
 of contigs (60 + 70 = 130), but the contig of length 70 does
 not account for 50% of the total length of contigs. Thus, N50
 is equal to 60.

NG50 statistic: The NG50 length is a modified version of N50
 that is defined when the length of the genome is known (or 
 can be estimated). It is defined as the maximal contig length
 for which all contigs of at least that length comprise at 
 least half of the length of the genome. NG50 allows for 
 meaningful comparisons between different assemblies for the
 same genome. For example, consider the five toy contigs we 
 considered previously: [10, 20, 30, 60, 70]. These contigs 
 only add to 190 nucleotides, but say that we know that the 
 genome from which they have been generated has length 300. 
 In this example, the contigs of length 30, 60, and 70 account
 for at least 50% of the genome length (30 + 60 + 70 = 160);
 but the contigs of length 60 and 70 no longer account for 
 at least 50% of the genome length (60 + 70 = 130). Thus, 
 NG50 is equal to 30.

NGA50 statistic: If we already know a reference genome 
for a species, then we can test the accuracy of a newly 
assembled genome against this reference. The NGA50 statistic
 is a modified version of NG50 accounting for assembly errors
 (called misassemblies). To compute NGA50, errors in the 
 contigs are accounted for by comparing contigs to a reference
 genome. All of the misassembled contigs are broken at 
 misassembly breakpoints, resulting in a larger number of 
 contigs with the same total length. For example, if there
 is a missasembly breakpoint at position 10 in a contig of
 length 30, this contig will be broken into contigs of length
 10 and 20.

NGA50 is calculated as the NG50 statistic for the set of
 contigs resulting after breaking at misassembly breakpoints.
 For example, consider our example before, for which the genome
 length is 300. If the largest contig in [10, 20, 30, 60, 70]
 is broken into two contigs of length 20 and 50 (resulting in
 the set of contigs [10, 20, 20, 30, 50, 60]), then. contigs
 of length 20, 30, 50, and 60 account for at least 50% of the
 genome length (20 + 30 + 50 + 60 = 160). But contigs of
 length 30, 50, and 60 do not account for at least 50% of the
 genome length (30 + 50 + 60 = 140). Thus, NGA50 is equal to 20.

"""


def NX_statistic(contig_lengths, X): # must be sorted list
    contig_lengths = list(reversed(sorted(contig_lengths)))
    cumulative_length = []
    for i in range(len(contig_lengths)):
        cumulative_length.append(sum(contig_lengths[:i+1]))
    for i in range(len(cumulative_length)):
        if cumulative_length[i] >= (X/100)*cumulative_length[-1]:
            NX = contig_lengths[i]
            LX = i+1
            break
    return (NX, LX)


def NGX_statistic(contig_lengths, genome_size, X):

    cumulative_length = []
    for i in range(len(contig_lengths)):
        cumulative_length.append(sum(contig_lengths[:i+1]))
    for i in range(len(cumulative_length)):
        if cumulative_length[i] >= (X/100)*genome_size:
            NGX = contig_lengths[i]
            LGX = i+1
            break
    return (NGX, LGX)


contig_lengths = [200,1,100, 80, 60, 60, 50, 50, 30, 30, 20, 20]
staph_genome = 1000 
Nstats = NX_statistic(contig_lengths, 50)
NGstats = NGX_statistic(contig_lengths, staph_genome, 50)
    
    