import jellyfish
import numpy as np
import Bio.SeqIO


def count_kmers(seq, k, data):
    # Extracts k-mers for a given sequence seq and saves it in a dictionary data with associated frequency
    size = len(seq)
    for i in range(size - k + 1):
        kmer = seq[i: i + k]
        try:
            data[kmer] += 1
        except KeyError:
            data[kmer] = 1


def count_all(gen, k):
    # Extracts k-mers with associated frequency for a list of sequence gem
    data = {}
    print(f'Counting {k}mers ...')
    for i in range(len(gen)):
        count_kmers(gen[i].seq, k, data)
    print(f'Process completed: {len(data)} {k}-mers found')
    return data

def clear_dict(dic,f):
    # Clear the dictionaries eliminating k-mers with less than f_min occurrencies
    return {key:val for key, val in dic.items() if val >= f}

def delete_copies(dict1, dict2):
    # Deletes element that are common in both dictionary
    lst = []
    print(f'Deleting common k-mers ...')
    for key in dict1:
        if key in dict2.keys():
            lst.append(key)
    for key in lst:
        del dict1[key]
        del dict2[key]


def concatenate(dic, k):
    # Concatenates k-mers (explained in report)
    concatenated = []
    to_delete = set()
    stop = False
    print(f'Concatenating {k}mers ...')
    for key1 in dic:
        string = str(key1)
        if key1 not in to_delete:
            for key2 in dic:
                if key2 != key1 and key2 not in to_delete:
                    last = string[-k + 1:]
                    first = string[0:k - 1]
                    string2 = str(key2)
                    last2 = string2[-k + 1:]
                    first2 = string2[0:k - 1]
                    if last == first2:
                        string += string2[-1]
                        to_delete.add(key2)
                        to_delete.add(key1)
                    elif first == last2:
                        string = string2[0] + string
                        to_delete.add(key2)
                        to_delete.add(key1)
            concatenated.append(string)
            print(f'Process completed')
    return concatenated


def calc_dist(seq1, seq2):
    # Computes Levenshtein distance between sequences in the different list seq1 and seq2
    print(f'Computing distance between sequences ...')
    matrix = np.zeros((len(seq1), len(seq2)))
    i = 0
    for read1 in seq1:
        j = 0
        for read2 in seq2:
            matrix[i, j] = jellyfish.levenshtein_distance(str(read1), str(read2))
            j += 1
        i += 1
    return matrix

def print_results(ind,seq1,seq2):
    l=len(ind[0])
    i=0
    while i<l:
        index=(ind[0][i],ind[1][i])
        print(seq1(ind[0][i])
        print(seq2(ind[1][i])
        i+=1
    print(i[0],a[i[0]])
    print(i[1],b[i[1]])
    
