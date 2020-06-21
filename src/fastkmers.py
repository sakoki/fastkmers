import os
import gzip
import argparse
import functools
import numpy as np
import multiprocessing as mp
from collections import defaultdict
from utils.decorator import timer


def make_N_kmers(sequence, N):
    """Count kmers of N size and their frequencies
    
    Parameters
    ----------
    N : int
        kmer size
        
    Returns
    -------
    defaultdict
        kmers and their frequencies 
    
    """
    kmers = defaultdict(int)    
    seq_length = len(sequence)
    
    for i in range(seq_length+1):
        if (N + i) <= seq_length:  # right index boundary
            kmers[sequence[i:N+i]] += 1 
        else:
            break
            
    return kmers


def merge_kmers(kmers):
    """Merge values of overlapping keys between two dictionaries
    
    Parameters
    ----------
    kmers : list of dict of {str: int}
    
    Returns
    -------
    dict of {str: int}
    
    """
    merged = {}
    for i in kmers:
        merged = {k: merged.get(k, 0) + i.get(k, 0) for k in merged.keys() | i.keys()}
    return merged


@timer
def get_sequences(file_path):
    """Extract sequences from a fastq file
    
    Parameters
    ----------
    file_path : str
        path to fastq file
        
    Returns
    -------
    list of nucleotide sequence fragments
        
    Notes
    -----
    May take considerably long for large fastq files
    
    """
    
    sequences = []

    with gzip.open(file_path, 'r') as fastq:
        for _ in fastq:
            try:
                sequences.append(next(fastq).decode("utf-8").rstrip())
                next(fastq)  # separator
                next(fastq)  # phred_score
            except StopIteration:
                print("Sequence with missing fields found. A sequence in a fastq file should contain 4 lines")

    return sequences


@timer
def quantify_kmers(sequences, N):
    """Generate kmers for multiple sequences, then merge results
    
    Parameters
    ----------
    sequences : array 
    N : int
    
    Returns
    -------
    kmers and their corresponding frequencies
    
    """
    return merge_kmers((make_N_kmers(sequence, N) for sequence in sequences))


def parallelize(array, func, n_cores, *args, **kwargs):
    """Split array amongst n_cores and parallelize opertation"""
    
    print(f'Parallelizing {func.__name__} on {n_cores} cores...')
    split_array = np.array_split(array, n_cores)
    pool = mp.Pool(n_cores)
    func = functools.partial(func, *args, **kwargs)
    results = pool.map(func, split_array)
    pool.close()
    pool.join()
    return results
            

def sort_results(kmers, descending=True):
    """Sort kmers by their frequencies
    
    Parameters
    ----------
    kmers : dict of {str : int}
    descending : bool
    
    
    Returns
    -------
    list of tuple of (str, int)
        Sorted kmers
    
    """
    
    return sorted(kmers.items(), key=lambda x: x[1], reverse=True)


def write_kmers(fname, kmers):
    """Save kmers to file"""
    with open(fname, 'w') as f:
        for i in kmers:
            f.write(f'{i[0]}\t{i[1]}\n')


def main(args):
    sequences = get_sequences(args.fastq)
    if args.cores > 1:
        if args.cores > os.cpu_count():
            raise ValueError(f"Too many cores requested: {args.cores} request cores > {os.cpu_count()} system cores")
        results = sort_results(merge_kmers(parallelize(array=sequences, 
                                                       func=quantify_kmers,
                                                       n_cores=args.cores,
                                                       N=args.N)))
    else:
        results = sort_results(quantify_kmers(args.fastq, args.N))
    write_kmers(args.output, results)

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq', type=str, nargs='+', help='Input fastq file (zipped)')
    parser.add_argument('--N', type=int, help='Size of kmers to quantify')
    parser.add_argument('--output', type=str, help='Output file name')
    parser.add_argument('--cores', type=int, help='Number of cpu cores to use for parallelization', default=1)
    args = parser.parse_args()
    main(args)
