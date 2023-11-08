import os 
import argparse
import pysam 
import re 
from Bio import SeqIO 
import pandas as pd 
import numpy as np 
import time
from mpi4py import MPI


def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("ref", help = "reference fasta")
    parser.add_argument("bam", help = "bam file")
    parser.add_argument("-bed", "--bed", required = False, help = "bed file used as template to count reads in regions. ")
    parser.add_argument("-o", "--output", help = "output name")
    args = parser.parse_args()
    return args

def seq_dict(args):
    seq = SeqIO.parse(args.ref, "fasta")
    ref_dict = SeqIO.to_dict(seq)
    return ref_dict

def pileup_mismatch(args, ref_dict, sub_loc):
    """use pysam.pileup to walk through each bp position and check mismatches in read column compared to reference column."""
    ### bam handle
    bam_handle = pysam.AlignmentFile(args.bam, "rb")
    ### for loop each sequence in ref_dict
    mismatch_list = list()
    for row in sub_loc.itertuples():
        rref, rstart, rend = row.ref, row.start, row.end
        for column in bam_handle.pileup(rref, rstart, rend, truncate = True):
            ref_pos = column.reference_pos
            ref_bp = ref_dict[rref].seq[ref_pos]
            for read in column.pileups:
                if read.query_position != None and read.alignment.query_sequence[read.query_position] != ref_bp:
                    mismatch_list.append((read.alignment.query_name, "X", row.ref, ref_pos, ref_bp, read.alignment.query_sequence[read.query_position]))
    mismatch_df = pd.DataFrame(mismatch_list, columns = ["read", "cigar", "ref", "pos", "ref_bp", "read_bp"])
    return mismatch_df

def fetch_indel(args, ref_dict, sub_loc):
    """"use pysam.fetch to get all insertion and deletion in the alignment file. """
    bam_handle = pysam.AlignmentFile(args.bam, "rb")
    indel_list = list()
    for row in sub_loc.itertuples():
        rref, rstart, rend = row.ref, row.start, row.end
        for read in bam_handle.fetch(rref, rstart, rend):
            read_indel = re.findall("[0-9]+[I|D]", read.cigarstring) # get I/D cigar
            aln_block = read.get_blocks() # alignment blocks in reference coordinates
            read_block = np.array(aln_block) - read.reference_start + read.qstart # alignmment blocks in read (initial version not corrected for insertion/deletion)
            ## update read_block
            if len(aln_block) > 1:
                for rr in range(0, len(read_indel)):
                    rname = read.query_name
                    cigar_tag = read_indel[rr]
                    ref_pos = aln_block[rr][-1]
                    ref_bp = str(ref_dict[rref].seq[ref_pos-1:aln_block[rr+1][0]])
                    if cigar_tag.endswith("I"): # insertion
                        read_block[rr+1:] = read_block[rr+1:] + int(cigar_tag.strip("I"))
                    if cigar_tag.endswith("D"): # deletion
                        read_block[rr+1:] = read_block[rr+1:] - int(cigar_tag.strip("D"))
                    ### get read_bp based on updated read_block
                    read_bp = read.query_sequence[read_block[rr][-1]-1:read_block[rr+1][0]]
                    indel_list.append((rname, cigar_tag, rref, ref_pos, ref_bp, read_bp))
    indel_df = pd.DataFrame(indel_list, columns = ["read", "cigar", "ref", "pos", "ref_bp", "read_bp"])
    return indel_df

def bam_count(args, sub_loc):
    """count reads in each of reference sequences. """
    bam_handle = pysam.AlignmentFile(args.bam, "rb")
    count_list = list()
    for row in sub_loc.itertuples():
        region_read = {read.query_name for read in bam_handle.fetch(row.ref, row.start, row.end)}
        count_list.append(len(region_read))
    return count_list

def endtime(start):
    """elapsed time"""
    end = time.time()
    t = end-start
    if t < 60:
        print('{:.2f} seconds elapsed'.format(t))
    elif t < 3600:
        print('{:.2f} minutes elapsed'.format(t/60))
    else:
        print('{:.2f} hours elapsed'.format(t/3600))
    return end

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    args = args_parser()
    out = args.output
    ref_dict = seq_dict(args)
    if rank == 0:
        start = time.time()
        split_task = {w:[] for w in range(size)}
        t_idx = 0
        bed = args.bed
        if bed != None:
            loc_df = pd.read_table(bed, sep = "\t", header = None, names = ["ref", "start", "end"])
            for row in loc_df.itertuples():
                split_task[t_idx].append(row.Index)
                t_idx = (t_idx+1)%size
        else:
            loc_df = pd.DataFrame([(r, 0, len(ref_dict[r].seq)) for r in ref_dict], columns = ["ref", "start", "end"])
            for row in loc_df.itertuples():
                split_task[t_idx].append(row.Index)
                t_idx = (t_idx+1)%size
    else:
        split_task = None
        loc_df = None
    loc_df = comm.bcast(loc_df, root = 0)
    split_task = comm.bcast(split_task, root = 0)
    ### 
    if rank == 0:
        print("Mismatches ......")
    sub_loc = loc_df[loc_df.index.isin(split_task[rank])]
    mismatch_rank = pileup_mismatch(args, ref_dict, sub_loc)
    mismatch_all = comm.gather(mismatch_rank, root = 0)
    # 
    if rank == 0:
        print("INDELs ......")
    indel_rank = fetch_indel(args, ref_dict, sub_loc)
    indel_all = comm.gather(indel_rank, root = 0)
    # 
    if args.bed != None:
        if rank == 0:
            print("Count ...... ")
        count_df = sub_loc.copy()
        count_df["count"] = bam_count(args, sub_loc)
        count_all = comm.gather(sub_loc, root = 0)
    # 
    if rank == 0:
        pd.concat(mismatch_all, axis = 0).to_csv(f"{out}_mismatch.txt", sep = "\t", index = False)
        pd.concat(indel_all, axis = 0).to_csv(f"{out}_indel.txt", sep = "\t", index = False)
        if args.bed != None:
            pd.concat(count_all, axis = 0).to_csv(f"{out}_count.txt", sep = "\t", index = False)
        test = 1
        endtime(start)
    else:
        test = None
    comm.bcast(test, root = 0)

if __name__ == "__main__":
    main()
