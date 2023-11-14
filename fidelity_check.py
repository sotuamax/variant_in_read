import os 
import argparse
import pysam 
import re 
from Bio import SeqIO 
import pandas as pd 
import numpy as np 
import time
from mpi4py import MPI
from collections import Counter
import datetime
import sys 

VERSION = "1.0"

class ArgFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

def args_parser():
    '''parser the argument from terminal command'''
    header = f"fidelity_check.py: a multiple-processing program for the detection of variants in reads alignment to reference. "
    example = "Example: \nmpiexec -n CORES REF BAM -o OUT [-bed BED]"
    parser=argparse.ArgumentParser(epilog = example, formatter_class=lambda prog: ArgFormatter(prog,max_help_position=100,width=150), description="")
    parser.add_argument("--version", action = "version", version = f"software version, {VERSION}")
    parser.add_argument("ref", metavar = "REF", help = "reference fasta")
    parser.add_argument("bam", metavar = "BAM", help = "bam file")
    parser.add_argument("-bed", "--bed", metavar = "BED", required = False, help = "bed file used as template to count reads in regions. ")
    parser.add_argument("-o", "--output", metavar = "OUT", help = "output name", required = True)
    args = parser.parse_args()
    args.start_date = datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S") 
    args.version = "fidelity_check.v"+VERSION
    args.command=" ".join(sys.argv)
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

def parse_out(out):
    out_df = pd.read_table(f"{out}.txt", sep = "\t", header = 0)
    out_df["CIGAR"] = out_df["cigar"].str.replace("[0-9]", "", regex = True)
    out_df["temp"] = np.where(out_df["CIGAR"] == "D", out_df["cigar"].str.strip("D"), 0)
    out_df["temp"] = out_df["temp"].astype("int")
    out_df["end"] = out_df["pos"] + out_df["temp"] + 1
    out_freq = pd.DataFrame(out_df.groupby(by = ["ref", "pos", "end", "CIGAR"], as_index = False).size())
    # 
    ref_seq = sorted(set(out_freq["ref"]))
    cigar_set = sorted(set(out_freq["CIGAR"]))
    # 
    for c in cigar_set:
        c_list = list()
        for r in ref_seq:
            r_list = []
            sub_s = out_freq[out_freq["ref"] == r]
            for row in sub_s.itertuples():
                r_list += list(range(row.pos, row.end))
            r_df = pd.DataFrame.from_dict(dict(Counter(r_list)), orient = "index", columns = ["freq"])
            r_df["ref"] = r 
            r_df["pos"] = r_df.index
            c_list.append(r_df)
        C_df = pd.concat(c_list, axis = 0)
        C_df.to_csv(f"{out}_{c}_freq.txt", sep = "\t", index = False)

def update_freq(out, count_df):
    for c in ["D", "X", "I"]:
        freq_df = pd.read_table(f"{out}_{c}_freq.txt", sep = "\t", header = 0)
        freq_df_new = list()
        for row in count_df.itertuples():
            freq_df_sub = freq_df[(freq_df["ref"] == row.ref) & (freq_df["pos"] >= row.start) & (freq_df["pos"] < row.end)].copy()
            freq_df_sub["count"] = row.count
            freq_df_new.append(freq_df_sub)
        freq_df_new = pd.concat(freq_df_new, axis = 0)
        freq_df_new.to_csv(f"{out}_{c}_freq.txt", sep = "\t", index = False)

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
        print(f"Running {args.version}")
        print(f"  Start on: {args.start_date}")
        print(f"  Used command: {args.command}")
        print("==============================")
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
        count_all = comm.gather(count_df, root = 0)
    # 
    if rank == 0:
        mismatch_all_df = pd.concat(mismatch_all, axis = 0)
        indel_all_df = pd.concat(indel_all, axis = 0)
        all_df = pd.concat([mismatch_all_df, indel_all_df], axis = 0)
        all_df.to_csv(f"{out}.txt", sep = "\t", index = False)
        if args.bed != None:
            count_df = pd.concat(count_all, axis = 0)
            count_df.to_csv(f"{out}_count.txt", sep = "\t", index = False)
            parse_out(out)
            update_freq(out, count_df)
        test = 1
        endtime(start)
    else:
        test = None
    comm.bcast(test, root = 0)

if __name__ == "__main__":
    main()
