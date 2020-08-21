#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
from wiggletools.wiggle import Wiggle
from wiggletools.wiggle_matrix import WiggleMatrix
import os
import glob
from Bio import SeqIO
import multiprocessing as mp
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refseq_files", required=True, help="", type=str, nargs="+")
    parser.add_argument("--wiggle_files", required=True, help="", type=str, nargs="+")
    parser.add_argument("--processes", default=1, help="", type=int)
    parser.add_argument("--output_dir", default=None, help="", type=str)
    subparsers = parser.add_subparsers(help="commands", dest="command")

    step_height_parser = subparsers.add_parser("to_step_height")
    step_height_parser.add_argument("--step_range", default=3, help="", type=int)
    step_height_parser.add_argument("--step_direction", choices=["start_end", "end_start"],
                                    required=True, help="", type=str)
    step_height_parser.add_argument("--output_prefix", default="STEP_HIEGHT_", help="", type=str)

    percentile_parser = subparsers.add_parser("to_percentile")
    percentile_parser.add_argument("--nth", required=True, type=float, choices=range(0, 101, 1))
    percentile_parser.add_argument("--scope", required=True, type=str, choices=["global", "stretch"])
    percentile_parser.add_argument("--output_prefix", default=f"nth_PERCENTILE_", help="", type=str)

    log2_parser = subparsers.add_parser("to_log2")
    log2_parser.add_argument("--output_prefix", default="LOG2_", help="", type=str)

    log10_parser = subparsers.add_parser("to_log10")
    log10_parser.add_argument("--output_prefix", default="LOG10_", help="", type=str)

    combine_parser = subparsers.add_parser("agg_merge")
    combine_parser.add_argument("--by", choices=["max", "min", "average"], help="", type=str)
    combine_parser.add_argument("--output_prefix", default="by_AGG_", help="", type=str)

    split_parser = subparsers.add_parser("split")
    split_parser.add_argument("--by", choices=["seqid", "fasta"], help="", type=str)

    arth_parser = subparsers.add_parser("arithmetic")
    arth_parser.add_argument("--opt", choices=["div", "mul", "sub", "add"], help="", type=str)
    arth_parser.add_argument("--value", required=True, help="", type=float)
    arth_parser.add_argument("--output_prefix", default="", help="", type=str)

    args = parser.parse_args()
    # Handling main arguments
    fasta_pathes = []

    for item in args.refseq_files:
        for sub_item in glob.glob(item):
            fasta_pathes.append(os.path.abspath(sub_item))

    wiggle_pathes = []
    for item in args.wiggle_files:
        for sub_item in glob.glob(item):
            wiggle_pathes.append(os.path.abspath(sub_item))

    chrom_sizes = get_chrom_sizes(fasta_pathes)
    pool = mp.Pool(processes=args.processes)
    processes = []
    if len(wiggle_pathes) == 0:
        print("Error: No wiggles loaded")
        exit(1)
    for wiggle_path in wiggle_pathes:
        processes.append(pool.apply_async(call_functions, (os.path.abspath(wiggle_path), chrom_sizes, args, )))
    if args.command == "agg_merge":
        parsed_wiggles = [p.get() for p in processes]
        wig_matrix = WiggleMatrix(parsed_wiggles, chrom_sizes, processes=args.processes)
        wig_matrix.agg_merge(by=args.by)
        wig_matrix.write_matrix_to_wiggle_files(
            wig_matrix.f_wiggle_matrix_df.loc[:, ["seqid", "location", "agg_col_forward"]],
            os.path.abspath(args.output_dir),
            f'{args.output_prefix.replace("by", args.by.upper())}forward')
        wig_matrix.write_matrix_to_wiggle_files(
            wig_matrix.r_wiggle_matrix_df.loc[:, ["seqid", "location", "agg_col_reverse"]],
            os.path.abspath(args.output_dir),
            f'{args.output_prefix.replace("by", args.by.upper())}reverse')
    else:
        for p in processes:
            p.get()
    pool.close()
    exit(0)


def call_functions(wig_path, chrom_sizes, args):
    wiggle = Wiggle(wig_path, chrom_sizes)
    writer_flag = True
    if args.command == "to_log2":
        wiggle.to_log2(inplace=True)
        writer_flag = True
    elif args.command == "to_log10":
        wiggle.to_log10(inplace=True)
        writer_flag = True
    elif args.command == "to_percentile":
        wiggle.to_percentile(args.nth, scope=args.scope, inplace=True)
        writer_flag = True
        args.output_prefix = args.output_prefix.replace("nth", f"{args.nth}th")
    elif args.command == "to_step_height":
        wiggle.to_step_height(args.step_range, args.step_direction, inplace=True)
        writer_flag = True
    elif args.command == "agg_merge":
        return wiggle.get_wiggle(is_full=False)
    elif args.command == "split":
        wiggle.split_wiggle(args.by, args.output_dir)
        writer_flag = False
    elif args.command == "arithmetic":
        wiggle.arithmethic(args.opt, args.value, inplace=True)
        writer_flag = True
    else:
        print(f"Error: unexpected command {args.command}.")
        exit(1)
    if writer_flag:
        if args.output_dir is None:
            wiggle.write_wiggle(
                f"{os.path.dirname(wiggle.file_path)}/{args.output_prefix}"
                f"{os.path.basename(wiggle.file_path)}")
        else:
            wiggle.write_wiggle(
                f"{os.path.abspath(args.output_dir)}/{args.output_prefix}{os.path.basename(wiggle.file_path)}")


def get_chrom_sizes(fasta_pathes):
    ret_list = []
    for fasta_path in fasta_pathes:
        print(f"==> Parsing reference sequence: {os.path.basename(fasta_path)}")
        for seq_record in SeqIO.parse(fasta_path, "fasta"):
            ret_list.append({"seqid": seq_record.id,
                             "size": len(seq_record.seq),
                             "fasta": os.path.basename(fasta_path)})
    return ret_list


main()
exit(0)
