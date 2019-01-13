#! /usr/bin/env python

"""
:Abstract: Format gene counts table from STAR for compatibility with voom-limma DGE
           workflow.
:Author: Akshay Paropkari
"""

import sys
import argparse
from os import walk
from os.path import join, realpath


def prog_options():
    parser = argparse.ArgumentParser(description="Format gene counts table from STAR for "
                                     "compatibility with voom-limma DGE workflow. This "
                                     "script will aggregate all counts files and produce "
                                     "one output file.")
    parser.add_argument("input_dir", nargs="*", help="Path to directory containing gene"
                        " counts tables from STAR run. The filenames must end with "
                        "ReadsPerGene.out.tab, which is default STAR behavior.")
    parser.add_argument("-o", "--output_prefix", default="./",
                        help="Path to directory where the output count table will be "
                        "saved. By default, the output count table will be saved in the "
                        "current working directory.")
    return parser.parse_args()


def main():
    args = prog_options()

    # Initialize variables to store gene counts as {gene_id: {sample_id: count}}
    counts = dict()
    samples = set()

    # Read in all ReadsPerGene.out.tab files
    for root, directory, files in walk(args.input_dir[0]):
        for f in files:
            try:
                assert f.endswith("ReadsPerGene.out.tab")
            except AssertionError:
                # not gene count files, skip
                continue
            else:
                # valid file encountered, process
                star_file = join(root, f)
                sample_id = f.replace("ReadsPerGene.out.tab", "")
                samples.add(sample_id)
                with open(star_file) as inf:
                    for line in inf:
                        try:
                            assert line.startswith("N_")
                        except AssertionError:
                            # gene counts line
                            line = line.strip().split("\t")
                            try:
                                assert line[0] in counts.keys()
                            except AssertionError:
                                # gene_id not initialized, create entry
                                counts[line[0]] = {sample_id: line[2]}
                            else:
                                # gene_id entry exists, update with sample_id and count
                                counts[line[0]][sample_id] = line[2]
                        else:
                            # stats line encountered, skip
                            continue

    # Write out counts to a file
    samples = sorted(samples)
    out_fnh = join(realpath(args.output_prefix), "gene_raw_counts.txt")
    with open(out_fnh, "w") as outf:
        outf.write("gene_id\t{}\n".format("\t".join(samples)))
        for gene_id in counts.keys():
            line_content = [gene_id]
            for sid in samples:
                try:
                    line_content.append(counts[gene_id][sid])
                except KeyError:
                    # sample doesn't have counts entry, initialize to 0
                    line_content.append("0")
                else:
                    continue
            outf.write("{}\n".format("\t".join(line_content)))

    print("Gene counts collected in {}".format(realpath(out_fnh)))


if __name__ == "__main__":
    sys.exit(main())
