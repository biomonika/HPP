#!/usr/bin/env python3
import argparse, gzip, sys
from collections import namedtuple
import pysam

FQEntry = namedtuple("FQEntry", ["seq", "qual"])

def open_text_maybe_gz(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")

def norm_name(s: str) -> str:
    """
    Normalize read names so FASTQ <-> BAM match:
    - drop leading '@' if present
    - trim
    - keep only first whitespace-delimited token
    """
    n = s.strip()
    if n.startswith("@"):
        n = n[1:]
    n = n.split()[0]
    return n

def load_fastq_to_dict(fq_path):
    fq = {}
    with open_text_maybe_gz(fq_path) as fh:
        while True:
            header = fh.readline()
            if not header: break

            seq  = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not (seq and plus and qual):
                raise ValueError(f"Truncated FASTQ at {fq_path}")
            if not header.startswith("@"):
                raise ValueError(f"Bad FASTQ header: {header.strip()}")

            # normalize read name
            name = norm_name(header)
            seq, qual = seq.strip(), qual.strip()
            if len(seq) != len(qual):
                raise ValueError(f"FASTQ {name} has seq/qual length mismatch: {len(seq)} vs {len(qual)}")
            fq[name] = FQEntry(seq, qual)
            #print(f"Name in fastq: {name} length: {len(seq)}")
    
    print(f"Number of entries in fastq: {len(fq)}")
    return fq

def replace_seq_qual(bam_in, bam_out, fq_dict, max_rel_diff=0.10):
    replaced = skipped_no_fastq = wrong_relative_sequence_length = total = 0

    no_fastq_log = open("skipped_no_fastq.txt", "w")
    len_gate_log = open("wrong_relative_sequence_length.txt", "w")

    for read in bam_in.fetch(until_eof=True):
        print("\n")
        total += 1
        qname = norm_name(read.query_name)
        fq_entry = fq_dict.get(qname)

        print(f"Name in bam file: {qname}")

        if fq_entry is None:
            print("skipped_no_fastq")
            skipped_no_fastq += 1
            no_fastq_log.write(qname + "\n")
            continue

        len_bam = len(read.query_sequence or "")
        len_fq = len(fq_entry.seq)

        print(f"length in bam: {len_bam}; length in fastq: {len_fq}")

        rel_diff = abs(len_bam - len_fq) / max(len_bam, len_fq, 1)
        print(f"rel_diff: {rel_diff}")

        if rel_diff <= max_rel_diff:
            #the lengths are not too far off, so the replacement can be considered
            replaced += 1 #we are decided we will do the replacement, either with padding or "as is""

            if len_fq < len_bam:
                #if the sequence from fastq is shorter, it needs to be padded
                pad_len = len_bam - len_fq
                fq_entry = FQEntry(
                    fq_entry.seq + ("N" * pad_len),
                    fq_entry.qual + ("I" * pad_len)  # "I" as a quality score default
                )
                #replace the sequence and the quality from bam with a padded sequence
                read.query_sequence = fq_entry.seq
                read.query_qualities = pysam.qualitystring_to_array(fq_entry.qual)
             
                bam_out.write(read) #ONLY WRITE INTO THE BAM FILE NOW
                print(f"Wrote padded entry.")
                print(f"length of the written query: {len(read.query_sequence)}")
            else:
                #fastq sequence is longer, so it can be just swapped "as is"
                #we know the total lengts are not too far off so this is okay

                read.query_sequence = fq_entry.seq
                read.query_qualities = pysam.qualitystring_to_array(fq_entry.qual)

                bam_out.write(read) #ONLY WRITE INTO THE BAM FILE NOW
                print(f"Wrote entry \"as is\".")
        else:
            print("wrong_relative_sequence_length")
            wrong_relative_sequence_length += 1
            len_gate_log.write(qname + "\n")

    no_fastq_log.close()
    len_gate_log.close()

    print(
        f"[summary] total={total} | replaced={replaced} | "
        f"skipped_no_fastq={skipped_no_fastq} | "
        f"wrong_relative_sequence_length={wrong_relative_sequence_length}",
        file=sys.stderr,
    )

def main():
    p = argparse.ArgumentParser("Replace BAM seq/qual from FASTQ by read name (Nanopore, single-end).")
    p.add_argument("-b","--bam", required=True)
    p.add_argument("-f","--fastq", required=True)
    p.add_argument("-o","--out", required=True)
    p.add_argument("--max-rel-diff", type=float, default=0.10)
    args = p.parse_args()

    fq_dict = load_fastq_to_dict(args.fastq)

    # Allow BAMs without @SQ
    with pysam.AlignmentFile(args.bam, "rb", check_sq=False) as bam_in:
        with pysam.AlignmentFile(args.out, "wb", header=bam_in.header) as bam_out:
            replace_seq_qual(bam_in, bam_out, fq_dict, max_rel_diff=args.max_rel_diff)

    try:
        pysam.index(args.out)
    except Exception as e:
        print(f"[warn] Could not index output BAM: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
