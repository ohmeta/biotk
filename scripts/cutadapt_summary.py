#!/usr/bin/env python3
import subprocess
import glob
import os
import argparse
import pandas as pd

# ========== 函数定义 ==========

def run_cutadapt(r1, r2, method, seq1, seq2, threads=4):
    """
    运行 cutadapt 检测 adapter/primer，返回解析后的 DataFrame 格式的结果
    method: 'aA' or 'gG'
    """
    params = ["cutadapt", "-j", str(threads),
              "--report=minimal", "-o", "/dev/null", "-p", "/dev/null"]
    if method == "aA":
        params += ["-a", seq1, "-A", seq2]
    elif method == "gG":
        params += ["-g", seq1, "-G", seq2]
    else:
        raise ValueError("method must be 'aA' or 'gG'")
    params += [r1, r2]

    try:
        result = subprocess.run(params, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ cutadapt failed: {r1}")
        return None

    lines = [line.strip() for line in result.stdout.splitlines()]
    if len(lines) < 2:
        print(f"⚠️ cutadapt report missing for {r1}")
        return None

    header = lines[0].split()
    values = lines[1].split()
    df = pd.DataFrame([values], columns=header)
    return df


def extract_summary(df, prefix):
    """
    从 cutadapt TSV 结果中提取关键信息，并加上前缀
    """
    status = str(df["status"].iloc[0])
    in_reads = int(df["in_reads"].iloc[0])
    in_bp = int(df["in_bp"].iloc[0])
    too_long = int(df["too_long"].iloc[0])
    too_short = int(df["too_short"].iloc[0])
    too_many_n = int(df["too_many_n"].iloc[0])

    out_reads = int(df["out_reads"].iloc[0])
    out_bp = int(df["out_bp"].iloc[0])
    out2_bp = int(df["out2_bp"].iloc[0])
    qualtrim_bp = int(df["qualtrim_bp"].iloc[0])
    qualtrim2_bp = int(df["qualtrim2_bp"].iloc[0])
    w_adapters = int(df["w/adapters"].iloc[0])
    w_adapters2 = int(df["w/adapters2"].iloc[0])

    percent_r1 = (w_adapters / in_reads * 100) if in_reads > 0 else 0.0
    percent_r2 = (w_adapters2 / in_reads * 100) if in_reads > 0 else 0.0

    return {
        f"{prefix}_status": status,
        f"{prefix}_in_reads": in_reads,
        f"{prefix}_out_reads": out_reads,

        f"{prefix}_w_adapters": w_adapters,
        f"{prefix}_w_adapters2": w_adapters2,

        f"{prefix}_too_long": too_long,
        f"{prefix}_too_short": too_short,
        f"{prefix}_too_many_n": too_many_n,

        f"{prefix}_in_bp": in_bp,
        f"{prefix}_out_bp": out_bp,
        f"{prefix}_out2_bp": out2_bp,
        f"{prefix}_qualtrim_bp": qualtrim_bp,
        f"{prefix}_qualtrim2_bp": qualtrim2_bp,

        f"{prefix}_adapter_percent_r1": percent_r1,
        f"{prefix}_adapter_percent_r2": percent_r2,
    }


def main():
    parser = argparse.ArgumentParser(description="Detect adapters/primers via cutadapt --report=tsv (datavzrd-ready)")
    parser.add_argument("-i", "--input_dir", required=True, help="Input directory containing *_1.fastq.gz and *_2.fastq.gz files")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file for visualization")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for cutadapt")
    parser.add_argument("--adapter_fwd", default="AGRGTTYGATYMTGGCTCAG", help="Forward adapter/primer sequence")
    parser.add_argument("--adapter_rev", default="CTGCWGCCHCCCGTAGG", help="Reverse adapter/primer sequence")
    args = parser.parse_args()

    records = []
    fastq_files = sorted(glob.glob(os.path.join(args.input_dir, "*_1.fastq.gz")))
    if not fastq_files:
        print("⚠️ No FASTQ files found.")
        return

    for r1_path in fastq_files:
        r2_path = r1_path.replace("_1.fastq.gz", "_2.fastq.gz")
        if not os.path.exists(r2_path):
            print(f"⚠️ Missing pair for {r1_path}")
            continue

        sample = os.path.basename(r1_path).split("_1.fastq.gz")[0]
        print(f"🔍 Analyzing sample: {sample}")

        df_aA = run_cutadapt(r1_path, r2_path, "aA", args.adapter_fwd, args.adapter_rev, args.threads)
        df_gG = run_cutadapt(r1_path, r2_path, "gG", args.adapter_fwd, args.adapter_rev, args.threads)

        if df_aA is None or df_gG is None:
            continue

        rec = {"sample": sample}
        rec.update(extract_summary(df_aA, "three_prime"))
        rec.update(extract_summary(df_gG, "five_prime"))
        records.append(rec)

    df_out = pd.DataFrame(records)
    df_out.to_csv(args.output, sep="\t", index=False)
    print(f"✅ Report written to: {args.output}")
    print("💡 Visualize with datavzrd:")
    print(f"   datavzrd --input {args.output} --output adapter_report.html")


if __name__ == "__main__":
    main()
