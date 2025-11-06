## Bioinformatics Toolkits

### Question: When should you use -a, -A, -g, or -G with Cutadapt?

In Cutadapt, adapter or primer sequences can be trimmed from sequencing reads. The choice of parameter depends on the read type and where the adapter appears:


| Parameter | Use case                                                                         |
| --------- | -------------------------------------------------------------------------------- |
| `-a`      | Trim a **3' adapter** from **read 1** (forward read)                             |
| `-A`      | Trim a **3' adapter** from **read 2** (reverse read) – only for paired-end reads |
| `-g`      | Trim a **5' adapter** from **read 1** (forward read)                             |
| `-G`      | Trim a **5' adapter** from **read 2** (reverse read) – only for paired-end reads |


Tip: Use -a/-A for standard 3' adapter trimming, and -g/-G if you know your primers/adapters are at the 5' end.

You can also use cutadapt_summary.py to detect adapters automatically:

```shell
python scripts/cutadapt_summary.py --help
```

Usage:


```shell
usage: cutadapt_summary.py [-h] -i INPUT_DIR -o OUTPUT [--threads THREADS]
                           [--adapter_fwd ADAPTER_FWD]
                           [--adapter_rev ADAPTER_REV]

Detect adapters/primers via cutadapt --report=minimal (datavzrd-ready)

options:
  -h, --help            Show this help message and exit
  -i, --input_dir INPUT_DIR
                        Input directory containing *_1.fastq.gz and *_2.fastq.gz files
  -o, --output OUTPUT   Output TSV file for visualization
  --threads THREADS     Number of threads for cutadapt
  --adapter_fwd ADAPTER_FWD
                        Forward adapter/primer sequence
  --adapter_rev ADAPTER_REV
                        Reverse adapter/primer sequence
```

Recommendation:

If you are unsure which adapters are present, first run cutadapt_summary.py without specifying adapters. It will suggest the most likely adapters, which you can then trim using Cutadapt with the appropriate -a/-A/-g/-G options.
