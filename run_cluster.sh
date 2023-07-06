#!/bin/bash
num_jobs=250

snakemake --cluster-config cluster.yaml --cluster-sync "sbatch --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus} --partition {cluster.partition} --output {cluster.output} --error {cluster.error} --job-name {cluster.name} --wait" -j $num_jobs --latency-wait 40 --max-jobs-per-second 2 --max-status-checks-per-second 5 -T 3

