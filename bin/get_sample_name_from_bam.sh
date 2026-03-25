#!/bin/bash

BAM_FILE=$1

if [ -z "$BAM_FILE" ]; then
    echo "Usage: $0 <bam_file>"
    exit 1
fi

if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file '$BAM_FILE' does not exist."
    exit 1
fi

samtools view -H "$BAM_FILE" | grep "^@RG" | sed -n 's/.*SM:\([^[:space:]]\+\).*/\1/p' | head -n 1
