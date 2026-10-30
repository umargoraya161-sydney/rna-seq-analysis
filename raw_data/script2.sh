

#!/bin/bash

INPUT_DIR=$(pwd)
OUTPUT_DIR="${INPUT_DIR}"

for FASTQ_FILE in "${INPUT_DIR}"/*_trimmed.fastq
do
    FILENAME=$(basename "${FASTQ_FILE}" .fastq)

    hisat2 -q \
        --rna-strandness R \
        -x "${INPUT_DIR}/grch38/genome" \
        -U "${FASTQ_FILE}" \
        -p 4 \
        | samtools sort -o "${OUTPUT_DIR}/${FILENAME}.bam"

done