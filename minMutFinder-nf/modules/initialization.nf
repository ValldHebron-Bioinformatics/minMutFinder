#!/usr/bin/env nextflow

// This file is part of minMutFinder.
//
// minMutFinder is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// minMutFinder is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with minMutFinder. If not, see <https://www.gnu.org/licenses/>.
//
// Copyright (C) 2024 Ignasi Prats Méndez

process dirCreator {
    errorStrategy 'terminate'
    input:
    path out_path

    output:
    stdout

    script:
    """
    #!/bin/bash
    echo "Creating directories..."
    QC_DIR="${out_path}/qc"
    # FASTQ="${out_path}/fastq"
    ASSEMBLY="${out_path}/assembly"
    VARIANT_CALLING="${out_path}/variant_calling"
    MUTATIONS="${out_path}/mutations"
    PLOTS="${out_path}/plots"

    # Create directories if they do not exist
    for DIR in "\${QC_DIR}" "\${ASSEMBLY}" "\${VARIANT_CALLING}" "\${MUTATIONS}" "\${PLOTS}"; do
        if [[ -d "\${DIR}" ]]; then
            echo -e "\${DIR} output folder already exists"
        else
            mkdir -p "\${DIR}"
        fi
    done

    # Create stopper file in mutations directory
    touch stopper
    """
}

process refCheck {
    errorStrategy 'terminate'

    input:
    path out_path
    file ref_seq

    output:
    tuple file("ref_seq.fasta"), file("ref_seq.fasta.fai")

    script:
    """
    ASSEMBLY="${out_path}/assembly"

    # Wait until the assembly directory is created
    while [[ ! -d "\${ASSEMBLY}" ]]; do
        echo -e "\${ASSEMBLY} output folder not created"
    done

    # Run reference sequence check and prepare reference sequence
    python3 ${params.project_data}/bin/ref_check.py --ref-seq "${ref_seq}"
    cp ${ref_seq} ref_seq.fasta
    samtools faidx ref_seq.fasta
    """
}
