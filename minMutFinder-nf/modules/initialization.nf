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

    input:
    path out_path

    output:
    path "${out_path}"
    path "${out_path}/fastq"
    path "${out_path}/assembly"
    path "${out_path}/assembly/references"
    path "${out_path}/variant_calling"
    path "${out_path}/mutations"
    path "${out_path}/plots"
    path "${out_path}/qc"

    script:
    """
    #!/bin/bash
    echo Creating directories..
    DIR_SAMPLE="${out_path}"
    if [[ -d "\${DIR_SAMPLE}" ]]; then echo -e ""\${DIR_SAMPLE}" output folder already exists"; else mkdir "\${DIR_SAMPLE}"; fi
    QC_DIR="${out_path}/qc"
    if [[ -d "\${QC_DIR}" ]]; then echo -e ""\${QC_DIR}" output folder already exists"; else mkdir "\${QC_DIR}"; fi
    FASTQ="${out_path}/fastq"
    if [[ -d "\${FASTQ}" ]]; then echo -e ""\${FASTQ}" output folder already exists"; else mkdir "\${FASTQ}"; fi
    ASSEMBLY="${out_path}/assembly"
    REF_PATH="${out_path}/assembly/references"
    if [[ -d "\${ASSEMBLY}" ]]; then 
        echo -e ""\${ASSEMBLY}" output folder already exists"
    else 
        mkdir "\${ASSEMBLY}";
        mkdir "\${REF_PATH}";
    fi
    VARIANT_CALLING="${out_path}/variant_calling"
    if [[ -d "\$VARIANT_CALLING" ]]; then echo -e ""\${VARIANT_CALLING}" output folder already exists"; else mkdir "\${VARIANT_CALLING}"; fi
    MUTATIONS="${out_path}/mutations"
    if [[ -d "\${MUTATIONS}" ]]; then echo -e ""\${MUTATIONS}" output folder already exists"; else mkdir \${MUTATIONS}; fi
    touch ${out_path}/mutations/stopper
    PLOTS="${out_path}/plots"
    if [[ -d "\${PLOTS}" ]]; then echo -e ""\${PLOTS}" output folder already exists"; else mkdir "\${PLOTS}"; fi
    """
}


process refCheck {
    input:
    file ref_seq

    output:
    stdout

    script:
    """
    python3 $projectDir/bin/ref_check.py --ref-seq "${params.ref_seq}"
    """
}
