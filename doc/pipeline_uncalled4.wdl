version 1.2

workflow ont_uncalled4 {
    input {
        File bam_file
        Directory pod5_dir
        String sample_id
        File ref_genome
        Int cpus
    }

    call Uncalled4Align {
        input:
            bam_file = bam_file,
            pod5_dir = pod5_dir,
            sample_id = sample_id,
            ref_genome = ref_genome,
            cpus = cpus
    }

    output {
        File uncalled4_bam = Uncalled4Align.uncalled4_bam
        File uncalled4_bai = Uncalled4Align.uncalled4_bai
        File uncalled4_log = Uncalled4Align.uncalled4_log
    }
}

task Uncalled4Align {
    input {
        File bam_file
        Directory pod5_dir
        String sample_id
        File ref_genome
        Int cpus
    }

    command <<<
    set -euo pipefail
    mkdir -p uncalled4_bam
    output_base="~{sample_id}"
    unsorted_bam="${output_base}.unsorted.bam"
    output_bam="${output_base}.bam"
    log_file="${output_base}.log"

    {
        echo $(date): running uncalled4

        echo pod5 dir contains:
        ls "~{pod5_dir}"
        echo

        time uncalled4 align --rna \
            --ref "~{ref_genome}" \
            --reads "~{pod5_dir}" \
            --recursive \
            --bam-in "~{bam_file}" \
            --bam-out "${unsorted_bam}" \
            -p ~{cpus}

        echo $(date): sorting BAM
        samtools sort --output-fmt BAM --write-index -o "${output_bam}" "${unsorted_bam}"
        rm ${unsorted_bam}

        echo $(date): finished running uncalled4
    } 2>&1 | tee "${log_file}"
    >>>

    output {
        File uncalled4_bam = sample_id + ".bam"
        File uncalled4_bai = sample_id + ".bam.csi"
        File uncalled4_log = sample_id + ".log"
    }

    runtime {
        cpu: cpus
        gpu: false
        memory: "16GB"
        maxRunTime: 172800 #48 hours (48 * 3600 seconds)
        runtime_minutes: 1 #47 hours (47 * 60 minutes)
        docker: "joshtburdick/uncalled4@sha256:c130f85cb2e7151dddb849d882ef3ac16b8f8283573bf53f3fdec051de2bcb01"
    }
}
