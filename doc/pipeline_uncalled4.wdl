version 1.0

workflow ont_uncalled4 {
    input {
        File bam_file,
        File pod5_dir,
        String sample_id,
        File ref_genome,
    }

    call Uncalled4Align {
        input:
            bam_file = bam_file,
            pod5_dir = pod5_dir,
            sample_id = sample,
            ref_genome = ref_genome
    }

    output {
        File uncalled4_bam = Uncalled4Align.bam
        File uncalled4_bai = Uncalled4Align.bai
        File uncalled4_log = Uncalled4Align.log
    }

task Uncalled4Align {
    input {
        File bam_file,
        File pod5_dir,
        String sample_id,
        File ref_genome,
        Int cpus
    }

    command <<<
    set -euo pipefail
    filename=$(basename "~{pod5_file}" .pod5)
    
    mkdir -p uncalled4_bam

    echo $(date): running uncalled4
    /usr/bin/time --verbose \
    uncalled4 align --rna \
        --ref ${ref_genome} \
        --reads ${pod5_dir} \
        --recursive \
        --bam-in ${input_bam} \
        --bam-out ${unsorted_bam} \
        -p ${cpus}

    echo $(date): sorting BAM
    module load samtools
    samtools sort --output-fmt BAM --write-index \
    -o ${output_bam} \
    ${unsorted_bam}

    echo $(date): finished running uncalled4
    >>>

    output {
        File uncalled4_bam = Uncalled4Align.bam
        File uncalled4_bai = Uncalled4Align.bai
        File uncalled4_log = Uncalled4Align.log
    }

    runtime {
        cpu: cpus
        gpu: false
        memory: "16GB"
        maxRunTime: 172800 #48 hours (48 * 3600 seconds)
        runtime_minutes: 1 #47 hours (47 * 60 minutes)
        docker: "FIXME"
    }
}
