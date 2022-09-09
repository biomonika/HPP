version 1.0


workflow runHicHifiasm{
    input {
        String sampleID
        Array[File] hicForwardArray
        Array[File] hicReverseArray
        Array[File] bam_files
        Int threadCount
        Int memSize
        Int diskSizeGB=512
    }

    call formatToFastq {
        input:
            sampleID=sampleID,
            bam_files=bam_files,
            threadCount=threadCount,
            memSize=memSize,
            diskSizeGB=diskSizeGB  
    }

    call concatenateAndCompressHiFiReads {
        input:
            sampleID=sampleID,
            fastqArray=formatToFastq.fastqArray,
            threadCount=threadCount,
            memSize=memSize,
            diskSizeGB=diskSizeGB  
    }

    call concatenateHiCReads {
        input:
            sampleID=sampleID,
            hicForwardArray=hicForwardArray,
            hicReverseArray=hicReverseArray,
            threadCount=threadCount,
            memSize=memSize,
            diskSizeGB=diskSizeGB  
    }

    call trioHifiasm {
        input:
            hicForward=concatenateHiCReads.hicForward,
            hicReverse=concatenateHiCReads.hicReverse,
            HiFiReads=concatenateAndCompressHiFiReads.HiFiReads,
            threadCount=threadCount,
            memSize=memSize,
            diskSizeGB=diskSizeGB
    }

    call merylDiff {
        input:
            sampleID=sampleID,
            hap1=trioHifiasm.hap1,
            hap2=trioHifiasm.hap2,
            threadCount=threadCount,
            memSize=memSize,
            diskSizeGB=diskSizeGB
    }

    output {
        File only_hap1 = merylDiff.only_hap1
        File only_hap2 = merylDiff.only_hap2
    }

    parameter_meta {
        sampleID: " NIST ID (or Coriell ID) of the sample whose reads are going to be assembled"
        hicForward: "The Hi-C reads used for the haplotype phasing, paired forward reads in .gz format"
        hicReverse: "The Hi-C reads used for the haplotype phasing, paired reverse reads in .gz format"
        bam_files: "(Array of) HiFi reads in bam format used for the assembly with hifiasm"
        threadCount: "(default=48) The number of cores for running hifiasm"
        memSize: "(default=256) The memory size (GB) for running hifiasm"
        diskSizeGB: "(default=256) The local storage (GB) for running hifiasm"

    }
    meta {
        author: "Monika Cechova"
        email: "mcechova@ucsc.edu"
    }
}

task formatToFastq {
    input{
        String sampleID
        Array[File] bam_files
        Int threadCount
        Int memSize
        Int diskSizeGB
    }
    command <<<

        # exit when a command fails
        set -eux -o pipefail

        for x in ~{sep=' ' bam_files}
        do
            samtools bam2fq -@ ~{threadCount} ${x} > $(basename "${x}" .bam).fastq
        done;
   
    >>>

    output {
         Array[File] fastqArray = glob("*.fastq")
    }

    runtime {
        memory: memSize + "GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/samtools:1.6--hcd7b337_9"
        preemptible: 1
    }

}

task concatenateAndCompressHiFiReads {
    input{
        String sampleID
        Array[File] fastqArray
        Int threadCount
        Int memSize
        Int diskSizeGB
    }
    command <<<

        # exit when a command fails
        set -eux -o pipefail

        for x in ~{sep=' ' fastqArray}
        do
            cat ${x} >> ~{sampleID}.concHiFiReads.fq
        done;

        gzip ~{sampleID}.concHiFiReads.fq
   
    >>>

    output {
         File HiFiReads = "~{sampleID}.concHiFiReads.fq.gz"
    }

    runtime {
        memory: 16 + "GB"
        cpu: 1
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "docker.io/library/ubuntu:latest"
        preemptible: 1
    }

}

task concatenateHiCReads {
    input{
        String sampleID
        Array[File] hicForwardArray
        Array[File] hicReverseArray
        Int threadCount
        Int memSize
        Int diskSizeGB
    }
    command <<<

        # exit when a command fails
        set -eux -o pipefail

        #concatenate forward HiC reads
        for x in ~{sep=' ' hicForwardArray}
        do
            cat ${x} >> ~{sampleID}.concHicArrayForward.fq.gz
        done;

        #concatenate reverse HiC reads
        for x in ~{sep=' ' hicReverseArray}
        do
            cat ${x} >> ~{sampleID}.concHicArrayReverse.fq.gz
        done;
   
    >>>

    output {
         File hicForward = "~{sampleID}.concHicArrayForward.fq.gz"
         File hicReverse = "~{sampleID}.concHicArrayReverse.fq.gz"
    }

    runtime {
        memory: memSize + "GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "docker.io/library/ubuntu:latest"
        preemptible: 1
    }

}


task trioHifiasm {
    input{
        File hicForward
        File hicReverse
        File HiFiReads
        Int threadCount
        Int memSize
        Int diskSizeGB
    }
    command <<<

        # exit when a command fails
        set -eux -o pipefail

        hifiasm -o "hifiasm_output" -t ~{threadCount} --h1 ~{hicForward} --h2 ~{hicReverse} ~{HiFiReads}
        awk '/^S/{print ">"$2;print $3}' *.hap1.p_ctg.gfa > hap1.fa
        awk '/^S/{print ">"$2;print $3}' *.hap2.p_ctg.gfa > hap2.fa
       
    >>>

    output {
        File hap1 = "hap1.fa"
        File hap2 = "hap2.fa"
    }

    runtime {
        memory: memSize + "GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/hifiasm:0.16.1--h5b5514e_1"
        preemptible: 1
    }
}

task merylDiff {
    input{
        String sampleID
        File hap1
        File hap2
        Int threadCount
        Int memSize
        Int diskSizeGB
    }
    command <<<

        # exit when a command fails
        set -eux -o pipefail

        meryl threads=~{threadCount} count compress k=21 ~{hap1} output hap1.meryl
        meryl threads=~{threadCount} count compress k=21 ~{hap2} output hap2.meryl

        meryl difference hap1.meryl hap2.meryl output only_hap1.meryl
        meryl difference hap2.meryl hap1.meryl output only_hap2.meryl
       
        tar cvf ~{sampleID}.only_hap1.meryl.tar only_hap1.meryl
        tar cvf ~{sampleID}.only_hap2.meryl.tar only_hap2.meryl
       
    >>>

    output {
        File only_hap1 = "~{sampleID}.only_hap1.meryl.tar"
        File only_hap2 = "~{sampleID}.only_hap2.meryl.tar"
    }

    runtime {
        memory: memSize + "GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/meryl:1.3--h1b792b2_0"
        preemptible: 1
    }

}