version 1.0


workflow longestReadsUpToCoverage{
    input {
        Array[File] fastq_files
        String genome_size
        String desired_coverage
        Int preemptible = 1
    }

    scatter (fastq in fastq_files) {
        call runBioAwk {
            input:
                fastq=fastq,
                memSizeGB=32,
                preemptible=preemptible
        }

        call findSubsetOfReads {
            input:
                read_lengths=runBioAwk.read_lengths,
                genome_size=genome_size,
                desired_coverage=desired_coverage,
                memSizeGB=32,
                preemptible=preemptible
        }

        call subsampleFastq {
            input:
                fastq=fastq,
                read_names=findSubsetOfReads.read_names,
                memSizeGB=32,
                preemptible=preemptible
        }

        call compressFastq {
            input:
                fastqSumbsampled=subsampleFastq.subsampledLongestReads,
                memSizeGB=32,
                threadCount=32,
                preemptible=preemptible
        }
    }

    call concatenate {
        input:
            subsampled_files=compressFastq.compressedSubsampledLongestReads,
            memSizeGB=32,
            preemptible=1
    }

    output {
        File compressFastq = concatenate.subsampled_final
        Array[File] summary = findSubsetOfReads.summary
    }

    parameter_meta {
        fastq_files: "The fastq reads to be subsampled"
        genome_size: "The genome size of the reference is basepairs [bps]"
        desired_coverage: "The final desired coverage for the subset of the longest reads. Please use an integer value."
    }
    meta {
        author: "Monika Cechova"
        email: "mcechova@ucsc.edu"
    }
}

task runBioAwk {
    input{
        File fastq
        Int memSizeGB
        Int diskSizeGB = ceil(size(fastq, "GB")) + 64
        Int preemptible
    }

    String fastq_name=basename(sub(sub(sub(fastq, "\\.gz$", ""), "\\.fasta$", ""), "\\.fa$", "")) #remove the file extension

    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        bioawk -c fastx '{ print $name, length($seq)}' ~{fastq} | sort -rgk2 >~{fastq_name}.read_lengths.txt

    >>>

    output {
        File read_lengths="${fastq_name}.read_lengths.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible : preemptible
        docker: "quay.io/biocontainers/bioawk:1.0--hed695b0_5"
    }
}

task findSubsetOfReads {
    input{
        File read_lengths
        String genome_size
        String desired_coverage
        Int memSizeGB
        Int diskSizeGB = 64
        Int preemptible
    }

    command <<<
        R --no-save --args ~{read_lengths} ~{genome_size} ~{desired_coverage} <<Rscript
        args <- commandArgs(trailingOnly = TRUE)
        print(args)
        filename<-as.character(args[1])
        genome_size<-as.numeric(as.character(args[2]))
        desired_coverage<-as.numeric(as.character(args[3]))
        read_lengths <-as.data.frame(read.table(filename, header=FALSE, sep="\t"))
        colnames(read_lengths)<-c("name","length")

        lengths<-read_lengths["length"]
        names<-read_lengths["name"]

        if ((genome_size*desired_coverage) > sum(lengths/1)) {

            summary<-"The desired coverage for this genome size exceeds the number of available basepairs. The subsampling is not possible."
            sufficient_subset<-names

        } else {
            sufficient_reads<-(cumsum(lengths/1)/genome_size)>desired_coverage #cummulative sum of reads sorted from the longest to the shortest
            smallest_number_of_reads_needed<-min(which(sufficient_reads==TRUE))
            summary<-paste("With genome size",genome_size,"the total number of the longest reads needed is", smallest_number_of_reads_needed, "with the total cummulative length of",sum(head(lengths,smallest_number_of_reads_needed)),"bps to achieve the minimum coverage of",desired_coverage,"in the subsampled dataset.")

            sufficient_subset<-head(names,smallest_number_of_reads_needed)
            
        }

        write.table(summary,file="summary.txt",quote=FALSE,col.names=FALSE, row.names=FALSE)
        write.table(sufficient_subset,file="longest_reads.ids.txt",quote=FALSE,col.names=FALSE, row.names=FALSE) 

        Rscript

    >>>

    output {
        File read_names="longest_reads.ids.txt"
        File summary="summary.txt"
    }

    runtime {
        docker: "rocker/tidyverse"
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
    }
}


task subsampleFastq {
    input{
        File fastq
        File read_names
        Int memSizeGB
        Int diskSizeGB = ceil(size(fastq, "GB")) * 2 + 64
        Int preemptible
    }

    String fastq_name=basename(sub(sub(sub(fastq, "\\.gz$", ""), "\\.fasta$", ""), "\\.fa$", "")) #remove the file extension

    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        seqtk subseq ~{fastq} ~{read_names} >~{fastq_name}.subsampledLongestReads.fastq
        

    >>>

    output {
        File subsampledLongestReads = "${fastq_name}.subsampledLongestReads.fastq"
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible : preemptible
        docker: "quay.io/biocontainers/seqtk:1.3--hed695b0_2"
    }
}

task compressFastq {
    input{
        File fastqSumbsampled
        Int memSizeGB
        Int threadCount
        Int diskSizeGB = ceil(size(fastqSumbsampled, "GB")) * 2 + 64
        Int preemptible
    }

    String fastq_name=basename(fastqSumbsampled)

    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        cat ~{fastqSumbsampled} | pigz -p ~{threadCount} >~{fastq_name}.gz        

    >>>

    output {
        File compressedSubsampledLongestReads = "${fastq_name}.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible : preemptible
        docker: "quay.io/biocontainers/pigz:2.3.4"
    }
}

task concatenate {
    input{
        Array[File] subsampled_files
        Int memSizeGB
        Int preemptible
    }

    Int file_size = ceil(size(subsampled_files, "GB"))
    Int diskSizeGB = 2 * file_size + 64

    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        cat ~{sep=' ' subsampled_files} >concatenated.fastq.gz

    >>>

    output {
        File subsampled_final="concatenated.fastq.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible : preemptible
        docker: "ubuntu:18.04"
    }
}



