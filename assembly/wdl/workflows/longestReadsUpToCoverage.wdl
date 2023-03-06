version 1.0


workflow longestReadsUpToCoverage{
    input {
        File fastq
        String fastq_name=basename(sub(sub(sub(fastq, "\\.gz$", ""), "\\.fasta$", ""), "\\.fa$", "")) #remove the file extension
        String genome_size
        String desired_coverage
        Int preemptible = 1
    }

    call runBioAwk {
        input:
            fastq=fastq,
            fastq_name=fastq_name,
            preemptible=preemptible
    }

    call findSubsetOfReads {
        input:
            read_lengths=runBioAwk.read_lengths,
            genome_size=genome_size,
            desired_coverage=desired_coverage,
            preemptible=preemptible
    }

    call subsampleFastq {
        input:
            fastq=fastq,
            fastq_name=fastq_name,
            read_names=findSubsetOfReads.read_names,
            preemptible=preemptible
    }

    output {
        File subsampledLongestReads = subsampleFastq.subsampledLongestReads
        File summary = findSubsetOfReads.summary
    }

    parameter_meta {
        fastq: "The fastq reads to be subsampled"
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
        String fastq_name
        Int memSizeGB = 32
        Int preemptible
    }

    Int file_size = ceil(size(fastq, "GB"))
    Int diskSizeGB = 3 * file_size

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
        Int memSizeGB=8
        Int diskSizeGB=8
        Int preemptible
    }

    command <<<
        R --no-save --args ~{read_lengths} ~{genome_size} ~{desired_coverage} <<Rscript
        args <- commandArgs(trailingOnly = TRUE)
        print(args)
        filename<-as.character(args[1])
        genome_size<-as.numeric(as.character(args[2]))
        desired_coverage<-as.numeric(as.character(args[3]))
        read_lengths <-as.data.frame(read.table(filename))
        colnames(read_lengths)<-c("name","length")

        lengths<-read_lengths["length"]
        names<-read_lengths["name"]

        sufficient_reads<-(cumsum(lengths)/genome_size)>desired_coverage #cummulative sum of reads sorted from the longest to the shortest
        smallest_number_of_reads_needed<-min(which(sufficient_reads==TRUE))

        if ((genome_size*desired_coverage) > sum(lengths)) {

            summary<-"The desired coverage for this genome size exceeds the number of available basepairs. The subsampling is not possible."

            sufficient_subset<-names

        } else {

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
        String fastq_name
        File read_names
        Int memSizeGB = 32
        Int preemptible
    }

    Int file_size = ceil(size(fastq, "GB"))
    Int diskSizeGB = 3 * file_size

    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        seqtk subseq ~{fastq} ~{read_names} | gzip >~{fastq_name}.subsampledLongestReads.fastq.gz
        

    >>>

    output {
        File subsampledLongestReads = "${fastq_name}.subsampledLongestReads.fastq.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible : preemptible
        docker: "quay.io/biocontainers/seqtk:1.3--hed695b0_2"
    }
}



