version 1.0


workflow estimateTelomeres {
  input {
    Array[File] bam_files
    Int flank = 100000
    Int minReadLength = 15000
    Int preemptible = 1
  }

  scatter (bam_file in bam_files) {

      call ExtractFasta {
        input:
            bam_file=bam_file,
            flank=flank,
            memSizeGB=2,
            threadCount=1,
            preemptible=preemptible
      }

      call runSeqtk {
        input:
            fasta_file_p_forward=ExtractFasta.fasta_file_p_forward,
            fasta_file_p_reverse=ExtractFasta.fasta_file_p_reverse,
            fasta_file_q_forward=ExtractFasta.fasta_file_q_forward,
            fasta_file_q_reverse=ExtractFasta.fasta_file_q_reverse,
            preemptible=preemptible
      }

      call filterTelomeres {
        input:
            bed_file_p_forward=runSeqtk.bed_file_p_forward,
            bed_file_p_reverse=runSeqtk.bed_file_p_reverse,
            bed_file_q_forward=runSeqtk.bed_file_q_forward,
            bed_file_q_reverse=runSeqtk.bed_file_q_reverse,
            minReadLength=minReadLength,
            preemptible=preemptible
      }

      call calculateMedian {
        input:
            telomere_lengths_p_forward=filterTelomeres.telo_p_forward,
            telomere_lengths_p_reverse=filterTelomeres.telo_p_reverse,
            telomere_lengths_q_forward=filterTelomeres.telo_q_forward,
            telomere_lengths_q_reverse=filterTelomeres.telo_q_reverse,
            preemptible=preemptible
      }

  }

    output {
        Array[File] summary = calculateMedian.summary
    }

    parameter_meta {
        bam_file: " The input bam file."
        minReadLength: "The minimal length for a read to be considered."
        flank: "The maximal distance from either end of a chromosome, in which a read will be used if overlapping such coordinates."
    }
    meta {
        author: "Monika Cechova"
        email: "mcechova@ucsc.edu"
    }

}

task ExtractFasta {
  input {
    File bam_file
    Int flank
    Int memSizeGB
    Int threadCount
    Int diskSizeGB=ceil(size(bam_file, "GB")) * 4
    Int preemptible
    }

  command <<<

    #index the bam file so that random retrieval is possible
    samtools index ~{bam_file}

    #derive chromosome name from the first entry, this assumes bam files to be correctly separated by chromosome
    chromosome=`samtools view ~{bam_file} | head -n 1 | cut -f 3`
    chromosome_length=`samtools view -H ~{bam_file} | grep "^@SQ" | sed 's/SN://;s/LN://' | sed -e 's/\t/ /g' | grep "${chromosome} " | cut -d' ' -f3`

    echo "chromosome:" $chromosome
    echo "chromosome_length:" $chromosome_length

    start=1
    last_start=$((${chromosome_length} - ~{flank} + 1))

    #add headers first
    samtools view -H ~{bam_file} > p_forward.sam
    samtools view -H ~{bam_file} > p_reverse.sam
    samtools view -H ~{bam_file} > q_forward.sam
    samtools view -H ~{bam_file} > q_reverse.sam

    echo "Coordinates to use for subsetting:"
    echo "${chromosome}:$start-~{flank}"
    echo "${chromosome}:$last_start-"

    # Extract telomeric regions at the start and end of the chromosome

    #forward
    samtools view -F 16 ~{bam_file} "${chromosome}:$start-~{flank}" >> p_forward.sam
    #reverse
    samtools view -f 16 ~{bam_file} "${chromosome}:$start-~{flank}" >> p_reverse.sam

    #forward
    samtools view -F 16 ~{bam_file} "${chromosome}:$last_start-" >> q_forward.sam
    #reverse
    samtools view -f 16 ~{bam_file} "${chromosome}:$last_start-" >> q_reverse.sam


    sam_files=("p_forward.sam" "p_reverse.sam" "q_forward.sam" "q_reverse.sam")

    # Loop over each SAM file
    for sam_file in "${sam_files[@]}"; do
        base_name="${sam_file%.sam}"

        # Convert to BAM, filter for MAPQ >= 10, keep only primary alignments
        samtools view -q 10 -F 2308 -Sb "$sam_file" > "${base_name}.bam"
        samtools sort "${base_name}.bam" -o "${base_name}_sorted.bam"
        samtools index "${base_name}_sorted.bam"

        # Convert to FASTA format
        samtools fasta "${base_name}_sorted.bam" > "${base_name}.fasta"

        #remove unnecessary files
        rm ${sam_file} ${base_name}.bam ${base_name}_sorted.bam
    done


  >>>

  output {
    File fasta_file_p_forward = "p_forward.fasta"
    File fasta_file_p_reverse = "p_reverse.fasta"
    File fasta_file_q_forward = "q_forward.fasta"
    File fasta_file_q_reverse = "q_reverse.fasta"
  }

  runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible : preemptible
        docker: "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
  }
}

task runSeqtk {
    input{
        File fasta_file_p_forward
        File fasta_file_p_reverse
        File fasta_file_q_forward
        File fasta_file_q_reverse
        Int memSizeGB = 2
        Int preemptible
    }
    
    command <<<

      #search the fasta_file for the telomeric repeat
      seqtk telo -s 10 ~{fasta_file_p_forward} > reads.telo.p_forward.bed 2> reads.telo.p_forward.count
      seqtk telo -s 10 ~{fasta_file_p_reverse} > reads.telo.p_reverse.bed 2> reads.telo.p_reverse.count
      seqtk telo -s 10 ~{fasta_file_q_forward} > reads.telo.q_forward.bed 2> reads.telo.q_forward.count
      seqtk telo -s 10 ~{fasta_file_q_reverse} > reads.telo.q_reverse.bed 2> reads.telo.q_reverse.count
        
    >>>

    output {
        File bed_file_p_forward = "reads.telo.p_forward.bed"
        File bed_file_p_reverse = "reads.telo.p_reverse.bed"
        File bed_file_q_forward = "reads.telo.q_forward.bed"
        File bed_file_q_reverse = "reads.telo.q_reverse.bed"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "quay.io/biocontainers/seqtk:1.4--he4a0461_2"
    }

}

task filterTelomeres {
    input{
        File bed_file_p_forward
        File bed_file_p_reverse
        File bed_file_q_forward
        File bed_file_q_reverse
        Int minReadLength
        Int memSizeGB = 2
        Int preemptible
    }
    
    command <<<

    #filter for minimum read length; also, the telomeric repeat must extend to the very edge of the read
   
    #p-arms
    cat ~{bed_file_p_forward} | awk -v minReadLength="$minReadLength" '$4 >= minReadLength && ($2 == 0 || $3 == $4)' | awk '{print ($3-$2)}' > "bed_file_p_forward_filtered.lengths.txt"
    
    cat ~{bed_file_p_reverse} | awk -v minReadLength="$minReadLength" '$4 >= minReadLength && ($2 == 0 || $3 == $4)' | awk '{print ($3-$2)}' > "bed_file_p_reverse_filtered.lengths.txt"

    #q-arms
    cat ~{bed_file_q_forward} | awk -v minReadLength="$minReadLength" '$4 >= minReadLength && ($2 == 0 || $3 == $4)' | awk '{print ($3-$2)}' > "bed_file_q_forward_filtered.lengths.txt"
    
    cat ~{bed_file_q_reverse} | awk -v minReadLength="$minReadLength" '$4 >= minReadLength && ($2 == 0 || $3 == $4)' | awk '{print ($3-$2)}' > "bed_file_q_reverse_filtered.lengths.txt"
    
    >>>

    output {
        File telo_p_forward = "bed_file_p_forward_filtered.lengths.txt"
        File telo_p_reverse = "bed_file_p_reverse_filtered.lengths.txt"
        File telo_q_forward = "bed_file_q_forward_filtered.lengths.txt"
        File telo_q_reverse = "bed_file_q_reverse_filtered.lengths.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "quay.io/biocontainers/bioawk:1.0--hed695b0_5"
    }
}

task calculateMedian {
    input{
        File telomere_lengths_p_forward
        File telomere_lengths_p_reverse
        File telomere_lengths_q_forward
        File telomere_lengths_q_reverse
        Int preemptible
    }

    command <<<
        R --no-save --args ~{telomere_lengths_p_forward} ~{telomere_lengths_p_reverse} ~{telomere_lengths_q_forward} ~{telomere_lengths_q_reverse} <<Rscript
        args <- commandArgs(trailingOnly = TRUE)
        print(args)

        lengths_p_forward<-scan(as.character(args[1]), what = numeric())
        lengths_p_reverse<-scan(as.character(args[2]), what = numeric())
        lengths_q_forward<-scan(as.character(args[3]), what = numeric())
        lengths_q_reverse<-scan(as.character(args[4]), what = numeric())

        median_telomere_length_lengths_p_forward<-median(lengths_p_forward)
        median_telomere_length_lengths_p_reverse<-median(lengths_p_reverse)
        median_telomere_length_lengths_q_forward<-median(lengths_q_forward)
        median_telomere_length_lengths_q_reverse<-median(lengths_q_reverse)

        #percentile_90_telomere_length_starts <- quantile(lengths_starts, 0.9)
        #percentile_90_telomere_length_ends <- quantile(lengths_ends, 0.9)

        output_df <- data.frame(
          median_telomere_length_p_forward = median_telomere_length_lengths_p_forward, 
          median_telomere_length_p_reverse = median_telomere_length_lengths_p_reverse, 
          median_telomere_length_q_forward = median_telomere_length_lengths_q_forward, 
          median_telomere_length_q_reverse = median_telomere_length_lengths_q_reverse
        )

        write.table(output_df,file="telomeric.summary.txt",quote=FALSE,col.names=TRUE, row.names=FALSE)

        Rscript

    >>>

    output {
        File summary="telomeric.summary.txt"
    }

    runtime {
        docker: "rocker/tidyverse"
    }

}
