version 1.0


workflow findAssemblyBreakpoints{
    input {
        File assembly
        String assembly_name=basename(sub(sub(sub(assembly, "\\.gz$", ""), "\\.fasta$", ""), "\\.fa$", "")) #remove the file extension
        File reference
        String reference_name=basename(sub(sub(sub(assembly, "\\.gz$", ""), "\\.fasta$", ""), "\\.fa$", "")) #remove the file extension
        File annotationBed
        File annotationSD
        File annotationCENSAT
        Int telomericMinLength = 400
        Int flankLength = 1000
        Int threadCount = 32
        Int preemptible = 1
    }

    call formatAssembly {
        input:
            assembly=assembly,
            assembly_name=assembly_name,
            preemptible=preemptible
    }

    call generateAssemblyEdges {
        input:
            assembly=formatAssembly.formattedAssembly,
            assembly_name=assembly_name,
            preemptible=preemptible
    }

    call runBioawk {
        input:
            assembly=assembly,
            assembly_name=assembly_name,
            edges=generateAssemblyEdges.edges,
            preemptible=preemptible
    }
    
    call combineIntoFasta {
        input:
            assembly_name=assembly_name,
            headers=runBioawk.headers,
            edges=generateAssemblyEdges.edges,
            preemptible=preemptible
    }
    
    call runNCRF {
        input:
            edgesFasta=combineIntoFasta.edgesFasta,
            assembly_name=assembly_name,
            telomericMinLength=telomericMinLength,
            preemptible=preemptible
    }


    call runMashMap {
        input:
            assembly=assembly,
            assembly_name=assembly_name,
            reference=reference,
            threadCount=threadCount,
            preemptible=preemptible
    }

    call assessCompletness {
        input:
            assembly_name=assembly_name,
            telomericEnds=runNCRF.ends,
            lengths=runBioawk.lengths,
            unknown=runBioawk.unknown,
            mashmap=runMashMap.mashmap,
            preemptible=preemptible
    }

    call createGenomeFile {
        input:
            reference=reference,
            reference_name=reference_name,
            preemptible=preemptible
    }

    call mapToCHM13 {
        input:
            assembly=formatAssembly.formattedAssembly,
            assembly_name=assembly_name,
            reference=reference,
            preemptible=preemptible
    }

    call evaluate {
        input:
            mashmap_tmp=mapToCHM13.mashmap_tmp,
            scaffolds=assessCompletness.scaffolds,
            assembly_name=assembly_name,
            preemptible=preemptible
    }

    call createBedFiles {
        input:
            assembly_name=assembly_name,
            genomeFile=createGenomeFile.genomeFile,
            assemblyBed=evaluate.assemblyBed,
            flankLength=flankLength,
            preemptible=preemptible
    }

    call filterFlanks {
        input:
            mashmap=evaluate.mappedAssembly,
            flanksBedStart=createBedFiles.flanksBedStart,
            flanksBedEnd=createBedFiles.flanksBedEnd,
            telomericEnds=runNCRF.ends,
            assembly_name=assembly_name,
            flankLength=flankLength,
            preemptible=preemptible
    }

    call intersectBed {
        input:
            assembly_name=assembly_name,
            assemblyBed=filterFlanks.filteredFlanksBed,
            annotationBed=annotationBed,
            annotationSD=annotationSD,
            annotationCENSAT=annotationCENSAT,
            preemptible=preemptible
    }

    call formatBreakAnnotation {
        input:
            assembly_name=assembly_name,
            breakAnnotation_region=intersectBed.breakAnnotation_region,
            breakAnnotation_SD=intersectBed.breakAnnotation_SD,
            breakAnnotation_CENSAT=intersectBed.breakAnnotation_CENSAT,
            preemptible=preemptible
    }


    output {
        File T2Tcontigs = assessCompletness.contigs
        File T2Tscaffolds = assessCompletness.scaffolds
        File breakAnnotation_region = formatBreakAnnotation.formattedBreakAnnotation_region
        File breakAnnotation_SD = formatBreakAnnotation.formattedBreakAnnotation_SD
        File breakAnnotation_CENSAT = formatBreakAnnotation.formattedBreakAnnotation_CENSAT
        File assembly_CHM13 = evaluate.mappedAssembly
        File filteredFlanksBed = filterFlanks.filteredFlanksBed
    }

    parameter_meta {
        reference: " Reference which to use in order to derive chromosome names"
        assembly: " Assembly for which breakpoints should be evaluated"
        annotationBed: " Annotation in the bed format that is used to describe the assembly breaks"
    }
    meta {
        author: "Monika Cechova"
        email: "mcechova@ucsc.edu"
    }
}

task formatAssembly {
    input{
        File assembly
        String assembly_name
        Int memSizeGB = 32
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        #unify line endings and filter out sequences shorter than 100kb
        seqtk seq -l0 -L 200000 ~{assembly} >~{assembly_name}.formatted.fa

    >>>

    output {
        File formattedAssembly="~{assembly_name}.formatted.fa"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "quay.io/biocontainers/seqtk:1.3--hed695b0_2"
    }
}

task generateAssemblyEdges {
    input{
        File assembly
        String assembly_name
        Int memSizeGB = 32
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        bases=1000; 
        #extract first and last bases from each fasta entry
        
        #extract bases from the start and end of each entry in fasta
        cat ~{assembly} | gawk -v len="$bases" -F "" '{if (NR % 2 == 0) {for (i=1; i<=NF; i++) {if (i <= len || (NF-i) < len) {printf $(i)}}; printf "\n"}}' >~{assembly_name}.edges.tmp.txt
        
        #split start and end into two rows
        cat ~{assembly_name}.edges.tmp.txt | sed -e "s/.\{$bases\}/&\n/g" | sed '/^$/d' >~{assembly_name}.edges.txt
        #rm ~{assembly_name}.edges.tmp.txt

        
    >>>

    output {
        File edges="${assembly_name}.edges.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "quay.io/biocontainers/gawk:5.1.0"
    }
}

task runBioawk {
    input{
        File assembly
        File edges
        String assembly_name
        Int memSizeGB = 32
        Int preemptible

    }

    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        #create file with the number of unknown bases
        bioawk -c fastx '{n=gsub(/N/, "", $seq); print $name "\t" n}' ~{assembly} >~{assembly_name}.unknown.txt

        #create headers with the contig names and annotating starts and ends
        bioawk -c fastx '{print ">"$name"_start";print ">"$name"_end";}' ~{assembly} >~{assembly_name}.headers.txt

        #create file with sequence lengths
        bioawk -c fastx '{ print $name, length($seq) }' < ~{assembly} >~{assembly_name}.lengths.txt


        
    >>>

    output {
        File unknown = "${assembly_name}.unknown.txt"
        File headers = "${assembly_name}.headers.txt"
        File lengths = "${assembly_name}.lengths.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "quay.io/biocontainers/bioawk:1.0--hed695b0_5"
    }
}

task combineIntoFasta {
    input{
        File headers
        File edges
        String assembly_name
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        #combine extracted sequences with headers to create .fasta file
        paste -d '\n' ~{headers} ~{edges} > ~{assembly_name}.combined.fasta
        
    >>>

    output {
        File edgesFasta = "${assembly_name}.combined.fasta"
    }

    runtime {
        preemptible : preemptible
        docker: "ubuntu:18.04"
    }

}

task runNCRF {
    input{
        File edgesFasta
        String assembly_name
        Int telomericMinLength
        Int memSizeGB = 32
        Int preemptible
    }
    command <<<

    #handle potential errors and quit early
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    #search this fasta for the telomeric repeat
    cat ~{edgesFasta} | NCRF --stats=events TTAGGG --minlength=~{telomericMinLength} | ncrf_summary >~{assembly_name}.telomeric.summary.txt
    cat ~{assembly_name}.telomeric.summary.txt | tail -n +2 | cut -f3 | sort | uniq >~{assembly_name}.telomeric.ends.txt
        
    >>>

    output {
        File ends = "${assembly_name}.telomeric.ends.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "quay.io/biocontainers/ncrf:1.01.02--hec16e2b_3"
    }

}

task runMashMap {
    input{
        File assembly
        String assembly_name
        File reference
        Int memSizeGB = 32
        Int threadCount
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mashmap --threads ~{threadCount} --perc_identity 95 --noSplit -r ~{reference} -q ~{assembly} -o ~{assembly_name}.mashmap.txt
        
    >>>

    output {
        File mashmap = "${assembly_name}.mashmap.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        preemptible : preemptible
        docker: "quay.io/biocontainers/mashmap:2.0--h543ed2d_4"
    }

}

task assessCompletness {
    input{
        String assembly_name
        File telomericEnds
        File lengths
        File unknown
        File mashmap
        Int memSizeGB = 32
        Int preemptible
    }
    command <<<

    #handle potential errors and quit early
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    #CREATE A SUMMARY FILE
    cat ~{telomericEnds} | while read line; do contig=`echo $line | sed 's/_.*//'`; egrep --color ${contig} ~{lengths} | tr -d '\n'; echo -ne '\t'; egrep --color ${contig} ~{unknown} | tr -d '\n'; echo -ne '\t'; egrep --color ${contig} ~{mashmap} | tr -s "[:blank:]" "\t" ; done | cut -f1,2,4,10 | sort | uniq -c | sort -rgk1 | tr -s "[:blank:]" "\t" | sed 's/^\t//' >~{assembly_name}.SUMMARY.txt

    #write headers
    echo -e "name\tlength\tNs\tchromosome" >~{assembly_name}.T2T.scaffolds.txt
    echo -e "name\tlength\tNs\tchromosome" >~{assembly_name}.T2T.contigs.txt

    #check if there are telomeres on both sides
    if grep -q "^2" ~{assembly_name}.SUMMARY.txt; then
        #write T2T scaffolds, the number of Ns must be 0
        cat ~{assembly_name}.SUMMARY.txt | grep "^2" | awk '{if ($4>=0) print;}'  | cut -f2- | sort -k4,4 -V -s >>~{assembly_name}.T2T.scaffolds.txt || true

        #write T2T contigs, the number of Ns is greater than 0
        cat ~{assembly_name}.SUMMARY.txt | grep "^2" | awk '{if ($4==0) print;}'  | cut -f2- | sort -k4,4 -V -s >>~{assembly_name}.T2T.contigs.txt || true
    else
        echo "None of the contigs/scaffolds contained telomeres on both ends."
    fi
    
    echo "Done."    
    >>>

    output {
        File scaffolds = "${assembly_name}.T2T.scaffolds.txt"
        File contigs = "${assembly_name}.T2T.contigs.txt"
    }
    
    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "ubuntu:18.04"
    }
}

task createGenomeFile {
    input{
        File reference
        String reference_name
        Int memSizeGB = 32
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        bioawk -c fastx '{ print $name, length($seq) }' < ~{reference} > ~{reference_name}.genomeFile.txt

    >>>

    output {
        File genomeFile="${reference_name}.genomeFile.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "quay.io/biocontainers/bioawk:1.0--hed695b0_5"
    }
}



task mapToCHM13 {
    input{
        File assembly
        String assembly_name
        File reference
        Int memSizeGB = 32
        Int threadCount = 100
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mashmap --threads ~{threadCount} --perc_identity 95 --filter_mode one-to-one --noSplit -r ~{reference} -q ~{assembly} -o ~{assembly_name}.mashmap.tmp

    >>>

    output {
        File mashmap_tmp = "${assembly_name}.mashmap.tmp"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        preemptible : preemptible
        docker: "quay.io/biocontainers/mashmap:2.0--h543ed2d_4"
    }

}

task evaluate {
    input{
        File mashmap_tmp
        File scaffolds
        String assembly_name
        Int memSizeGB = 1
        Int threadCount = 11
        Int preemptible
    }
    command <<<

    #handle potential errors and quit early
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    #exclude those rows that are too close to complete


    #EXCLUDE COMPLETE CHROMOSOMES
    cat ~{mashmap_tmp} | sort -k6,6 -k8,8 -V -s >~{assembly_name}.mashmap.txt       #coordinates in the reference/CHM13
    cat ~{scaffolds} | cut -f4 >~{assembly_name}.complete_chromosomes.lst           #list of complete chromosomes that should be exluded from the breakpoint analysis
    
    #create bed file from the mashmap alignment
    awk '{print $6 "\t" $8 "\t" $9 "\t" $1}' ~{mashmap_tmp} | sort -k1,1 -k2,2n -V -s >~{assembly_name}.tmp.bed     #create sorted bed file
    grep -v -f ~{assembly_name}.complete_chromosomes.lst ~{assembly_name}.tmp.bed >~{assembly_name}.bed || true         #exclude complete chromosomes/scaffolds

   
    >>>

    output {
        File mappedAssembly = "${assembly_name}.mashmap.txt"
        File assemblyBed = "${assembly_name}.bed"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "ubuntu:18.04"
    }
}

task createBedFiles {
    input{
        String assembly_name
        File genomeFile
        File assemblyBed
        Int flankLength
        Int memSizeGB = 1
        Int threadCount = 11
        Int preemptible
    }
    command <<<

    #handle potential errors and quit early
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    #1kb on each side of the breakpoint
    bedtools flank -i ~{assemblyBed} -g ~{genomeFile} -l ~{flankLength} -r 0 >~{assembly_name}.flanks.start.bed
    bedtools flank -i ~{assemblyBed} -g ~{genomeFile} -l 0 -r ~{flankLength} >~{assembly_name}.flanks.end.bed
    
    >>>

    output {
        File flanksBedStart = "${assembly_name}.flanks.start.bed"
        File flanksBedEnd = "${assembly_name}.flanks.end.bed"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6"
    }
}

task filterFlanks {
    input{
        File mashmap
        File flanksBedStart
        File flanksBedEnd
        File telomericEnds
        String assembly_name
        Int flankLength = 1000
        Int memSizeGB = 1
        Int threadCount = 11
        Int preemptible
    }
    command <<<

    #handle potential errors and quit early
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    #EXCLUDE BREAKS THAT ARE NOT REAL BREAKS, AS THEY ARE AT THE CHROMOSOME END AND CONTAIN TELOMERIC REPEAT
    awk '{if ($8<(1000000)) print;}' ~{mashmap} | cut -d' ' -f1 >~{assembly_name}.chromosome.starts.txt
    awk '{if (($7-$9)<(1000000)) print;}' ~{mashmap} | cut -d' ' -f1 >~{assembly_name}.chromosome.ends.txt

    while read -r contig; do grep "${contig}_start" ~{telomericEnds} | sed 's/_start//g' >>~{assembly_name}.chromosome.starts.toExclude.txt || true; done < ~{assembly_name}.chromosome.starts.txt
    
    while read -r contig; do grep "${contig}_end" ~{telomericEnds} | sed 's/_end//g' >>~{assembly_name}.chromosome.ends.toExclude.txt || true; done < ~{assembly_name}.chromosome.ends.txt

    #exclude starts that have telomeric repeat
    grep -v -f ~{assembly_name}.chromosome.starts.toExclude.txt ~{flanksBedStart}  >~{assembly_name}.NOstarts.bed || true

    #exclude ends that have telomeric repeat
    grep -v -f ~{assembly_name}.chromosome.ends.toExclude.txt ~{flanksBedEnd}  >~{assembly_name}.NOends.bed || true


    cat ~{assembly_name}.NOstarts.bed ~{assembly_name}.NOends.bed >~{assembly_name}.merged.bed
    awk '{if (($3-$2)==~{flankLength}) print;}' ~{assembly_name}.merged.bed >~{assembly_name}.filteredFlanks.bed
    
    >>>

    output {
        File filteredFlanksBed = "${assembly_name}.filteredFlanks.bed"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "ubuntu:18.04"
    }
}

task intersectBed {
    input{
        String assembly_name
        File assemblyBed
        File annotationBed
        File annotationSD
        File annotationCENSAT
        Int memSizeGB = 1
        Int threadCount = 11
        Int preemptible
    }
    command <<<

    #handle potential errors and quit early
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    #intersect flanks with the annotation of genomic regions
    echo -e "~{assembly_name}""\t""annotation" >~{assembly_name}.region.breakAnnotation.txt
    bedtools intersect -a ~{assemblyBed} -b ~{annotationBed} -loj | cut -f8 | sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' >>~{assembly_name}.region.breakAnnotation.txt
    
    echo -e "~{assembly_name}""\t""annotation" >~{assembly_name}.SD.breakAnnotation.txt
    bedtools intersect -a ~{assemblyBed} -b ~{annotationSD} -loj | cut -f8 | sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' >>~{assembly_name}.SD.breakAnnotation.txt
    
    echo -e "~{assembly_name}""\t""annotation" >~{assembly_name}.CENSAT.breakAnnotation.txt
    bedtools intersect -a ~{assemblyBed} -b ~{annotationCENSAT} -loj | cut -f8 | sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' >>~{assembly_name}.CENSAT.breakAnnotation.txt
    
    >>>

    output {
        File breakAnnotation_region = "${assembly_name}.region.breakAnnotation.txt"
        File breakAnnotation_SD = "${assembly_name}.SD.breakAnnotation.txt"
        File breakAnnotation_CENSAT = "${assembly_name}.CENSAT.breakAnnotation.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6"
    }
}

task formatBreakAnnotation {
    input{
        String assembly_name
        File breakAnnotation_region
        File breakAnnotation_SD
        File breakAnnotation_CENSAT
        Int memSizeGB = 1
        Int threadCount = 11
        Int preemptible
    }
    command <<<

    #handle potential errors and quit early
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    #intersect flanks with the annotation of genomic regions
    cat ~{breakAnnotation_region} | datamash reverse | datamash transpose --filler NA >>~{assembly_name}.formattedBreakAnnotation_region.txt
    cat ~{breakAnnotation_SD} | datamash reverse | datamash transpose --filler NA >>~{assembly_name}.formattedBreakAnnotation_SD.txt
    cat ~{breakAnnotation_CENSAT} | datamash reverse | datamash transpose --filler NA >>~{assembly_name}.formattedBreakAnnotation_CENSAT.txt
    
    >>>

    output {
        File formattedBreakAnnotation_region = "${assembly_name}.formattedBreakAnnotation_region.txt"
        File formattedBreakAnnotation_SD = "${assembly_name}.formattedBreakAnnotation_SD.txt"
        File formattedBreakAnnotation_CENSAT = "${assembly_name}.formattedBreakAnnotation_CENSAT.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "quay.io/biocontainers/datamash:1.1.0--0"
    }
}




