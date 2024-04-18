version 1.0


workflow evaluateHumanAssembly{
    input {
        File assembly
        String assembly_name=basename(sub(sub(sub(assembly, "\\.gz$", ""), "\\.fasta$", ""), "\\.fa$", "")) #remove the file extension
        Int telomericMinLength = 400
        Int telomericMaxDistance = 1000
        Int flankLength = 1000
        Int minContigLength = 200000
        Int minIdentity = 95
        Int threadCount = 32
        Int preemptible = 1
    }

    call getGenomeFile {
        input:
            preemptible=preemptible
    }

    call createGenomeFile {
        input:
            reference=getGenomeFile.reference,
            preemptible=preemptible
    }

    call formatAssembly {
        input:
            assembly=assembly,
            assembly_name=assembly_name,
            minContigLength=minContigLength,
            preemptible=preemptible
    }

    call trialMashmap {
        input:
            assembly=formatAssembly.formattedAssembly,
            assembly_name=assembly_name,
            reference=getGenomeFile.reference,
            minIdentity=minIdentity,
            preemptible=preemptible
    }

    call unifyAssembly {
        input:
            assembly=formatAssembly.formattedAssembly,
            assembly_name=assembly_name,
            contigsToBeReverseComplemented=trialMashmap.contigsToBeReverseComplemented,
            contigsToBeKeptAsIs=trialMashmap.contigsToBeKeptAsIs,
            preemptible=preemptible
    }

    call mapToCHM13 {
        input:
            assembly=unifyAssembly.unifiedAssembly,
            assembly_name=assembly_name,
            reference=getGenomeFile.reference,
            minIdentity=minIdentity,
            preemptible=preemptible
    }

    call runBioawk {
        input:
            assembly=unifyAssembly.unifiedAssembly,
            assembly_name=assembly_name,
            preemptible=preemptible
    }
    
    
    call runSeqtk {
        input:
            assembly=unifyAssembly.unifiedAssembly,
            assembly_name=assembly_name,
            telomericMinLength=telomericMinLength,
            telomericMaxDistance=telomericMaxDistance,
            preemptible=preemptible
    }

    call createAssemblyStatistics {
        input:
            mashmap=evaluate.mappedAssembly,
            preemptible=preemptible
    }

    call assessCompletness {
        input:
            assembly_name=assembly_name,
            telomericEnds=runSeqtk.ends,
            lengths=runBioawk.lengths,
            unknown=runBioawk.unknown,
            mashmap=mapToCHM13.mashmap,
            preemptible=preemptible
    }

    call evaluate {
        input:
            mashmap=mapToCHM13.mashmap,
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
            telomericEnds=runSeqtk.ends,
            assembly_name=assembly_name,
            flankLength=flankLength,
            preemptible=preemptible
    }

    call intersectBed {
        input:
            assembly_name=assembly_name,
            assemblyBed=filterFlanks.filteredFlanksBed,
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
        File unifiedAssembly = unifyAssembly.unifiedAssembly
        File chromosome_names_from_CHM13 = evaluate.mappedAssembly
        File filteredFlanksBed = filterFlanks.filteredFlanksBed
        File bed_region = intersectBed.bed_region
        File bed_SD = intersectBed.bed_SD
        File bed_CENSAT = intersectBed.bed_CENSAT
        File assemblyStatistics = createAssemblyStatistics.assemblyStatistics
    }

    parameter_meta {
        assembly: " Assembly for which breakpoints should be evaluated"
        telomericMinLength: "The minimal length for a telomeric repeat to be considered a telomere"
        telomericMaxDistance: "The maximal distance from either end of a chromosome for a telomeric repeat to be considered a telomere"
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
        Int minContigLength
        Int memSizeGB = 32
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        #unify line endings and filter out sequences shorter than 200kb
        seqtk seq -l0 -L ~{minContigLength} ~{assembly} >~{assembly_name}.formatted.fa

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

task trialMashmap {
    input{
        File assembly
        String assembly_name
        File reference
        Int minIdentity
        Int memSizeGB = 32
        Int threadCount = 32
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mashmap --threads ~{threadCount} --perc_identity ~{minIdentity} --noSplit -r ~{reference} -q ~{assembly} -o ~{assembly_name}.mashmap.tmp

        cat ~{assembly_name}.mashmap.tmp | awk '{if ($5=="-") print $1}' | sed 's/>//g' >~{assembly_name}.contigsToBeReverseComplemented.lst
        cat ~{assembly_name}.mashmap.tmp | awk '{if ($5=="+") print $1}' | sed 's/>//g' >~{assembly_name}.contigsToBeKeptAsIs.lst

    >>>

    output {
        File contigsToBeReverseComplemented = "${assembly_name}.contigsToBeReverseComplemented.lst"
        File contigsToBeKeptAsIs = "${assembly_name}.contigsToBeKeptAsIs.lst"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        preemptible : preemptible
        docker: "quay.io/biocontainers/mashmap:2.0--h543ed2d_4"
    }
}

task unifyAssembly {
    input{
        File assembly
        String assembly_name
        File contigsToBeReverseComplemented
        File contigsToBeKeptAsIs
        Int memSizeGB = 32
        Int threadCount = 32
        Int diskSizeGB = 25
        Int preemptible
    }

    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        seqtk subseq ~{assembly} ~{contigsToBeReverseComplemented} >contigsToBeReverseComplemented.tmp
        seqtk seq -r contigsToBeReverseComplemented.tmp >contigsToBeReverseComplemented.fa
        rm contigsToBeReverseComplemented.tmp
        seqtk subseq ~{assembly} ~{contigsToBeKeptAsIs} >contigsToBeKeptAsIs.fa
        
        cat contigsToBeReverseComplemented.fa contigsToBeKeptAsIs.fa >~{assembly_name}.unifiedAssembly.fa

        #remove unnecessary files
        rm -f contigsToBeReverseComplemented.fa
        rm -f contigsToBeKeptAsIs.fa

    >>>

    output {
        File unifiedAssembly = "${assembly_name}.unifiedAssembly.fa"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible : preemptible
        docker: "quay.io/biocontainers/seqtk:1.3--hed695b0_2"
    }
}

task mapToCHM13 {
    input{
        File assembly
        String assembly_name
        File reference
        Int minIdentity
        Int memSizeGB = 32
        Int threadCount = 32
        Int preemptible
    }

    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mashmap --threads ~{threadCount} --perc_identity ~{minIdentity} --noSplit -r ~{reference} -q ~{assembly} -o ~{assembly_name}.mashmap.tmp

    >>>

    output {
        File mashmap = "${assembly_name}.mashmap.tmp"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        preemptible : preemptible 
        docker: "quay.io/biocontainers/mashmap:2.0--h543ed2d_4"
    }
}

task runBioawk {
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


task runSeqtk {
    input{
        File assembly
        String assembly_name
        Int telomericMinLength
        Int telomericMaxDistance
        Int memSizeGB = 32
        Int preemptible
    }
    
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        #search the assembly for the telomeric repeat
        seqtk telo ~{assembly} > ~{assembly_name}.telo.bed 2> ~{assembly_name}.telo.count

        #check telomeres at the beginning of chromosomes
        awk -v maxdist=~{telomericMaxDistance} -v minlen=~{telomericMinLength} '{
                if (($2 <= maxdist) && ($3 - $2) >= minlen) {
                    print $1"_start";
                }
            }' ~{assembly_name}.telo.bed >> ~{assembly_name}.telomeric.ends.txt

        #check telomeres at the ends of chromosomes
        awk -v maxdist=~{telomericMaxDistance} -v minlen=~{telomericMinLength} '{
                if (($3 >= ($4-maxdist)) && ($3 - $2) >= minlen) {
                    print $1"_end";
                }
            }' ~{assembly_name}.telo.bed >> ~{assembly_name}.telomeric.ends.txt

        #sort the final file
        sort ~{assembly_name}.telomeric.ends.txt > sorted.~{assembly_name}.telomeric.ends.txt && mv sorted.~{assembly_name}.telomeric.ends.txt ~{assembly_name}.telomeric.ends.txt
        
    >>>

    output {
        File ends = "${assembly_name}.telomeric.ends.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "quay.io/biocontainers/seqtk:1.4--he4a0461_2"
    }

}

task createAssemblyStatistics {
    input{
        File mashmap
        Int preemptible
    }

    command <<<
        R --no-save --args ~{sep=" " mashmap} <<Rscript
        args <- commandArgs(trailingOnly = TRUE)
        filename<-args[1]

        mashmap<-as.data.frame(read.table(filename))
        colnames(mashmap)<-c("contig","contig_length","cstart","cend","strand","chr","chr_length","chstart","chend","identity")
        coef<-1000000

        library("dplyr")

        #generate useful statistics about the assembly, e.g. how many Mbs of sequence were assembled
        filtered_mashmap<-mashmap %>% 
            mutate(chm13_length = chend - chstart) %>%
            mutate(assembly_length = cend - cstart) %>%
            group_by(chr) %>% 
            summarise(reference_alignment_Mb = sum(chm13_length)/coef, 
            assembly_alignment_Mb = sum(assembly_length)/coef, 
            chromosome_length_Mb=max(chr_length)/coef,
            largest_assembled_sequence_Mb=max(assembly_length)/coef)%>%
            arrange(chr,.by_group=TRUE)
        write.table(filtered_mashmap,file="assembly.statistics.tmp",quote=FALSE,row.names=FALSE,sep="\t")
        Rscript

        cat assembly.statistics.tmp | sort -k1,1 -V -s >assembly.statistics.txt

    >>>

    output {
        File assemblyStatistics="assembly.statistics.txt"
    }

    runtime {
        docker: "rocker/tidyverse"
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

task getGenomeFile {
    input{
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        cp /custom_annotations/chm13v2.0.fa.gz chm13v2.0.fa.gz

    >>>

    output {
        File reference="chm13v2.0.fa.gz"
    }

    runtime {
        preemptible : preemptible
        docker: "biomonika/custom_annot:latest"
    }
}

task createGenomeFile {
    input{
        File reference
        Int memSizeGB = 32
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        bioawk -c fastx '{ print $name, length($seq) }' < ~{reference} > genomeFile.txt

    >>>

    output {
        File genomeFile="genomeFile.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "quay.io/biocontainers/bioawk:1.0--hed695b0_5"
    }
}

task evaluate {
    input{
        File mashmap
        File scaffolds
        String assembly_name
        Int memSizeGB = 32
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
        cat ~{mashmap} | sort -k6,6 -k8,8 -V -s >~{assembly_name}.mashmap.txt       #coordinates in the reference/CHM13
        cat ~{scaffolds} | cut -f4 >~{assembly_name}.complete_chromosomes.lst           #list of complete chromosomes that should be exluded from the breakpoint analysis
        
        #create bed file from the mashmap alignment
        awk '{print $6 "\t" $8 "\t" $9 "\t" $1}' ~{mashmap} | sort -k1,1 -k2,2n -V -s >~{assembly_name}.tmp.bed     #create sorted bed file
        grep -w -v -f ~{assembly_name}.complete_chromosomes.lst ~{assembly_name}.tmp.bed >~{assembly_name}.bed || true         #exclude complete chromosomes/scaffolds

   
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
        Int memSizeGB = 4
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
        Int memSizeGB = 32
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
        awk '{if ($8<(1000000)) print;}' ~{mashmap} | cut -d' ' -f1 | sort >~{assembly_name}.chromosome.starts.txt
        awk '{if (($7-$9)<(1000000)) print;}' ~{mashmap} | cut -d' ' -f1 | sort >~{assembly_name}.chromosome.ends.txt

        #if both the chromosome start _and_ has telomeric repeat, then exclude
        cat ~{telomericEnds} | grep "start" | sed 's/_start//g' | sort >starts_with_telomeric_repeat.txt
        comm -12 ~{assembly_name}.chromosome.starts.txt starts_with_telomeric_repeat.txt >~{assembly_name}.chromosome.starts.toExclude.txt

        #if both the chromosome end _and_ has telomeric repeat, then exclude
        cat ~{telomericEnds} | grep "end" | sed 's/_end//g' | sort >ends_with_telomeric_repeat.txt
        comm -12 ~{assembly_name}.chromosome.ends.txt ends_with_telomeric_repeat.txt >~{assembly_name}.chromosome.ends.toExclude.txt

        #exclude starts that have telomeric repeat
        grep -w -v -f ~{assembly_name}.chromosome.starts.toExclude.txt ~{flanksBedStart}  >~{assembly_name}.NOstarts.bed || true

        #exclude ends that have telomeric repeat
        grep -w -v -f ~{assembly_name}.chromosome.ends.toExclude.txt ~{flanksBedEnd}  >~{assembly_name}.NOends.bed || true


        cat ~{assembly_name}.NOstarts.bed ~{assembly_name}.NOends.bed >~{assembly_name}.filteredFlanks.bed

        #if the filtering below is not performed, some of the flanks will be potentially as short as 1bp 
        #awk '{if (($3-$2)==~{flankLength}) print;}' ~{assembly_name}.merged.bed >~{assembly_name}.filteredFlanks.bed
    
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
        Int memSizeGB = 4
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
        echo -e "~{assembly_name}""\t""annotation" >~{assembly_name}.SD.breakAnnotation.txt
        echo -e "~{assembly_name}""\t""annotation" >~{assembly_name}.CENSAT.breakAnnotation.txt

        cp /custom_annotations/chm13v2.0_GenomeFeature_v1.0.bed chm13v2.0_GenomeFeature_v1.0.bed
        cp /custom_annotations/chm13v2.0_SD.flattened.bed chm13v2.0_SD.flattened.bed
        cp /custom_annotations/chm13v2.0_censat_v2.0.bed chm13v2.0_censat_v2.0.bed

        annotationBed="chm13v2.0_GenomeFeature_v1.0.bed"
        annotationSD="chm13v2.0_SD.flattened.bed"
        annotationCENSAT="chm13v2.0_censat_v2.0.bed"

        bedtools intersect -a ~{assemblyBed} -b ${annotationBed} -loj | cut -f8 | sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' >>~{assembly_name}.region.breakAnnotation.txt
        bedtools intersect -a ~{assemblyBed} -b ${annotationSD} -loj | cut -f8 | sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' >>~{assembly_name}.SD.breakAnnotation.txt
        bedtools intersect -a ~{assemblyBed} -b ${annotationCENSAT} -loj | cut -f8 | sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' >>~{assembly_name}.CENSAT.breakAnnotation.txt

        #GENERATE BED FILES FOR PLOTTING
        echo -e "#chr""\t""start""\t""end""\t""annotation" >~{assembly_name}.region.breakAnnotation.bed
        echo -e "#chr""\t""start""\t""end""\t""annotation" >~{assembly_name}.region.breakAnnotation.bed
        echo -e "#chr""\t""start""\t""end""\t""annotation" >~{assembly_name}.region.breakAnnotation.bed
        
        bedtools intersect -a ~{assemblyBed} -b ${annotationBed} -loj | awk '{print $1 "\t" $2 "\t" $3 "\t" $8}' | sort -k1,1 -k2,2n -V -s >>~{assembly_name}.region.breakAnnotation.bed
        bedtools intersect -a ~{assemblyBed} -b ${annotationSD} -loj | awk '{print $1 "\t" $2 "\t" $3 "\t" $8}' | sort -k1,1 -k2,2n -V -s >>~{assembly_name}.SD.breakAnnotation.bed
        bedtools intersect -a ~{assemblyBed} -b ${annotationCENSAT} -loj | awk '{print $1 "\t" $2 "\t" $3 "\t" $8}' | sort -k1,1 -k2,2n -V -s >>~{assembly_name}.CENSAT.breakAnnotation.bed

    >>>

    output {
        File breakAnnotation_region = "${assembly_name}.region.breakAnnotation.txt"
        File breakAnnotation_SD = "${assembly_name}.SD.breakAnnotation.txt"
        File breakAnnotation_CENSAT = "${assembly_name}.CENSAT.breakAnnotation.txt"
        File bed_region = "${assembly_name}.region.breakAnnotation.bed"
        File bed_SD = "${assembly_name}.SD.breakAnnotation.bed"
        File bed_CENSAT = "${assembly_name}.CENSAT.breakAnnotation.bed"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        docker: "biomonika/custom_annot:latest"
    }
}

task formatBreakAnnotation {
    input{
        String assembly_name
        File breakAnnotation_region
        File breakAnnotation_SD
        File breakAnnotation_CENSAT
        Int memSizeGB = 4
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

