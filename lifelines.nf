#!/usr/bin/env nextflow

import Helper
import CollectInitialMetadata

// Pipeline version
if (workflow.commitId){
    version = "0.1 $workflow.revision"
} else {
    version = "0.1 (local version)"
}

params.help = false
if (params.help){
    Help.print_help(params)
    exit 0
}

def infoMap = [:]
if (params.containsKey("fastq")){
    infoMap.put("fastq", file(params.fastq).size())
}
if (params.containsKey("fasta")){
    if (file(params.fasta) instanceof LinkedList){
        infoMap.put("fasta", file(params.fasta).size())
    } else {
        infoMap.put("fasta", 1) 
    }
}
if (params.containsKey("accessions")){
    // checks if params.accessions is different from null
    if (params.accessions) {
        BufferedReader reader = new BufferedReader(new FileReader(params.accessions));
        int lines = 0;
        while (reader.readLine() != null) lines++;
        reader.close();
        infoMap.put("accessions", lines)
    }
}

Help.start_info(infoMap, "$workflow.start", "$workflow.profile")
CollectInitialMetadata.print_metadata(workflow)
    

// Placeholder for main input channels
if (params.fastq instanceof Boolean){exit 1, "'fastq' must be a path pattern. Provide value:'$params.fastq'"}
if (!params.fastq){ exit 1, "'fastq' parameter missing"}
IN_fastq_raw = Channel.fromFilePairs(params.fastq).ifEmpty { exit 1, "No fastq files provided with pattern:'${params.fastq}'" }

// Placeholder for secondary input channels


// Placeholder for extra input channels


// Placeholder to fork the raw input channel

IN_fastq_raw.set{ integrity_coverage_in_1_0 }


IN_genome_size_1_1 = Channel.value(params.genomeSize_1_1)
    .map{it -> it.toString().isNumber() ? it : exit(1, "The genomeSize parameter must be a number or a float. Provided value: '${params.genomeSize__1_1}'")}

IN_min_coverage_1_1 = Channel.value(params.minCoverage_1_1)
    .map{it -> it.toString().isNumber() ? it : exit(1, "The minCoverage parameter must be a number or a float. Provided value: '${params.minCoverage__1_1}'")}

process integrity_coverage_1_1 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_1 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_1 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_1 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId integrity_coverage_1_1 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set sample_id, file(fastq_pair) from integrity_coverage_in_1_0
    val gsize from IN_genome_size_1_1
    val cov from IN_min_coverage_1_1
    // This channel is for the custom options of the integrity_coverage.py
    // script. See the script's documentation for more information.
    val opts from Channel.value('')

    output:
    set sample_id,
        file(fastq_pair),
        file('*_encoding'),
        file('*_phred'),
        file('*_coverage'),
        file('*_max_len') into MAIN_integrity_1_1
    file('*_report') optional true into LOG_report_coverage1_1_1
    set sample_id, val("1_1_integrity_coverage"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_integrity_coverage_1_1
set sample_id, val("integrity_coverage_1_1"), val("1_1"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_integrity_coverage_1_1
file ".versions"

    script:
    template "integrity_coverage.py"

}

// TRIAGE OF CORRUPTED SAMPLES
LOG_corrupted_1_1 = Channel.create()
MAIN_PreCoverageCheck_1_1 = Channel.create()
// Corrupted samples have the 2nd value with 'corrupt'
MAIN_integrity_1_1.choice(LOG_corrupted_1_1, MAIN_PreCoverageCheck_1_1) {
    a -> a[2].text == "corrupt" ? 0 : 1
}

// TRIAGE OF LOW COVERAGE SAMPLES
integrity_coverage_out_1_0 = Channel.create()
SIDE_phred_1_1 = Channel.create()
SIDE_max_len_1_1 = Channel.create()

MAIN_PreCoverageCheck_1_1
// Low coverage samples have the 4th value of the Channel with 'fail'
    .filter{ it[4].text != "fail" }
// For the channel to proceed with FastQ in 'sample_good' and the
// Phred scores for each sample in 'SIDE_phred'
    .separate(integrity_coverage_out_1_0, SIDE_phred_1_1, SIDE_max_len_1_1){
        a -> [ [a[0], a[1]], [a[0], a[3].text], [a[0], a[5].text]  ]
    }

/** REPORT_COVERAGE - PLUG-IN
This process will report the expected coverage for each non-corrupted sample
and write the results to 'reports/coverage/estimated_coverage_initial.csv'
*/
process report_coverage_1_1 {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/coverage_1_1/'

    input:
    file(report) from LOG_report_coverage1_1_1.filter{ it.text != "corrupt" }.collect()

    output:
    file 'estimated_coverage_initial.csv'

    """
    echo Sample,Estimated coverage,Status >> estimated_coverage_initial.csv
    cat $report >> estimated_coverage_initial.csv
    """
}

/** REPORT_CORRUPT - PLUG-IN
This process will report the corrupted samples and write the results to
'reports/corrupted/corrupted_samples.txt'
*/
process report_corrupt_1_1 {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/corrupted_1_1/'

    input:
    val sample_id from LOG_corrupted_1_1.collect{it[0]}

    output:
    file 'corrupted_samples.txt'

    """
    echo ${sample_id.join(",")} | tr "," "\n" >> corrupted_samples.txt
    """

}


SIDE_phred_1_1.set{ SIDE_phred_1_3 }


SIDE_max_len_1_1.set{ SIDE_max_len_3_5 }


IN_index_files_1_2 = Channel.value(params.refIndex_1_2)

clear = params.clearInput_1_2 ? "true" : "false"
checkpointClear_1_2 = Channel.value(clear)

process remove_host_1_2 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_2 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_2 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_2 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId remove_host_1_2 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir 'results/mapping/remove_host_1_2/', pattern: '*_bowtie2.log', mode: 'copy'

    input:
    set sample_id, file(fastq_pair) from integrity_coverage_out_1_0
    val bowtie2Index from IN_index_files_1_2
    val clear from checkpointClear_1_2

    output:
    set sample_id , file("${sample_id}*.headersRenamed_*.fq.gz") into remove_host_out_1_1
    set sample_id, file("*_bowtie2.log") into into_json_1_2
    set sample_id, val("1_2_remove_host"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_remove_host_1_2
set sample_id, val("remove_host_1_2"), val("1_2"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_remove_host_1_2
file ".versions"

    script:
    """
    {
        bowtie2 -x ${bowtie2Index} -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -p $task.cpus 1> ${sample_id}.bam 2> ${sample_id}_bowtie2.log

        samtools view -buh -f 12 -o ${sample_id}_samtools.bam -@ $task.cpus ${sample_id}.bam

        rm ${sample_id}.bam

        samtools fastq -1 ${sample_id}_unmapped_1.fq -2 ${sample_id}_unmapped_2.fq ${sample_id}_samtools.bam

        rm ${sample_id}_samtools.bam

        renamePE_samtoolsFASTQ.py -1 ${sample_id}_unmapped_1.fq -2 ${sample_id}_unmapped_2.fq

        gzip *.headersRenamed_*.fq
        rm *.fq

        if [ "$clear" = "true" ];
        then
            work_regex=".*/work/.{2}/.{30}/.*"
            file_source1=\$(readlink -f \$(pwd)/${fastq_pair[0]})
            file_source2=\$(readlink -f \$(pwd)/${fastq_pair[1]})
            if [[ "\$file_source1" =~ \$work_regex ]]; then
                rm \$file_source1 \$file_source2
            fi
        fi

    } || {
        echo fail > .status
    }
    """
}



process report_remove_host_1_2 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_2 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_2 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_2 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId remove_host_1_2 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    input:
    set sample_id, file(bowtie_log) from into_json_1_2

    output:
    set sample_id, val("1_2_report_remove_host"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_report_remove_host_1_2
set sample_id, val("report_remove_host_1_2"), val("1_2"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_report_remove_host_1_2
file ".versions"

    script:
    template "process_mapping.py"

}


// Check sliding window parameter
if ( params.trimSlidingWindow_1_3.toString().split(":").size() != 2 ){
    exit 1, "'trimSlidingWindow_1_3' parameter must contain two values separated by a ':'. Provided value: '${params.trimSlidingWindow_1_3}'"
}
if ( !params.trimLeading_1_3.toString().isNumber() ){
    exit 1, "'trimLeading_1_3' parameter must be a number. Provide value: '${params.trimLeading_1_3}'"
}
if ( !params.trimTrailing_1_3.toString().isNumber() ){
    exit 1, "'trimTrailing_1_3' parameter must be a number. Provide value: '${params.trimTrailing_1_3}'"
}
if ( !params.trimMinLength_1_3.toString().isNumber() ){
    exit 1, "'trimMinLength_1_3' parameter must be a number. Provide value: '${params.trimMinLength_1_3}'"
}

IN_trimmomatic_opts_1_3 = Channel.value([params.trimSlidingWindow_1_3,params.trimLeading_1_3,params.trimTrailing_1_3,params.trimMinLength_1_3])
IN_adapters_1_3 = Channel.value(params.adapters_1_3)

clear = params.clearInput_1_3 ? "true" : "false"
checkpointClear_1_3 = Channel.value(clear)

process fastqc_1_3 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_3 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId fastqc_trimmomatic_1_3 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir "reports/fastqc_1_3/", pattern: "*.html"

    input:
    set sample_id, file(fastq_pair) from remove_host_out_1_1
    val ad from Channel.value('None')

    output:
    set sample_id, file(fastq_pair), file('pair_1*'), file('pair_2*') into MAIN_fastqc_out_1_3
    file "*html"
    set sample_id, val("1_3_fastqc"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_fastqc_1_3
set sample_id, val("fastqc_1_3"), val("1_3"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_fastqc_1_3
file ".versions"

    script:
    template "fastqc.py"
}

/** FASTQC_REPORT - MAIN
This process will parse the result files from a FastQC analyses and output
the optimal_trim information for Trimmomatic
*/
process fastqc_report_1_3 {

    // Send POST request to platform
    
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_3 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId fastqc_trimmomatic_1_3 \"$params.platformSpecies\" false"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }
    

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/fastqc_1_3/run_1/', pattern: '*summary.txt', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out_1_3
    val opts from Channel.value("--ignore-tests")

    output:
    set sample_id, file(fastq_pair), 'optimal_trim', ".status" into _MAIN_fastqc_trim_1_3
    file '*_trim_report' into LOG_trim_1_3
    file "*_status_report" into LOG_fastqc_report_1_3
    file "${sample_id}_*_summary.txt" optional true
    set sample_id, val("1_3_fastqc_report"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_fastqc_report_1_3
set sample_id, val("fastqc_report_1_3"), val("1_3"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_fastqc_report_1_3
file ".versions"

    script:
    template "fastqc_report.py"

}

MAIN_fastqc_trim_1_3 = Channel.create()
_MAIN_fastqc_trim_1_3
        .filter{ it[3].text == "pass" }
        .map{ [it[0], it[1], file(it[2]).text] }
        .into(MAIN_fastqc_trim_1_3)


/** TRIM_REPORT - PLUG-IN
This will collect the optimal trim points assessed by the fastqc_report
process and write the results of all samples in a single csv file
*/
process trim_report_1_3 {

    publishDir 'reports/fastqc_1_3/', mode: 'copy'

    input:
    file trim from LOG_trim_1_3.collect()

    output:
    file "FastQC_trim_report.csv"

    """
    echo Sample,Trim begin, Trim end >> FastQC_trim_report.csv
    cat $trim >> FastQC_trim_report.csv
    """
}


process compile_fastqc_status_1_3 {

    publishDir 'reports/fastqc_1_3/', mode: 'copy'

    input:
    file rep from LOG_fastqc_report_1_3.collect()

    output:
    file 'FastQC_1run_report.csv'

    """
    echo Sample, Failed? >> FastQC_1run_report.csv
    cat $rep >> FastQC_1run_report.csv
    """

}


/** TRIMMOMATIC - MAIN
This process will execute trimmomatic. Currently, the main channel requires
information on the trim_range and phred score.
*/
process trimmomatic_1_3 {

    // Send POST request to platform
    
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_3 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId fastqc_trimmomatic_1_3 \"$params.platformSpecies\" false"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }
    

    tag { sample_id }
    publishDir "results/trimmomatic_1_3", pattern: "*.gz"

    input:
    set sample_id, file(fastq_pair), trim_range, phred from MAIN_fastqc_trim_1_3.join(SIDE_phred_1_3)
    val opts from IN_trimmomatic_opts_1_3
    val ad from IN_adapters_1_3
    val clear from checkpointClear_1_3

    output:
    set sample_id, "${sample_id}_*trim.fastq.gz" into _fastqc_trimmomatic_out_1_2
    file 'trimmomatic_report.csv'
    set sample_id, val("1_3_trimmomatic"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_trimmomatic_1_3
set sample_id, val("trimmomatic_1_3"), val("1_3"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_trimmomatic_1_3
file ".versions"

    script:
    template "trimmomatic.py"

}


_fastqc_trimmomatic_out_1_2.into{ fastqc_trimmomatic_out_1_2;mash_screen_in_1_3;megahit_in_1_4;_LAST_fastq_3_6 }


if (binding.hasVariable("SIDE_mashSketchOutChannel_2_4")){
    IN_reference_file_2_4 = SIDE_mashSketchOutChannel_2_4
} else {
    IN_reference_file_2_4 = Channel.value(params.refFile_2_4)
}

// check if noWinner is provided or not
winnerVar = (params.noWinner_2_4 == false) ? "-w" : ""

// process to run mashScreen and sort the output into
// sortedMashScreenResults_{sampleId}.txt
process mashScreen_2_4 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 2_4 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 2_4 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 2_4 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId mash_screen_2_4 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    input:
    set sample_id, file(reads) from mash_screen_in_1_3
    val refFile from IN_reference_file_2_4

    output:
    set sample_id, file("sortedMashScreenResults*.txt") into mashScreenResults_2_4
    set sample_id, val("2_4_mashScreen"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_mashScreen_2_4
set sample_id, val("mashScreen_2_4"), val("2_4"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_mashScreen_2_4
file ".versions"

    """
    mash screen -i ${params.identity_2_4} -v ${params.pValue_2_4} -p \
    ${task.cpus} ${winnerVar} ${refFile} ${reads} > mashScreenResults_${sample_id}.txt
    sort -gr mashScreenResults_${sample_id}.txt > sortedMashScreenResults_${sample_id}.txt
    """
}

// process to parse the output to json format
process mashOutputJson_2_4 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 2_4 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 2_4 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 2_4 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId mash_screen_2_4 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    publishDir 'results/mashscreen/mashscreen_json_2_4', mode: 'copy'

    input:
    set sample_id, file(mashtxt) from mashScreenResults_2_4

    output:
    set sample_id, file("sortedMashScreenResults*.json") optional true into mashScreenOutputChannel_2_4
    set sample_id, val("2_4_mashOutputJson"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_mashOutputJson_2_4
set sample_id, val("mashOutputJson_2_4"), val("2_4"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_mashOutputJson_2_4
file ".versions"

    script:
    template "mashscreen2json.py"
}


if ( params.megahitKmers_3_5.toString().split(" ").size() <= 1 ){
    if (params.megahitKmers_3_5.toString() != 'auto'){
        exit 1, "'megahitKmers_3_5' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.megahitKmers_3_5}"
    }
}
IN_megahit_kmers_3_5 = Channel.value(params.megahitKmers_3_5)

clear = params.clearInput_3_5 ? "true" : "false"
checkpointClear_3_5 = Channel.value(clear)

process megahit_3_5 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 3_5 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 3_5 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 3_5 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId megahit_3_5 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir 'results/assembly/megahit_3_5/', pattern: '*_megahit*.fasta', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), max_len from megahit_in_1_4.join(SIDE_max_len_3_5)
    val kmers from IN_megahit_kmers_3_5
    val clear from checkpointClear_3_5

    output:
    set sample_id, file('*megahit*.fasta') into megahit_out_3_4
    set sample_id, file('megahit/intermediate_contigs/k*.contigs.fa') into IN_fastg3_5
    set sample_id, val("3_5_megahit"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_megahit_3_5
set sample_id, val("megahit_3_5"), val("3_5"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_megahit_3_5
file ".versions"

    script:
    template "megahit.py"

}

fastg = params.fastg_3_5 ? "true" : "false"
process megahit_fastg_3_5{

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 3_5 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 3_5 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 3_5 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId megahit_3_5 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir "results/assembly/megahit_3_5/$sample_id", pattern: "*.fastg"

    input:
    set sample_id, file(kmer_files) from IN_fastg3_5
    val run_fastg from fastg

    output:
    file "*.fastg" optional true
    set sample_id, val("3_5_megahit_fastg"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_megahit_fastg_3_5
set sample_id, val("megahit_fastg_3_5"), val("3_5"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_megahit_fastg_3_5
file ".versions"

    script:
    """
    if [ ${run_fastg} == "true" ]
    then
        for kmer_file in ${kmer_files};
        do
            echo \$kmer_file
            k=\$(echo \$kmer_file | cut -d '.' -f 1);
            echo \$k
            megahit_toolkit contig2fastg \$k \$kmer_file > \$kmer_file'.fastg';
        done
    fi
    """
}


IN_min_contig_lenght_3_6 = Channel.value(params.min_contig_lenght_3_6)
IN_max_iteration_3_6 = Channel.value(params.max_iteration_3_6)
IN_prob_threshold_3_6 = Channel.value(params.prob_threshold_3_6)

clear = params.clearInput_3_6 ? "true" : "false"
checkpointClear_3_6 = Channel.value(clear)

process maxbin2_3_6 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 3_6 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 3_6 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 3_6 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId maxbin2_3_6 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    publishDir "results/maxbin2_3_6/${sample_id}/"

    input:
    set sample_id, file(assembly), file(fastq) from megahit_out_3_4.join(_LAST_fastq_3_6)
    val minContigLenght from IN_min_contig_lenght_3_6
    val maxIterations from IN_max_iteration_3_6
    val probThreshold from IN_prob_threshold_3_6
    val clear from checkpointClear_3_6

    output:
    set sample_id, file(assembly), file ('*_maxbin.*.fasta'), file ('bin_status.txt') into binCh_3_6
    file '*_maxbin.{abundance,log,summary}'
    set sample_id, file("*_maxbin.summary") into intoReport_3_6

    set sample_id, val("3_6_maxbin2"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_maxbin2_3_6
set sample_id, val("maxbin2_3_6"), val("3_6"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_maxbin2_3_6
file ".versions"

    script:
    """
    {
        run_MaxBin.pl -contig ${assembly} -out ${sample_id}_maxbin -reads ${fastq[0]} -reads2 ${fastq[1]} \
        -thread $task.cpus -min_contig_length ${minContigLenght} -max_iteration ${maxIterations} \
        -prob_threshold ${probThreshold}

        echo pass > .status

        #in case maxbin fails to bin sequences for a sample:
        if ls *_maxbin.*.fasta 1> /dev/null 2>&1; then echo "true" > bin_status.txt; else echo "false" \
        > false_maxbin.0.fasta; echo "false" > bin_status.txt; fi


        if [ "$clear" = "true" ];
        then
            work_regex=".*/work/.{2}/.{30}/.*"
            file_source1=\$(readlink -f \$(pwd)/${fastq[0]})
            file_source2=\$(readlink -f \$(pwd)/${fastq[1]})
            assembly_file=\$(readlink -f \$(pwd)/${assembly})
            if [[ "\$file_source1" =~ \$work_regex ]]; then
                rm \$file_source1 \$file_source2 \$assembly_file
            fi
        fi
    } || {
        echo fail > .status
    }
    """
}

process report_maxbin2_3_6{

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 3_6 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 3_6 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 3_6 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId maxbin2_3_6 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    input:
    set sample_id, file(tsv) from  intoReport_3_6

    output:
    set sample_id, val("3_6_report_maxbin2"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_report_maxbin2_3_6
set sample_id, val("report_maxbin2_3_6"), val("3_6"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_report_maxbin2_3_6
file ".versions"

    script:
    template "process_tsv.py"

}

// If maxbin fails to obtain bins for a sample, the workflow continues with the original assembly
_maxbin2_out_3_5 = Channel.create()

OUT_binned = Channel.create()
OUT_unbinned = Channel.create()

failedBinning = Channel.create()
successfulBinning = Channel.create()

binCh_3_6.choice(failedBinning, successfulBinning){ it -> it[3].text == "false\n" ? 0 : 1 }

failedBinning.map{ it -> [it[0], it[1]] }.into(OUT_unbinned)

successfulBinning.map{ it -> [it[2].toString().tokenize('/').last().tokenize('.')[0..-2].join('.'), it[2]]}
    .transpose()
    .map{it -> [it[1].toString().tokenize('/').last().tokenize('.')[0..-2].join('.'),it[1]]}
    .into(OUT_binned)

OUT_binned.mix(OUT_unbinned).set{ _maxbin2_out_3_5 }



_maxbin2_out_3_5.into{ maxbin2_out_3_5;mash_dist_in_3_6;mlst_in_3_7;kraken2_fasta_in_3_8;abricate_in_3_9 }

IN_shared_hashes_4_7 = Channel.value(params.shared_hashes_4_7)

IN_mash_dist_input = Channel.create()
// If the side channel with the sketch exists, join the corresponding .msh file
// with the appropriate sample_id
if (binding.hasVariable("SIDE_mashSketchOutChannel_4_7")){
    mash_dist_in_3_6
        .join(SIDE_mashSketchOutChannel_4_7)
        .into(IN_mash_dist_input)
// Otherwise, always use the .msh file provided in the docker image
} else {
    mash_dist_in_3_6
        .map{ it -> [it[0], it[1], params.refFile_4_7] }
        .into(IN_mash_dist_input)
}

// runs mash dist
process runMashDist_4_7 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 4_7 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 4_7 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 4_7 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId mash_dist_4_7 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    publishDir 'results/mashdist/mashdist_txt_4_7/'

    input:
    set sample_id, file(fasta), refFile from IN_mash_dist_input

    output:
    set sample_id, file(fasta), file("*_mashdist.txt") into mashDistOutChannel_4_7
    set sample_id, val("4_7_runMashDist"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_runMashDist_4_7
set sample_id, val("runMashDist_4_7"), val("4_7"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_runMashDist_4_7
file ".versions"

    """
    mash dist -i -p ${task.cpus} -v ${params.pValue_4_7} \
    -d ${params.mash_distance_4_7} ${refFile} ${fasta} > ${fasta}_mashdist.txt
    """

}

// parses mash dist output to a json file that can be imported into pATLAS
process mashDistOutputJson_4_7 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 4_7 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 4_7 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 4_7 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId mash_dist_4_7 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    publishDir 'results/mashdist/mashdist_json_4_7/'

    input:
    set sample_id, fasta, file(mashtxt) from mashDistOutChannel_4_7
    val shared_hashes from IN_shared_hashes_4_7

    output:
    set sample_id, file("*.json") optional true into mash_dist_out_4_6
    set sample_id, val("4_7_mashDistOutputJson"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_mashDistOutputJson_4_7
set sample_id, val("mashDistOutputJson_4_7"), val("4_7"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_mashDistOutputJson_4_7
file ".versions"

    script:
    template "mashdist2json.py"

}



process mlst_5_8 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 5_8 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 5_8 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 5_8 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId mlst_5_8 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set sample_id, file(assembly) from mlst_in_3_7

    output:
    file '*.mlst.txt' into LOG_mlst_5_8
    set sample_id, file(assembly), file(".status") into MAIN_mlst_out_5_8
    set sample_id, val("5_8_mlst"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_mlst_5_8
set sample_id, val("mlst_5_8"), val("5_8"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_mlst_5_8
file ".versions"

    script:
    """
    {
        expectedSpecies=${params.mlstSpecies_5_8}
        mlst $assembly >> ${sample_id}.mlst.txt
        mlstSpecies=\$(cat *.mlst.txt | cut -f2)
        json_str="{'expectedSpecies':\'\$expectedSpecies\',\
            'species':'\$mlstSpecies',\
            'st':'\$(cat *.mlst.txt | cut -f3)',\
            'tableRow':[{'sample':'${sample_id}','data':[\
                {'header':'MLST species','value':'\$mlstSpecies','table':'typing'},\
                {'header':'MLST ST','value':'\$(cat *.mlst.txt | cut -f3)','table':'typing'}]}]}"
        echo \$json_str > .report.json

        if [ ! \$mlstSpecies = \$expectedSpecies ];
        then
            printf fail > .status
        else
            printf pass > .status
        fi

    } || {
        printf fail > .status
    }
    """
}

process compile_mlst_5_8 {

    publishDir "results/annotation/mlst_5_8/"

    input:
    file res from LOG_mlst_5_8.collect()

    output:
    file "mlst_report.tsv"

    script:
    """
    cat $res >> mlst_report.tsv
    """
}

mlst_out_5_7 = Channel.create()
MAIN_mlst_out_5_8
    .filter{ it[2].text != "fail" }
    .map{ [it[0], it[1]] }
    .set{ mlst_out_5_7 }




IN_kraken2_DB_6_9 = Channel.value(params.kraken2DB_6_9)


//Process to run Kraken2
process kraken2_fasta_6_9 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 6_9 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 6_9 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 6_9 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId kraken2_fasta_6_9 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    publishDir "results/taxonomy/kraken2/", pattern: "*.txt"

    input:
    set sample_id, file(assembly) from kraken2_fasta_in_3_8
    val krakenDB from IN_kraken2_DB_6_9

    output:
    file("${sample_id}_kraken_report.txt")
    set sample_id, val("6_9_kraken2_fasta"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_kraken2_fasta_6_9
set sample_id, val("kraken2_fasta_6_9"), val("6_9"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_kraken2_fasta_6_9
file ".versions"

    script:
    """
    kraken2 --memory-mapping --threads $task.cpus --report ${sample_id}_kraken_report.txt --db ${krakenDB} ${assembly}
    """
}


if ( params.abricateDataDir_7_10 ){
    if ( !file(params.abricateDataDir_7_10).exists() ){
        exit 1, "'abricateDataDir_7_10' data directory was not found: '${params.abricateDatabases_7_10}'"
    }
    dataDirOpt = "--datadir ${params.abricateDataDir_7_10}"
} else {
    dataDirOpt = ""
}

if ( !params.abricateMinId_7_10.toString().isNumber() ){
    exit 1, "'abricateMinId_7_10' parameter must be a number. Provide value: '${params.abricateMinId_7_10}'"
}

if ( !params.abricateMinCov_7_10.toString().isNumber() ){
    exit 1, "'abricateMinCov_7_10' parameter must be a number. Provide value: '${params.abricateMinCov_7_10}'"
}


process abricate_7_10 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 7_10 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 7_10 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 7_10 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId abricate_7_10 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { "${sample_id} ${db}" }
    publishDir "results/annotation/abricate_7_10/${sample_id}"

    input:
    set sample_id, file(assembly) from abricate_in_3_9
    each db from params.abricateDatabases_7_10
    val min_id from Channel.value(params.abricateMinId_7_10)
    val min_cov from Channel.value(params.abricateMinCov_7_10)

    output:
    file '*.tsv' into abricate_out_7_10
    set sample_id, val("7_10_abricate_$db"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_abricate_7_10
set sample_id, val("abricate_7_10_$db"), val("7_10"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_abricate_7_10
file ".versions"

    script:
    """
    {
        # Run abricate
        abricate $dataDirOpt --minid $min_id --mincov $min_cov --db $db $assembly > ${sample_id}_abr_${db}.tsv
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}


process process_abricate_7_10 {

    tag "process_abricate_7_10"

    // Send POST request to platform
    
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh"
        afterScript "report_POST.sh $params.projectId $params.pipelineId 7_10 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId abricate_7_10 \"$params.platformSpecies\" false"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh"
        }
    

    input:
    file abricate_file from abricate_out_7_10.collect()

    output:
    set val('process_abricate'), val("7_10_process_abricate"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_process_abricate_7_10
set val('process_abricate'), val("process_abricate_7_10"), val("7_10"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_process_abricate_7_10
file ".versions"

    script:
    template "process_abricate.py"


}




/** STATUS
Reports the status of a sample in any given process.
*/
process status {

    tag { sample_id }
    publishDir "pipeline_status/$task_name"

    input:
    set sample_id, task_name, status, warning, fail, file(log) from STATUS_integrity_coverage_1_1.mix(STATUS_remove_host_1_2,STATUS_report_remove_host_1_2,STATUS_fastqc_1_3,STATUS_fastqc_report_1_3,STATUS_trimmomatic_1_3,STATUS_mashScreen_2_4,STATUS_mashOutputJson_2_4,STATUS_megahit_3_5,STATUS_megahit_fastg_3_5,STATUS_maxbin2_3_6,STATUS_report_maxbin2_3_6,STATUS_runMashDist_4_7,STATUS_mashDistOutputJson_4_7,STATUS_mlst_5_8,STATUS_kraken2_fasta_6_9,STATUS_abricate_7_10,STATUS_process_abricate_7_10)

    output:
    file '*.status' into master_status
    file '*.warning' into master_warning
    file '*.fail' into master_fail
    file '*.log'

    """
    echo $sample_id, $task_name, \$(cat $status) > ${sample_id}_${task_name}.status
    echo $sample_id, $task_name, \$(cat $warning) > ${sample_id}_${task_name}.warning
    echo $sample_id, $task_name, \$(cat $fail) > ${sample_id}_${task_name}.fail
    echo "\$(cat .command.log)" > ${sample_id}_${task_name}.log
    """
}

process compile_status_buffer {

    input:
    file status from master_status.buffer( size: 5000, remainder: true)
    file warning from master_warning.buffer( size: 5000, remainder: true)
    file fail from master_fail.buffer( size: 5000, remainder: true)

    output:
    file 'master_status_*.csv' into compile_status_buffer
    file 'master_warning_*.csv' into compile_warning_buffer
    file 'master_fail_*.csv' into compile_fail_buffer

    """
    cat $status >> master_status_${task.index}.csv
    cat $warning >> master_warning_${task.index}.csv
    cat $fail >> master_fail_${task.index}.csv
    """
}

process compile_status {

    publishDir 'reports/status'

    input:
    file status from compile_status_buffer.collect()
    file warning from compile_warning_buffer.collect()
    file fail from compile_fail_buffer.collect()

    output:
    file "*.csv"

    """
    cat $status >> master_status.csv
    cat $warning >> master_warning.csv
    cat $fail >> master_fail.csv
    """

}


/** Reports
Compiles the reports from every process
*/
process report {

    tag { sample_id }

    input:
    set sample_id,
            task_name,
            pid,
            report_json,
            version_json,
            trace from REPORT_integrity_coverage_1_1.mix(REPORT_remove_host_1_2,REPORT_report_remove_host_1_2,REPORT_fastqc_1_3,REPORT_fastqc_report_1_3,REPORT_trimmomatic_1_3,REPORT_mashScreen_2_4,REPORT_mashOutputJson_2_4,REPORT_megahit_3_5,REPORT_megahit_fastg_3_5,REPORT_maxbin2_3_6,REPORT_report_maxbin2_3_6,REPORT_runMashDist_4_7,REPORT_mashDistOutputJson_4_7,REPORT_mlst_5_8,REPORT_kraken2_fasta_6_9,REPORT_abricate_7_10,REPORT_process_abricate_7_10)

    output:
    file "*" optional true into master_report

    """
    prepare_reports.py $report_json $version_json $trace $sample_id $task_name 1 $pid $workflow.scriptId $workflow.runName
    """

}


process compile_reports {

    publishDir "pipeline_report/", mode: "copy"

    if ( params.reportHTTP != null ){
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH;"
        afterScript "metadata_POST.sh $params.projectId $params.pipelineId 0 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId 0 \"$params.platformSpecies\""
    }

    input:
    file report from master_report.collect()
    file forks from Channel.fromPath("${workflow.projectDir}/.forkTree.json")
    file dag from Channel.fromPath("${workflow.projectDir}/.treeDag.json")
    file js from Channel.fromPath("${workflow.projectDir}/resources/main.js.zip")

    output:
    file "pipeline_report.json"
    file "pipeline_report.html"
    file "src/main.js"

    script:
    template "compile_reports.py"
}


/**
* A process that creates a consensus from all the outputted json files
*/
process fullConsensus {

    tag { sample_id }

    publishDir 'results/consensus_/'

    input:
    set sample_id, file(infile_list) from mashScreenOutputChannel_2_4

    output:
    file "consensus_*.json"

    script:
    template "pATLAS_consensus_json.py"

}
workflow.onComplete {
  // Display complete message
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
}

workflow.onError {
  // Display error message
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
