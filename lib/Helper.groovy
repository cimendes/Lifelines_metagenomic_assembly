class Help {

    static def start_info(Map info, String time, String profile) {

        println ""
        println "============================================================"
        println "                F L O W C R A F T"
        println "============================================================"
        println "Built using flowcraft v1.4.0"
        println ""
        if (info.containsKey("fastq")){
        int nsamples = info.fastq / 2
        println " Input FastQ                 : $info.fastq"
        println " Input samples               : $nsamples"
        }
        if (info.containsKey("fasta")){
        println " Input Fasta                 : $info.fasta"
        }
        if (info.containsKey("accessions")){
        println " Input accessions            : $info.accessions"
        }
        println " Reports are found in        : ./reports"
        println " Results are found in        : ./results"
        println " Profile                     : $profile"
        println ""
        println "Starting pipeline at $time"
        println ""

    }

    static void complete_info(nextflow.script.WorkflowMetadata wf) {

        println ""
        println "Pipeline execution summary"
        println "=========================="
        println "Completed at                 : $wf.complete"
        println "Duration                     : $wf.duration"
        println "Success                      : $wf.success"
        println "Work directory               : $wf.workDir"
        println "Exit status                  : $wf.exitStatus"
        println ""

    }

    static def print_help(Map params) {

        println ""
        println "============================================================"
        println "                F L O W C R A F T"
        println "============================================================"
        println "Built using flowcraft v1.4.0"
        println ""
        println ""
        println "Usage: "
        println "    nextflow run lifelines.nf"
        println ""
        println "       --fastq                     Path expression to paired-end fastq files. (default: $params.fastq) (default: 'fastq/*_{1,2}.*')"
        println "       "
        println "       Component 'INTEGRITY_COVERAGE_1_1'"
        println "       ----------------------------------"
        println "       --genomeSize_1_1            Genome size estimate for the samples in Mb. It is used to estimate the coverage and other assembly parameters andchecks (default: 1)"
        println "       --minCoverage_1_1           Minimum coverage for a sample to proceed. By default it's setto 0 to allow any coverage (default: 0)"
        println "       "
        println "       Component 'REMOVE_HOST_1_2'"
        println "       ---------------------------"
        println "       --refIndex_1_2              Specifies the reference indexes to be provided to bowtie2. (default: '/index_hg19/hg19')"
        println "       --clearInput_1_2            Permanently removes temporary input files. This option is only useful to remove temporary files in large workflows and prevents nextflow's resume functionality. Use with caution. (default: false)"
        println "       "
        println "       Component 'FASTQC_TRIMMOMATIC_1_3'"
        println "       ----------------------------------"
        println "       --adapters_1_3              Path to adapters files, if any. (default: 'None')"
        println "       --trimSlidingWindow_1_3     Perform sliding window trimming, cutting once the average quality within the window falls below a threshold. (default: '5:20')"
        println "       --trimLeading_1_3           Cut bases off the start of a read, if below a threshold quality. (default: 3)"
        println "       --trimTrailing_1_3          Cut bases of the end of a read, if below a threshold quality. (default: 3)"
        println "       --trimMinLength_1_3         Drop the read if it is below a specified length. (default: 55)"
        println "       --clearInput_1_3            Permanently removes temporary input files. This option is only useful to remove temporary files in large workflows and prevents nextflow's resume functionality. Use with caution. (default: false)"
        println "       "
        println "       Component 'MASH_SCREEN_2_4'"
        println "       ---------------------------"
        println "       --noWinner_2_4              A variable that enables the use of -w option for mash screen. (default: false)"
        println "       --pValue_2_4                P-value cutoff for the distance estimation between two sequences to be included in the output. (default: 0.05)"
        println "       --identity_2_4              The percentage of identity between the reads input and the reference sequence (default: 0.9)"
        println "       --refFile_2_4               Specifies the reference file to be provided to mash. It can either be a fasta or a .msh reference sketch generated by mash. (default: '/ngstools/data/plasmid_db_reference.msh')"
        println "       "
        println "       Component 'MEGAHIT_3_5'"
        println "       -----------------------"
        println "       --megahitKmers_3_5          If 'auto' the megahit k-mer lengths will be determined from the maximum read length of each assembly. If 'default', megahit will use the default k-mer lengths. (default: $params.megahitKmers) (default: 'auto')"
        println "       --fastg_3_5                 Converts megahit intermediate contigs to fastg (default: false)"
        println "       --clearInput_3_5            Permanently removes temporary input files. This option is only useful to remove temporary files in large workflows and prevents nextflow's resume functionality. Use with caution. (default: false)"
        println "       "
        println "       Component 'ASSEMBLY_MAPPING_3_6'"
        println "       --------------------------------"
        println "       --minAssemblyCoverage_3_6   In auto, the default minimum coverage for each assembled contig is 1/3 of the assembly mean coverage or 10x, if the mean coverage is below 10x (default: 'auto')"
        println "       --AMaxContigs_3_6           A warning is issued if the number of contigs is overthis threshold. (default: 100)"
        println "       --genomeSize_3_6            Genome size estimate for the samples. It is used to check the ratio of contig number per genome MB (default: 2.1)"
        println "       "
        println "       Component 'PILON_3_7'"
        println "       ---------------------"
        println "       --clearInput_3_7            Permanently removes temporary input files. This option is only useful to remove temporary files in large workflows and prevents nextflow's resume functionality. Use with caution. (default: false)"
        println "       "
        println "       Component 'MASH_DIST_4_8'"
        println "       -------------------------"
        println "       --pValue_4_8                P-value cutoff for the distance estimation between two sequences to be included in the output. (default: 0.05)"
        println "       --mash_distance_4_8         Sets the maximum distance between two sequences to be included in the output. (default: 0.1)"
        println "       --shared_hashes_4_8         Sets a minimum percentage of hashes shared between two sequences in order to include its result in the output. (default: 0.8)"
        println "       --refFile_4_8               Specifies the reference file to be provided to mash. It can either be a fasta or a .msh reference sketch generated by mash. (default: '/ngstools/data/plasmid_db_reference.msh')"
        println "       "
        println "       Component 'MLST_5_9'"
        println "       --------------------"
        println "       --mlstSpecies_5_9           Specify the expected species for MLST checking. (default: null)"
        println "       "
        println "       Component 'KRAKEN2_FASTA_6_10'"
        println "       ------------------------------"
        println "       --kraken2DB_6_10            Specifies kraken2 database. Requires full path if database not on KRAKEN2_DB_PATH. (default: 'minikraken2_v1_8GB')"
        println "       "
        println "       Component 'ABRICATE_7_11'"
        println "       -------------------------"
        println "       --abricateDatabases_7_11    Specify the databases for abricate. (default: ['resfinder', 'card', 'vfdb', 'plasmidfinder', 'virulencefinder', 'bacmet'])"
        println "       --abricateDataDir_7_11      Specify the full path location of the database folders. (default: null)"
        println "       --abricateMinId_7_11        Minimum DNA %identity. (default: 75)"
        println "       --abricateMinCov_7_11       Minimum DNA %coverage. (default: 0)"
        
    }

}

class CollectInitialMetadata {

    public static void print_metadata(nextflow.script.WorkflowMetadata workflow){

        def treeDag = new File("${workflow.projectDir}/.treeDag.json").text
        def forkTree = new File("${workflow.projectDir}/.forkTree.json").text

        def metadataJson = "{'nfMetadata':{'scriptId':'${workflow.scriptId}',\
'scriptName':'${workflow.scriptName}',\
'profile':'${workflow.profile}',\
'container':'${workflow.container}',\
'containerEngine':'${workflow.containerEngine}',\
'commandLine':'${workflow.commandLine}',\
'runName':'${workflow.runName}',\
'sessionId':'${workflow.sessionId}',\
'projectDir':'${workflow.projectDir}',\
'launchDir':'${workflow.launchDir}',\
'startTime':'${workflow.start}',\
'dag':${treeDag},\
'forks':${forkTree}}}"

        def json = metadataJson.replaceAll("'", '"')

        def jsonFile = new File(".metadata.json")
        jsonFile.write json
    }
}