## ------------------------------------------------------------------------------------ ##
## To run using environment modules: 
## we need exactly this snakemake version in order to properly track i/o
## module load snakemake/7.8.3-foss-2020b
## snakemake --cores 40 --use-envmodules -npr
## (remove -npr to actually run the workflow)
## ------------------------------------------------------------------------------------ ##

## For a new dataset, need to provide the following manually:
## - a shell script that downloads the data into a folder named FASTQ
## - a text file 'samples.txt' with sample annotations and associations with FASTQ files
## - a general config file with info about the data set (read composition etc)
## - additional config files for each method

## If mutscan is updated, edit update_mutscan.txt
## to trigger a rerun of all rules that depend on mutscan

## ------------------------------------------------------------------------------------ ##
## Preamble
## ------------------------------------------------------------------------------------ ##
## Define modules
Rmod = "R-BioC/4.2-3.15-foss-2020b"
Rmoddevel = "R-BioC/devel-foss-2020b"  ## Only used for generating the final plots for the paper
popplermod = "poppler/21.06.1-GCC-10.2.0"  ## For pdftools
Pandocmod = "Pandoc/2.10"
dimsummod = "dimsum/1.2.11"
enrich2mod = "enrich2/1.3.1"

Rbin = "R CMD BATCH --no-restore --no-save"

## Define the datasets to include
datasets = ["Bolognesi_TDP43_290_331", "Diss_FOS", "Diss_FOS_JUN", "Li_tRNA_sel30"]
## For time reasons, run Enrich2 only for a subset of the datasets
datasets_enrich2 = ["Diss_FOS"]

## For each of the datasets, read the dataset/samples.txt file and the set 
## of configuration files. If a new method (with a config file) is added,
## read that here as well.  
import pandas as pd
import yaml
datasetsVec = []
samplesVec = []
configDict = {}
for ds in datasets:
    samples = pd.read_csv(ds + "/samples.txt", sep = '\t')
    samples_unique = list(set(samples.mutscan_name.values.tolist()))
    samplesVec.extend(samples_unique)
    datasetsVec.extend([ds] * len(samples_unique))
    
    configDict[ds] = {}
    ## Read dataset config
    with open(ds + "/dataset_config.yml", "r") as stream:
        try:
            configDict[ds]["dataset"] = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    ## Read mutscan config
    with open(ds + "/mutscan_config.yml", "r") as stream:
        try:
            configDict[ds]["mutscan"] = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    ## Read dimsum config
    with open(ds + "/dimsum_config.yml", "r") as stream:
        try:
            configDict[ds]["dimsum"] = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    ## Read Enrich2 config
    if ds in datasets_enrich2:
        with open(ds + "/enrich2_config.json", "r") as stream:
            try:
                configDict[ds]["enrich2"] = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)

# print(configDict)
# print(datasetsVec)
# print(samplesVec)

## ------------------------------------------------------------------------------------ ##
## Default rules
## ------------------------------------------------------------------------------------ ##
rule all:
    input:
        expand("{dataset}/get_data_checkpoint.txt", dataset = datasets),
        expand("{dataset}/output/mutscan/digestFastqs/{sample}.rds", zip, dataset = datasetsVec, sample = samplesVec),
        expand("{dataset}/output/mutscan/{dataset}_mutscan_se.rds", dataset = datasets),
        expand("{dataset}/output/mutscan/{dataset}_mutscan_testresults.rds", dataset = datasets),
        expand("{dataset}/output/mutscan/{dataset}_mutscan_nullcomparison.rds", dataset = datasets),
        expand("{dataset}/output/dimsum/report.html", dataset = datasets),
        expand("{dataset}/output/enrich2/out.log", dataset = datasets_enrich2),
        "benchmark/Rdata/benchmark_dataframe.rds",
        "benchmark/pdf/benchmark_results.pdf",
        "benchmark/pdf/nullcomparison.pdf",
        expand("{dataset}/output/summary_data/summary_table.rds", dataset = datasets),
        expand("{dataset}/plots/compare_output_upset.png", dataset = datasets),
        expand(["{dataset}/plots/compare_counts_{sample}_pairs_allvariants.png"], zip, dataset = datasetsVec, sample = samplesVec),
        ## scaling (ncores)
        expand("{dataset}/output/mutscan_scaling/digestFastqs_{sample}_{ncores}cores.rds", dataset = ["Li_tRNA_sel30"], sample = ["inputA"], ncores = [1, 2, 5, 10, 15, 20]),
        "benchmark/Rdata/benchmark_dataframe_scaling.rds",
        "benchmark/pdf/benchmark_results_scaling.pdf",
        ## scaling (nreads)
        expand("{dataset}/output/mutscan_scaling/digestFastqs_{sample}_{nreads}reads.rds", dataset = ["Li_tRNA_sel30"], sample = ["inputA"], nreads = [100000,500000,1000000,5000000,10000000,50000000,100000000,250000000,391420100]),
        "benchmark/Rdata/benchmark_dataframe_scaling_nreads.rds",
        "benchmark/pdf/benchmark_results_scaling_nreads.pdf",
        ## case study
        expand("{dataset}/CaseStudy/{dataset}_CaseStudy.html", dataset = ["Diss_FOS_JUN"])

rule paper:
    input:
        "paper/supplement/mutscan_supplement.pdf",
        "paper/figures/Fig6.png",
        "paper/R_package_versions.txt"

## ------------------------------------------------------------------------------------ ##
## Get data
## ------------------------------------------------------------------------------------ ##
rule get_data:
    input:
        script = "{dataset}/get_data.sh"
    output:
        checkpoint = "{dataset}/get_data_checkpoint.txt",
    shell:
        '''
        cd {wildcards.dataset} && \
        bash $(basename {input.script}) && \
        bash ../scripts/check_fastq_files.sh samples.txt && \
        touch $(basename {output.checkpoint})
        '''

## ------------------------------------------------------------------------------------ ##
## Run mutscan
## ------------------------------------------------------------------------------------ ##
## Quantification (separately for each sample)
rule mutscan_dms_quant:
    input: 
        "update_mutscan.txt",
        script = "scripts/mutscan_dms_quant.R",
        dataset_config = "{dataset}/dataset_config.yml",
        mutscan_config = "{dataset}/mutscan_config.yml",
        checkpoint = "{dataset}/get_data_checkpoint.txt",
        samplefile = "{dataset}/samples.txt"
    envmodules:
        f"{Rmod}"
    benchmark: 
        "{dataset}/output/benchmark/mutscan_digestFastqs_{sample}.txt"
    threads: 
        10
    log: 
        rout = "{dataset}/Rout/mutscan_digestFastqs_{sample}.Rout"
    output:
        rds = "{dataset}/output/mutscan/digestFastqs/{sample}.rds"
    shell:
        '''
        {Rbin} "--args dataset='{wildcards.dataset}' sample='{wildcards.sample}' nthreads={threads} nreads=-1 outrds='{output.rds}'" {input.script} {log.rout}
        '''

## Summarization
def get_sample_rds(wildcards):
    samplesVecFilt = [samplesVec[s] for s in range(len(samplesVec)) if datasetsVec[s] == wildcards.dataset]
    return [(wildcards.dataset + "/output/mutscan/digestFastqs/" + s + ".rds") for s in samplesVecFilt] 

rule mutscan_dms_summarize:
    input:
        "update_mutscan.txt",
        samplerds = get_sample_rds, 
        script = "scripts/mutscan_dms_summarize.R",
        annots = "{dataset}/samples.txt",
        mutscan_config = "{dataset}/mutscan_config.yml"
    envmodules:
        f"{Rmod}",
        f"{Pandocmod}"
    benchmark: 
        "{dataset}/output/benchmark/mutscan_summarize.txt"
    log: 
        rout = "{dataset}/Rout/mutscan_summarize.Rout"
    output:
        rds = "{dataset}/output/mutscan/{dataset}_mutscan_se.rds",
        rdsfilt = "{dataset}/output/mutscan/{dataset}_mutscan_se_filt.rds",
    shell:
        '''
        {Rbin} "--args dataset='{wildcards.dataset}' outrds='{output.rds}'" {input.script} {log.rout}
        '''

rule mutscan_dms_test:
    input:
        "update_mutscan.txt",
        inrds = "{dataset}/output/mutscan/{dataset}_mutscan_se_filt.rds", 
        script = "scripts/mutscan_dms_test.R",
        mutscan_config = "{dataset}/mutscan_config.yml"
    envmodules:
        f"{Rmod}",
    benchmark: 
        "{dataset}/output/benchmark/mutscan_test.txt"
    log: 
        rout = "{dataset}/Rout/mutscan_test.Rout"
    output:
        rds = "{dataset}/output/mutscan/{dataset}_mutscan_testresults.rds",
    shell:
        '''
        {Rbin} "--args dataset='{wildcards.dataset}' inrds='{input.inrds}' outrds='{output.rds}'" {input.script} {log.rout}
        '''

rule mutscan_dms_nullcomparison:
    input:
        "update_mutscan.txt",
        inrds = "{dataset}/output/mutscan/{dataset}_mutscan_se_filt.rds", 
        script = "scripts/mutscan_dms_nullcomparison.R",
        mutscan_config = "{dataset}/mutscan_config.yml"
    envmodules:
        f"{Rmod}",
    benchmark: 
        "{dataset}/output/benchmark/mutscan_nullcomparison.txt"
    log: 
        rout = "{dataset}/Rout/mutscan_nullcomparison.Rout"
    output:
        rds = "{dataset}/output/mutscan/{dataset}_mutscan_nullcomparison.rds",
    shell:
        '''
        {Rbin} "--args dataset='{wildcards.dataset}' inrds='{input.inrds}' outrds='{output.rds}'" {input.script} {log.rout}
        '''

## ------------------------------------------------------------------------------------ ##
## Run mutscan to check scaling with the number of cores/reads
## ------------------------------------------------------------------------------------ ##
rule mutscan_dms_quant_scaling:
    input: 
        "update_mutscan.txt",
        script = "scripts/mutscan_dms_quant.R",
        dataset_config = "{dataset}/dataset_config.yml",
        mutscan_config = "{dataset}/mutscan_config.yml",
        checkpoint = "{dataset}/get_data_checkpoint.txt",
        samplefile = "{dataset}/samples.txt"
    envmodules:
        f"{Rmod}"
    benchmark: 
        repeat("{dataset}/output/benchmark/mutscan_scaling_digestFastqs_{sample}_{ncores}cores.txt", 5)
    threads: 
        20  ## set to 20 as this is the highest number that will actually be used
    log: 
        rout = "{dataset}/Rout/mutscan_scaling_digestFastqs_{sample}_{ncores}cores.Rout"
    output:
        rds = "{dataset}/output/mutscan_scaling/digestFastqs_{sample}_{ncores}cores.rds"
    shell:
        '''
        {Rbin} "--args dataset='{wildcards.dataset}' sample='{wildcards.sample}' nthreads={wildcards.ncores} nreads=-1 outrds='{output.rds}'" {input.script} {log.rout}
        '''

rule mutscan_dms_quant_scaling_nreads:
    input: 
        "update_mutscan.txt",
        script = "scripts/mutscan_dms_quant.R",
        dataset_config = "{dataset}/dataset_config.yml",
        mutscan_config = "{dataset}/mutscan_config.yml",
        checkpoint = "{dataset}/get_data_checkpoint.txt",
        samplefile = "{dataset}/samples.txt"
    envmodules:
        f"{Rmod}"
    benchmark: 
        repeat("{dataset}/output/benchmark/mutscan_scaling_digestFastqs_{sample}_{nreads}reads.txt", 5)
    threads: 
        10
    log: 
        rout = "{dataset}/Rout/mutscan_scaling_digestFastqs_{sample}_{nreads}reads.Rout"
    output:
        rds = "{dataset}/output/mutscan_scaling/digestFastqs_{sample}_{nreads}reads.rds"
    shell:
        '''
        {Rbin} "--args dataset='{wildcards.dataset}' sample='{wildcards.sample}' nthreads={threads} nreads={wildcards.nreads} outrds='{output.rds}'" {input.script} {log.rout}
        '''

## ------------------------------------------------------------------------------------ ##
## Run DiMSum
## ------------------------------------------------------------------------------------ ##
rule dimsumwrap:
    input:
        checkpoint = "{dataset}/get_data_checkpoint.txt",
        samplefile = "{dataset}/samples.txt",
        dimsum_config = "{dataset}/dimsum_config.yml"
    output:
        report = "{dataset}/output/dimsum/report_wrap.html"
    benchmark: 
        "{dataset}/output/benchmark/dimsum_wrap.txt"
    envmodules:
        f"{dimsummod}"
    threads: 
        10
    params:
        fastqdir = lambda wildcards: wildcards.dataset + "/FASTQ",
        fastqFileExtension = lambda wildcards: configDict[wildcards.dataset]["dataset"]["fastqFileExtension"],
        stranded = lambda wildcards: configDict[wildcards.dataset]["dataset"]["stranded"],
        paired = lambda wildcards: configDict[wildcards.dataset]["dataset"]["paired"],
        gzipped = lambda wildcards: configDict[wildcards.dataset]["dataset"]["gzipped"],
        cutadaptCommand = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["cutadaptCommand"],
        cutadaptMinLength = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["cutadaptMinLength"],
        cutadaptErrorRate = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["cutadaptErrorRate"],
        vsearchMinQual = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["vsearchMinQual"],
        vsearchMaxee = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["vsearchMaxee"],
        vsearchMinovlen = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["vsearchMinovlen"],
        transLibrary = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["transLibrary"],
        transLibraryReverseComplement = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["transLibraryReverseComplement"],
        wildTypeSequence = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["wildTypeSequence"],
        sequenceType = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["sequenceType"],
        extraCommandsWrap = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["extraCommandsWrap"]
    shell:
        '''
        DiMSum --startStage 0 --stopStage 3 --numCores {threads} \
        --fastqFileDir {params.fastqdir} --fastqFileExtension {params.fastqFileExtension} \
        --gzipped {params.gzipped} --stranded {params.stranded} --paired {params.paired} \
        --experimentDesignPath {input.samplefile} --wildtypeSequence {params.wildTypeSequence} \
        {params.cutadaptCommand} --cutadaptMinLength {params.cutadaptMinLength} --cutadaptErrorRate {params.cutadaptErrorRate} \
        --vsearchMinQual {params.vsearchMinQual} --vsearchMaxee {params.vsearchMaxee} --vsearchMinovlen {params.vsearchMinovlen} \
        --outputPath {wildcards.dataset}/output --sequenceType {params.sequenceType} --projectName dimsum \
        --retainIntermediateFiles FALSE {params.extraCommandsWrap} \
        --transLibrary {params.transLibrary} --transLibraryReverseComplement {params.transLibraryReverseComplement} && 
        mv {wildcards.dataset}/output/dimsum/report.html {wildcards.dataset}/output/dimsum/report_wrap.html
        '''

rule dimsumsteam:
    input:
        wrapreport = "{dataset}/output/dimsum/report_wrap.html",
        samplefile = "{dataset}/samples.txt",
        dimsum_config = "{dataset}/dimsum_config.yml"
    output:
        report = "{dataset}/output/dimsum/report.html"
    benchmark:
        "{dataset}/output/benchmark/dimsum_steam.txt"
    envmodules:
        f"{dimsummod}"
    threads:
        10
    params:
        fastqdir = lambda wildcards: wildcards.dataset + "/FASTQ",
        maxSubstitutions = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["maxSubstitutions"],
        retainedReplicates = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["retainedReplicates"],
        fitnessMinInputCountAny = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["fitnessMinInputCountAny"],
        fitnessMinInputCountAll = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["fitnessMinInputCountAll"],
        fitnessMinOutputCountAny = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["fitnessMinOutputCountAny"],
        fitnessMinOutputCountAll = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["fitnessMinOutputCountAll"],
        wildTypeSequence = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["wildTypeSequence"],
        sequenceType = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["sequenceType"],
        extraCommandsSteam = lambda wildcards: configDict[wildcards.dataset]["dimsum"]["extraCommandsSteam"]
    shell:
        '''
        DiMSum --startStage 4 --stopStage 5 --numCores {threads} \
        --fastqFileDir {params.fastqdir} --outputPath {wildcards.dataset}/output --projectName dimsum \
        --experimentDesignPath {input.samplefile} --wildtypeSequence {params.wildTypeSequence} \
        --maxSubstitutions {params.maxSubstitutions} \
        --fitnessMinInputCountAny {params.fitnessMinInputCountAny} \
        --fitnessMinInputCountAll {params.fitnessMinInputCountAll} \
        --fitnessMinOutputCountAny {params.fitnessMinOutputCountAny} \
        --fitnessMinOutputCountAll {params.fitnessMinOutputCountAll} \
        --sequenceType {params.sequenceType} \
        --retainedReplicates {params.retainedReplicates} {params.extraCommandsSteam}
        '''

## ------------------------------------------------------------------------------------ ##
## Run Enrich2
## ------------------------------------------------------------------------------------ ##
## Before rerunning Enrich2, make sure that the whole output directory (especially the h5 files) is deleted
rule enrich2:
    input:
        checkpoint = "{dataset}/get_data_checkpoint.txt",
        enrich2_config = "{dataset}/enrich2_config.json"
    output:
        logfile = "{dataset}/output/enrich2/out.log",
        varcounts = "{dataset}/output/enrich2/tsv/{dataset}_exp/main_variants_counts.tsv"
    benchmark: 
        "{dataset}/output/benchmark/enrich2.txt"
    envmodules:
        f"{enrich2mod}"
    threads: 
        1
    params:
        workdir = "{dataset}"
    shell:
        '''
        enrich_cmd --no-plots --log {output.logfile} {input.enrich2_config} ratios wt
        '''

## ------------------------------------------------------------------------------------ ##
## Compare counts/scores across methods
## ------------------------------------------------------------------------------------ ##
def get_enrich2_count_file(wildcards):
    if wildcards.dataset in datasets_enrich2:
        return wildcards.dataset + "/output/enrich2/tsv/" + wildcards.dataset + "_exp/main_variants_counts.tsv"
    else:
        return []

def get_enrich2_score_file(wildcards):
    if wildcards.dataset in datasets_enrich2:
        return wildcards.dataset + "/output/enrich2/tsv/" + wildcards.dataset + "_exp/main_variants_scores_shared_full.tsv"
    else:
        return []

rule create_summary_table:
    input:
        "{dataset}/output/mutscan/{dataset}_mutscan_se_filt.rds",
        "{dataset}/output/mutscan/{dataset}_mutscan_testresults.rds",
        get_enrich2_count_file,
        get_enrich2_score_file, 
        "{dataset}/output/dimsum/report.html",
        script = "scripts/create_summary_table.R",
        conv1 = "scripts/enrich2_to_mutscan_ids.R",
        conv2 = "scripts/seqToName.cpp"
    output:
        rds = "{dataset}/output/summary_data/summary_table.rds"
    log:
        rout = "{dataset}/Rout/create_summary_table.Rout"
    envmodules:
        f"{Rmod}"
    shell:
        '''
        {Rbin} "--args dataset='{wildcards.dataset}' outrds='{output.rds}'" {input.script} {log.rout}
        '''

rule plot_comparison:
    input:
        rds = "{dataset}/output/summary_data/summary_table.rds",
        script = "scripts/plot_comparison_counts_scores.R"
    output:
        upset = "{dataset}/plots/compare_output_upset.png"
    log:
        rout = "{dataset}/Rout/plot_comparison_counts_scores.Rout"
    params:
        outbase = lambda wildcards, output: output["upset"].replace("_upset.png", "")
    envmodules:
        f"{Rmod}"
    shell:
        '''
        {Rbin} "--args sumtbl='{input.rds}' outbase='{params.outbase}'" {input.script} {log.rout}
        '''

rule plot_counts_singlesample:
    input:
        rds = "{dataset}/output/summary_data/summary_table.rds",
        script = "scripts/plot_counts_single_sample.R"
    output:
        plot = "{dataset}/plots/compare_counts_{sample}_pairs_allvariants.png"
    log:
        rout = "{dataset}/Rout/plot_counts_singlesample_{sample}.Rout"
    params:
        outbase = lambda wildcards, output: output["plot"].replace("_" + wildcards.sample + "_pairs_allvariants.png", "")
    envmodules:
        f"{Rmod}"
    shell:
        '''
        {Rbin} "--args sumtbl='{input.rds}' outbase='{params.outbase}' sname='{wildcards.sample}'" {input.script} {log.rout}
        '''

## ------------------------------------------------------------------------------------ ##
## Benchmark summary
## ------------------------------------------------------------------------------------ ##
## parse benchmark files and create data.frame rds
##   - input files are identified by {params.fileglob}
##      and assumed to have identical columns
##   - files are summarized into a data.frame with columns:
##        tool, step, sample (obtained by splitting the filename at '_'),
##        and all columns contained within the file (s, h.m.s, max_rss, ...)
rule create_benchmark_dataframe:
    input: 
        expand("{dataset}/output/mutscan/{dataset}_mutscan_testresults.rds", dataset = datasets),
        expand("{dataset}/output/dimsum/report.html", dataset = datasets),
        expand("{dataset}/output/enrich2/out.log", dataset = datasets_enrich2),
        script = "scripts/create_benchmark_dataframe.R",
    envmodules:
        f"{Rmod}"
    params:
        fileglob = "*/output/benchmark/*.txt"
    log: 
        rout = "benchmark/Rout/create_benchmark_dataframe.Rout"
    output:
        rds = "benchmark/Rdata/benchmark_dataframe.rds"
    shell:
        '''
        {Rbin} "--args outrds='{output.rds}' fileglob='{params.fileglob}'" {input.script} {log.rout}
        '''

rule create_benchmark_dataframe_scaling:
    input: 
        expand("{dataset}/output/mutscan_scaling/digestFastqs_{sample}_{ncores}cores.rds", dataset = ["Li_tRNA_sel30"], sample = ["inputA"], ncores = [1, 2, 5, 10, 15, 20]),
        script = "scripts/create_benchmark_dataframe_scaling.R",
    envmodules:
        f"{Rmod}"
    params:
        fileglob = "*/output/benchmark/*cores.txt"
    log: 
        rout = "benchmark/Rout/create_benchmark_dataframe_scaling.Rout"
    output:
        rds = "benchmark/Rdata/benchmark_dataframe_scaling.rds"
    shell:
        '''
        {Rbin} "--args outrds='{output.rds}' fileglob='{params.fileglob}'" {input.script} {log.rout}
        '''

rule create_benchmark_dataframe_scaling_nreads:
    input: 
        expand("{dataset}/output/mutscan_scaling/digestFastqs_{sample}_{nreads}reads.rds", dataset = ["Li_tRNA_sel30"], sample = ["inputA"], nreads = [100000,500000,1000000,5000000,10000000,50000000,100000000,250000000,391420100]),
        script = "scripts/create_benchmark_dataframe_scaling_nreads.R",
    envmodules:
        f"{Rmod}"
    params:
        fileglob = "*/output/benchmark/*reads.txt"
    log: 
        rout = "benchmark/Rout/create_benchmark_dataframe_scaling_nreads.Rout"
    output:
        rds = "benchmark/Rdata/benchmark_dataframe_scaling_nreads.rds"
    shell:
        '''
        {Rbin} "--args outrds='{output.rds}' fileglob='{params.fileglob}'" {input.script} {log.rout}
        '''

rule plot_benchmark_results:
    input: 
        rds = "benchmark/Rdata/benchmark_dataframe.rds",
        script = "scripts/plot_benchmark_results.R"
    envmodules:
        f"{Rmod}"
    log: 
        rout = "benchmark/Rout/plot_benchmark_results.Rout"
    output:
        pdf = "benchmark/pdf/benchmark_results.pdf"
    shell:
        '''
        {Rbin} "--args inrds='{input.rds}' outpdf='{output.pdf}'" {input.script} {log.rout}
        '''

rule plot_benchmark_results_scaling:
    input: 
        rds = "benchmark/Rdata/benchmark_dataframe_scaling.rds",
        script = "scripts/plot_benchmark_results_scaling.R"
    envmodules:
        f"{Rmod}"
    log: 
        rout = "benchmark/Rout/plot_benchmark_results_scaling.Rout"
    output:
        pdf = "benchmark/pdf/benchmark_results_scaling.pdf"
    shell:
        '''
        {Rbin} "--args inrds='{input.rds}' outpdf='{output.pdf}'" {input.script} {log.rout}
        '''

rule plot_benchmark_results_scaling_nreads:
    input: 
        rds = "benchmark/Rdata/benchmark_dataframe_scaling_nreads.rds",
        script = "scripts/plot_benchmark_results_scaling_nreads.R"
    envmodules:
        f"{Rmod}"
    log: 
        rout = "benchmark/Rout/plot_benchmark_results_scaling_nreads.Rout"
    output:
        pdf = "benchmark/pdf/benchmark_results_scaling_nreads.pdf"
    shell:
        '''
        {Rbin} "--args inrds='{input.rds}' outpdf='{output.pdf}'" {input.script} {log.rout}
        '''

rule plot_nullcomparison:
    input: 
        expand("{dataset}/output/mutscan/{dataset}_mutscan_nullcomparison.rds", dataset = datasets),
        script = "scripts/plot_nullcomparison.R",
    envmodules:
        f"{Rmod}"
    params:
        fileglob = "*/output/mutscan/*_mutscan_nullcomparison.rds"
    log: 
        rout = "benchmark/Rout/plot_nullcomparison.Rout"
    output:
        pdf = "benchmark/pdf/nullcomparison.pdf"
    shell:
        '''
        {Rbin} "--args outpdf='{output.pdf}' fileglob='{params.fileglob}'" {input.script} {log.rout}
        '''

## ------------------------------------------------------------------------------------ ##
## Run Diss_FOS_JUN case study
## ------------------------------------------------------------------------------------ ##
rule run_case_study:
    input:
        "{dataset}/get_data_checkpoint.txt",
        rmd = "{dataset}/CaseStudy/{dataset}_CaseStudy.Rmd"
    output:
        html = "{dataset}/CaseStudy/{dataset}_CaseStudy.html"
    envmodules:
        f"{Rmod}"
    shell:
        '''
        cd $(dirname {input.rmd}) && \
		Rscript -e "rmarkdown::render('$(basename {input.rmd})', clean = FALSE)"
        '''

## ------------------------------------------------------------------------------------ ##
## Generate figures for paper + supplement
## ------------------------------------------------------------------------------------ ##
rule make_paper_figures:
    input:
        "paper/figures/mutscan-overview/mutscan-overview.003.png",
        "Diss_FOS_JUN/CaseStudy/Diss_FOS_JUN_CaseStudy.html",
        "benchmark/pdf/nullcomparison.pdf",
        "benchmark/pdf/benchmark_results.pdf",
        "Diss_FOS/plots/compare_counts_CISOU1_pairs_allvariants.png",
        "Diss_FOS_JUN/plots/compare_counts_TRANSOU1_pairs_allvariants.png",
        "Bolognesi_TDP43_290_331/plots/compare_counts_bTDP1OU7_pairs_allvariants.png",
        "Li_tRNA_sel30/plots/compare_counts_sel30A_pairs_allvariants.png",
        expand("{dataset}/plots/compare_output_pairs.png", dataset = datasets),
        script = "scripts/assemble_figures_for_paper.R"
    envmodules:
        f"{Rmoddevel}",
        f"{popplermod}"
    output:
        expand("paper/figures/Fig{fig}.{ext}", fig=[1,2,3,4,5,6], ext=["png", "pdf"])
    log:
        rout = "paper/figures/assemble_figures_for_paper.Rout"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        '''
        {Rbin} "--args outdir='{params.outdir}'" {input.script} {log.rout}
        '''

rule compile_supplement:
    input:
        "benchmark/pdf/benchmark_results_scaling.pdf",
        "benchmark/pdf/benchmark_results_scaling_nreads.pdf",
        expand("{dataset}/plots/compare_output_upset.png", dataset = datasets),
        "Diss_FOS/plots/compare_counts_CISOU1_pairs_allvariants.png",
        "Diss_FOS_JUN/plots/compare_counts_TRANSOU1_pairs_allvariants.png",
        "Bolognesi_TDP43_290_331/plots/compare_counts_bTDP1OU7_pairs_allvariants.png",
        "Li_tRNA_sel30/plots/compare_counts_sel30A_pairs_allvariants.png",
        "Diss_FOS/plots/compare_output_pairs.png",
        "Bolognesi_TDP43_290_331/plots/compare_output_pairs.png",
        tex = "paper/supplement/mutscan_supplement.tex"
    output:
        pdf = "paper/supplement/mutscan_supplement.pdf"
    params:
        indir = lambda wildcards, input: os.path.dirname(input["tex"])
    shell:
        '''
        cd {params.indir} && \
        pdflatex $(basename {input.tex})
        '''

rule listpackages:
    input:
        expand("{dataset}/Rout/mutscan_nullcomparison.Rout", dataset = datasets),
        expand("{dataset}/Rout/plot_comparison_counts_scores.Rout", dataset = datasets),
        expand("{dataset}/Rout/plot_counts_singlesample_{sample}.Rout", zip, dataset = datasetsVec, sample = samplesVec),
        expand("benchmark/Rout/plot_benchmark_results{ext}.Rout", ext = ["", "_scaling", "_scaling_nreads"]),
        "benchmark/Rout/plot_nullcomparison.Rout",
        "Diss_FOS_JUN/CaseStudy/Diss_FOS_JUN_CaseStudy.html",
        script = "scripts/list_packages.R",
    output:
        outtxt = "paper/R_package_versions.txt",
        outrds = "paper/R_package_versions.rds",
    log:
        rout = "Rout/list_packages.Rout"
    params:
        Routdirs = ",".join([dataset + "/Rout" for dataset in datasets] + ["benchmark/Rout"]),
        Reportdirs = "Diss_FOS_JUN/CaseStudy",
    envmodules:
        f"{Rmod}"
    shell:
        '''
        {Rbin} "--args Routdirs='{params.Routdirs}' Reportdirs='{params.Reportdirs}' outtxt='{output.outtxt}' outrds='{output.outrds}'" {input.script} {log.rout}
        '''
