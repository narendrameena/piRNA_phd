

### refrance based trascriptome assemblies 

rule stringtie_transcriptome_rna_pacBio:
    localrule: True
    input:
        bam = config['resultDir']  + "/STAR_rna_strain_wise/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        gtf = config['resultDir'] + "/strains_genome_annotation/{sid}.gff3",
        sbam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
    output:
        gtf = protected(config['resultDir']  + "/refrancePrecursors/pacBio/{sid}/{sid}-{tp}.{rep}.transcripts.gtf"),
        #ctab = directory(config['resultDir']  + "/refrancePrecursors/pacBio/{sid}/{sid}-{tp}.{rep}.transcripts_ctab"),
        geneAbu = config['resultDir']  + "/refrancePrecursors/pacBio/{sid}/{sid}-{tp}.{rep}.transcripts.tab",
        covGtf = config['resultDir']  + "/refrancePrecursors/pacBio/{sid}/{sid}-{tp}.{rep}.transcripts.coverage.gtf"
    threads: 4
    conda: config['envsDir'] + "/refrancebasedprecursors.yaml"
    params:
        " -m 500 -l {sid}-{tp} "
    resources:
        mem = 200000,
        queue  = 'normal',
        tmpdir= config['tmpDir']
    priority:10
    log:
        config['logDir']  + "/cufflinks_rna_strain_wise/{sid}/{sid}-{tp}.{rep}.log"
    shell:
        """
        stringtie --mix  {input.sbam} {input.bam} -o {output.gtf}   -A  {output.geneAbu} -C  {output.covGtf}  {params} -p {threads} -G {input.gtf} > {log}
        """

rule stringtie_feature_counts_pacBio:
    localrule: True
    input:
        sam = config['resultDir']  + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam" ,# list of sam or bam files
        annotation = config['resultDir']  + "/refrancePrecursors/pacBio/{sid}/{sid}-{tp}.{rep}.transcripts.gtf",
    output:
        protected(multiext(config['resultDir']  + "/stringtieFeatureCounts/pacBio/{sid}/{sid}-{tp}.{rep}",
                 ".featureCounts",
                 ".featureCounts.summary",
                 ".featureCounts.jcounts"))
    threads:       2
    priority:10
    resources:
        mem = 20000,
        queue  = 'normal',
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra = config['params']['stringtie']['featureCounts'] 
    log:
        config['logDir']  +"/stringtieFeatureCounts/pacBio/list3/{sid}/{sid}-{tp}.{rep}.log"
    wrapper:
        "v5.8.2/bio/subread/featurecounts"

rule refrance_precursors_calculate_RPM_RPKM:
    localrule: True
    input:
        aligned_bam = config['resultDir'] + "/STAR_srna_strain_wise/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        annotation_gtf = config['resultDir'] + "/refrancePrecursors/pacBio/{sid}/{sid}-{tp}.{rep}.transcripts.gtf",
        script = config['scriptDir'] + "/sh/calculate_RPM_RPKM.sh",
        countF = config['resultDir']  + "/stringtieFeatureCounts/pacBio/{sid}/{sid}-{tp}.{rep}.featureCounts"
    output:
        output_file = config['resultDir'] + "/refrancePrecurosrsRPMRPKMCounts/pacBio/{sid}/{sid}-{tp}.{rep}.rpm_rpkm.txt",
    threads: 1
    conda: config['envsDir'] + "/refrancebasedprecursors.yaml"
    resources:
        mem = 200000,
        queue = 'normal',
        tmpdir = config['tmpDir']
    priority: 10
    log:
        config['logDir'] + "/rpm_rpkm_refrance_precursors_strain_wise/{sid}/{sid}-{tp}.{rep}.log"
    shell:
        """
        bash {input.script} {input.aligned_bam} {input.annotation_gtf} {input.countF} {output.output_file} > {log}
        """


rule filtred_refrance_precursors_rpm_100:
    localrule: True
    input:
        rpm = config['resultDir'] + "/refrancePrecurosrsRPMRPKMCounts/pacBio/{sid}/{sid}-{tp}.{rep}.rpm_rpkm.txt",
        gtf = config['resultDir']  + "/refrancePrecursors/pacBio/{sid}/{sid}-{tp}.{rep}.transcripts.gtf",
        script = config['scriptDir'] + "/perl/gtftobed.pl",
        script2 = config['scriptDir'] + "/sh/refrancefilter100tobed.sh",
    output:
        bed = config['resultDir'] + "/refrancePrecursors/pacBio/{sid}/{sid}-{tp}.{rep}.transcripts.bed",
        finalbed = config['resultDir'] + "/filtred_refranc_Precurosrs/pacBio/{sid}/{sid}-{tp}.{rep}.bed",
    threads: 1
    conda: config['envsDir'] + "/refrancebasedprecursors.yaml"
    resources:
        mem = 200000,
        queue = 'normal',
        tmpdir = config['tmpDir']
    priority: 10
    log:
        config['logDir'] + "/filtred_100_rpm_refrance_precursors_strain_wise/{sid}/{sid}-{tp}.{rep}.log"
    shell:
        """
        perl {input.script}  {input.gtf} >  {output.bed} 
        bash {input.script2} {input.rpm} {output.bed}  {output.finalbed} > {log}
        """


rule merge_precurosors:
    localrule: True
    input:
        refrance = config['resultDir'] + "/filtred_refranc_Precurosrs/pacBio/{sid}/{sid}-{tp}.{rep}.bed",
        trinity = config['resultDir']  + "/filter_precursors_bed/{sid}/{sid}-{tp}.{rep}.bed",
    output:
        first= config['resultDir'] + "/filtred_final_Precurosrs/pacBio/{sid}/{sid}-{tp}.{rep}_merge_info.bed",
        second = config['resultDir'] + "/filtred_final_Precurosrs/pacBio/{sid}/{sid}-{tp}.{rep}.bed",
    threads: 1
    conda: config['envsDir'] + "/refrancebasedprecursors.yaml"
    resources:
        mem = 200000,
        queue = 'normal',
        tmpdir = config['tmpDir']
    priority: 10
    log:
        config['logDir'] + "/filtred_precursors/{sid}/{sid}-{tp}.{rep}.log"
    shell:
        """
        cat {input.refrance} {input.trinity} | sort -k1,1 -k2,2n | bedtools merge -i -  -c 1,2,3,4,5,6,7,8,9,10,11,12 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > {output.first}
        cat {input.refrance} {input.trinity} | sort -k1,1 -k2,2n | bedtools merge -i -  -s -c 6 -o distinct  > {output.second}
        """

