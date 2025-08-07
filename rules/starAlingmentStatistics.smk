"""

creating tab delimited file of alingment statistics from STAR alingment 

"""

"sRNA seq log files"

rule STAR_alingment_stats_strains_srna:
    input:
        pacBio = expand(config['resultDir']  + f"/STAR_srna_strain_wise/{{sid}}/{{sid}}-{{tp}}.{{rep}}/Log.final.out",  sid=strains, tp=tps, rep=reps ),
        scriptDir = config['scriptDir']
    output:
        log_strain_stat = config['resultDir']  + "/log_stat/strain_srna_mapping_log_stat.tab",
    threads: 1
    priority:10
    resources:
        mem = 64000,
        queue  = 'normal',
    log:
        config['logDir']  + "/log_stat_srna_STAR_mapping.log"
    shell: """
        awk -f {input.scriptDir}/mergeLogFinal.awk  {input.pacBio} | awk -f {input.scriptDir}/sh/pivot_csv.awk > {output.log_strain_stat} 2> {log}
       
    """

#input
#grchm39 = expand(config['resultDir']  + f"/STAR_srna_GRCm39_106/{{sid}}-{{tp}}.{{rep}}/Log.final.out",  sid=strains, tp=tps, rep=reps ),
#       old_genomes = expand(config['resultDir']  + f"/STAR_srna_strain_wise_old_genomes/{{sid}}/{{sid}}-{{tp}}.{{rep}}/Log.final.out",  sid=strains, tp=tps, rep=reps ),
#       grchm38 = expand(config['resultDir']  + f"/STAR_srna_GRCm38_68/{{sid}}-{{tp}}.{{rep}}/Log.final.out",  sid=strains, tp=tps, rep=reps ),
#output
#  log_grc39_stat = config['resultDir']  + "/log_stat/grc39_srna_mapping_log_stat.tab",
#     log_old_genomes = config['resultDir']  + "/log_stat/old_genomes_srna_mapping_log_stat.tab",
#     log_grc38_stat = config['resultDir']  + "/log_stat/grc38_srna_mapping_log_stat.tab",
#commnad
# awk -f {input.scriptDir}/mergeLogFinal.awk  {input.grchm39} | awk -f {input.scriptDir}/sh/pivot_csv.awk > {output.log_grc39_stat} 2> {log}
#       awk -f {input.scriptDir}/mergeLogFinal.awk  {input.grchm38} | awk -f {input.scriptDir}/sh/pivot_csv.awk > {output.log_grc38_stat} 2> {log}
#       awk -f {input.scriptDir}/mergeLogFinal.awk  {input.old_genomes} | awk -f {input.scriptDir}/sh/pivot_csv.awk > {output.log_old_genomes} 2> {log}



"""
RNA mapping statistics 

"""


rule STAR_mapping_stats_strains_rna:
    input:
        log_files = expand(config['resultDir']  + f"/STAR_rna_strain_wise/{{sid}}-{{tp}}.{{rep}}/Log.final.out",  sid=strains, tp=tps, rep=reps ),
        scriptDir = config['scriptDir']
    output:
        log_strain_stat = config['resultDir']  + "/log_stat/strain_rna_mapping_log_stat.tab",
    threads: 1
    priority:10
    resources:
        mem = 64000,
        queue  = 'normal',
    log:
        config['logDir']  + f"/log_stat_"
    shell: """
        awk -f {input.scriptDir}/mergeLogFinal.awk  {input.log_files}  | awk -f {input.scriptDir}/sh/pivot_csv.awk > {output.log_strain_stat} 2> {log}_rna_pacBio_STAR_mapping.log
        """

rule STAR_mapping_stats_Old_genome_strains_rna:
    input:
        log_files = expand(config['resultDir']  + f"/STAR_rna_old_genome_strain_wise/{{sid}}-{{tp}}.{{rep}}/Log.final.out",  sid=strains, tp=tps, rep=reps ),
        scriptDir = config['scriptDir']
    output:
        log_strain_stat = config['resultDir']  + "/log_stat/Old_Genome_strain_rna_mapping_log_stat.tab",
    threads: 1
    priority:10
    resources:
        mem = 64000,
        queue  = 'normal',
    log:
        config['logDir']  + f"/log_stat_"
    shell: """
        awk -f {input.scriptDir}/mergeLogFinal.awk  {input.log_files}  | awk -f {input.scriptDir}/sh/pivot_csv.awk > {output.log_strain_stat} 2> {log}_rna_old_genome_STAR_mapping.log
        """





