
"""
Plots
"""



rule raw_12_5dpp:
    output:
        protected(config['workdir']  + "analysis/seqCount/raw/InitialPheatMap_zScore_12.5dpp.pdf")
    log:
        stderr=config['logDir']  + "/seqCount/rScript_raw_12.5dpp.err"
    threads: 10
    priority:4
    resources:
        mem = 650000,
        queue = 'normal',
    conda: config['envsDir'] + "/snakmake.yaml"
    shell:
        """
        cd /lustre/scratch123/hgi/projects/mouse_epi/users/naru/inProgres/mice_PiRNA/workflow/scripts/R/seqCount/
        jupyter nbconvert --ExecutePreprocessor.kernel_name='ir' --to notebook --execute raw_DESeq_PCA_12.5dpp.ipynb 
        """


rule raw_16_5dpc:
    output:
        protected(config['workdir']  + "analysis/seqCount/raw/InitialPheatMap_zScore_16.5dpc.pdf")
    log:
        stderr=config['logDir']  + "/seqCount/rScript_raw_16.5dpc.err"
    threads: 10
    priority:4
    resources:
        mem = 650000,
        queue = 'normal',
    conda: config['envsDir'] + "/snakmake.yaml"
    shell:
        """
        cd /lustre/scratch123/hgi/projects/mouse_epi/users/naru/inProgres/mice_PiRNA/workflow/scripts/R/seqCount/
        jupyter nbconvert --ExecutePreprocessor.kernel_name='ir' --to notebook --execute raw_DESeq_PCA_16.5dpc.ipynb
        """
        


rule raw_20_5dpp:
    output:
        protected(config['workdir']  + "analysis/seqCount/raw/InitialPheatMap_zScore_20.5dpp.pdf")
    log:
        stderr=config['logDir']  + "/seqCount/rScript_raw_20.5dpp.err"
    threads: 10
    priority:4
    resources:
        mem = 650000,
        queue = 'normal',
    conda: config['envsDir'] + "/snakmake.yaml"
    shell:
        """
        cd /lustre/scratch123/hgi/projects/mouse_epi/users/naru/inProgres/mice_PiRNA/workflow/scripts/R/seqCount/
        jupyter nbconvert --ExecutePreprocessor.kernel_name='ir' --to notebook --execute raw_DESeq_PCA_20.5dpp.ipynb
        """


rule alinged_12_5dpp:
    output:
        protected(config['workdir']  + "analysis/seqCount/alinged/InitialPheatMap_zScore_12.5dpp.pdf")
    log:
        stderr=config['logDir']  + "/seqCount/rScript_alinged_12.5dpp.err"
    threads: 10
    priority:4
    resources:
        mem = 650000,
        queue = 'normal',
    conda: config['envsDir'] + "/snakmake.yaml"
    shell:
        """
        cd /lustre/scratch123/hgi/projects/mouse_epi/users/naru/inProgres/mice_PiRNA/workflow/scripts/R/seqCount/
        jupyter nbconvert --ExecutePreprocessor.kernel_name='ir' --to notebook --execute alinged_DESeq_PCA_12.5dpp.ipynb
        """


rule alinged_16_5dpc:
    output:
        protected(config['workdir']  + "analysis/seqCount/alinged/InitialPheatMap_zScore_16.5dpc.pdf")
    log:
        stderr=config['logDir']  + "/seqCount/rScript_alinged_16.5dpc.err"
    threads: 10
    priority:4
    resources:
        mem = 650000,
        queue = 'normal',
    conda: config['envsDir'] + "/snakmake.yaml"
    shell:
        """
        cd /lustre/scratch123/hgi/projects/mouse_epi/users/naru/inProgres/mice_PiRNA/workflow/scripts/R/seqCount/
        jupyter nbconvert  --ExecutePreprocessor.kernel_name='ir' --to notebook --execute alinged_DESeq_PCA_16.5dpc.ipynb
        """
        


rule alinged_20_5dpp:
    output:
        protected(config['workdir']  + "analysis/seqCount/alinged/InitialPheatMap_zScore_20.5dpp.pdf")
    log:
        stderr=config['logDir']  + "/seqCount/rScript_alinged_20.5dpp.err"
    threads: 10
    priority:4
    resources:
        mem = 650000,
        queue = 'normal',
    conda: config['envsDir'] + "/snakmake.yaml"
    shell:
        """
        cd /lustre/scratch123/hgi/projects/mouse_epi/users/naru/inProgres/mice_PiRNA/workflow/scripts/R/seqCount/
        jupyter nbconvert  --ExecutePreprocessor.kernel_name='ir' --to notebook --execute alinged_DESeq_PCA_20.5dpp.ipynb
        """

rule alinged:
    output:
        protected(config['workdir']  + "analysis/seqCount/alinged/InitialPheatMap_zScore.pdf")
    log:
        stderr=config['logDir']  + "/seqCount/rScript_aligned.err"
    threads: 10
    priority:4
    resources:
        mem = 650000,
        queue = 'normal',
    conda: config['envsDir'] + "/snakmake.yaml"
    shell:
        """
        cd /lustre/scratch123/hgi/projects/mouse_epi/users/naru/inProgres/mice_PiRNA/workflow/scripts/R/seqCount/
        jupyter nbconvert  --ExecutePreprocessor.kernel_name='ir' --to notebook --execute alinged_DESeq_PCA.ipynb
        """

rule raw:
    output:
        protected(config['workdir']  + "analysis/seqCount/raw/InitialPheatMap_zScore.pdf")
    log:
        stderr=config['logDir']  + "/seqCount/rScript_raw_.err"
    threads: 10
    priority:4
    resources:
        mem = 650000,
        queue = 'normal',
    conda: config['envsDir'] + "/snakmake.yaml"
    shell:
        """
        cd /lustre/scratch123/hgi/projects/mouse_epi/users/naru/inProgres/mice_PiRNA/workflow/scripts/R/seqCount/
        jupyter nbconvert --ExecutePreprocessor.kernel_name='ir' --to notebook --execute raw_DESeq_PCA.ipynb
        """