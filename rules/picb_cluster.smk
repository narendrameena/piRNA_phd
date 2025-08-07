rule STAR_index_strains_pacBio_genome_picb:
    localrule: True
    input:
        fasta=config['resources'] + "/PICB/refFasta/{sid}_chromosomes_MT.fasta",
        gtf=lambda wildcards: config['resources'] + f"/PICB/annotation/{wildcards.sid}_v3.3.gff3"
    output:
        directory(config['resources']  + "/PICB/index/{sid}")
    threads: 4
    params:
        extra = "--limitGenomeGenerateRAM  96000000000",
    message:
        "****STAR index****"
    priority:2
    resources:
        mem = 200000,
        queue  = '2004*',
    log:
        config['logDir']  + "/PICB/index/{sid}.log",
    wrapper:
        "v5.8.2/bio/star/index"


rule STAR_alingment_srna_strain_picb:
    localrule: True
    input:
        fq1 =  config['resultDir']  + "/cutadapt_se/{sid}-{tp}.{rep}.fastq",
        idx = config['resources'] + "/PICB/index/{sid}",
        #fq1 = config['resultDir']  + "/STAR_srna_strain_wise/tmp_fasta/{sid}-{tp}.{rep}.fasta"
    output:
        aln=temp(config['resultDir']  + "/STAR_srna_strain_wise_picb/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam"),
        log_final=config['resultDir']  + "/STAR_srna_strain_wise_picb/{sid}/{sid}-{tp}.{rep}/Log.final.out"
    threads: 4
    priority:3
    log:
        config['logDir']  + "/STAR_srna_strain_wise_picb/{sid}/{sid}-{tp}.{rep}.log"
    resources:
        mem = 200000,
        queue  = '2004*',
        tmpdir= config['tmpDir']
    params:
        idx = config['resources'] + "/PICB/index/{sid}",
        extra = config['params']['srna']['STAR'] 
    wrapper:
        "v5.8.2/bio/star/align"



rule STAR_strains_srna_bam_index_picb:
    localrule: True
    input:
        config['resultDir']  + "/STAR_srna_strain_wise_picb/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
    output:
        temp(config['resultDir']  + "/STAR_srna_strain_wise_picb/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam.bai"),
    log:
        config['logDir']  + "/STAR_srna_strain_wise_idx_picb/{sid}-{tp}.{rep}.log"
    resources:
        mem = 200000,
        queue  = '2004*',
    params:
        extra="",
    threads: 1
    priority:4
    shell:
        "samtools index {input}   > {log} 2>&1"

# "v1.7.0/bio/samtools/index"

rule run_picb_analysis:
    localrule: True
    input:
        fasta = config['resources'] + "/PICB/refFasta/{sid}_chromosomes_MT.fasta",
        bam = config['resultDir'] + "/STAR_srna_strain_wise_picb/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        bai = config['resultDir'] + "/STAR_srna_strain_wise_picb/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam.bai",
    output:
        excel = config['resultDir'] + "/picb_result/{sid}/{sid}-{tp}.{rep}.xlsx",
    log:
        config['logDir'] + "/picb_cluster/{sid}-{tp}.{rep}.log"
    resources:
        mem = 800000,
        queue = '2004*',
    params:
        script = config['scriptDir'] + "/R/picb_script.R",
        target_size_gb = 40,  # Target size in GB for downsampling
        downsampled_dir = config['resultDir'] + "/STAR_srna_strain_wise_picb/{sid}/{sid}-{tp}.{rep}"
    priority: 12
    shell:
        """
        set -euo pipefail

        # Define the directory and file paths for the downsampled BAM
        downsampled_dir="{params.downsampled_dir}"
        downsampled_bam="$downsampled_dir/Aligned.sortedByCoord.out.downsampled.bam"
        downsampled_bai="$downsampled_bam.bai"

        # Get the size of the original BAM file in bytes
        original_size=$(stat -c%s "{input.bam}")

        # Convert original size to GB
        original_size_gb=$(echo "scale=2; $original_size / (1024^3)" | bc)

        echo "Original BAM size: $original_size_gb GB" >&2

        # Check if the original BAM size is greater than target_size_gb
        if (( $(echo "$original_size_gb > {params.target_size_gb}" | bc -l) )); then
            echo "Original BAM size is greater than {params.target_size_gb} GB. Downsampling..." >&2

            # Calculate the target size in bytes
            target_size_bytes=$(echo "{params.target_size_gb} * 1024 * 1024 * 1024" | bc)

            # Count the total number of reads in the original BAM file
            echo "Counting total reads in the original BAM file..." >&2
            total_reads=$(samtools view -c "{input.bam}")

            # Calculate the average bytes per read
            avg_bytes_per_read=$(echo "scale=6; $original_size / $total_reads" | bc -l)

            # Calculate the target number of reads to get approximately 40GB
            target_reads=$(echo "scale=0; $target_size_bytes / $avg_bytes_per_read" | bc)

            # Calculate the fraction to keep
            fraction=$(echo "scale=6; $target_reads / $total_reads" | bc -l | sed 's/^\./0./')

            # Ensure the fraction is between reasonable limits
            min_fraction=0.0001  # Minimum fraction to avoid too small files
            if (( $(echo "$fraction > 1" | bc -l) )); then
                fraction=1
            elif (( $(echo "$fraction < $min_fraction" | bc -l) )); then
                echo "Calculated fraction is too small ($fraction). Skipping downsampling." >&2
                fraction=1
            fi

            echo "Total reads: $total_reads" >&2
            echo "Average bytes per read: $avg_bytes_per_read" >&2
            echo "Target reads: $target_reads" >&2
            echo "Fraction to keep: $fraction" >&2

            # Prepare the fraction for samtools (remove leading zeros and dot)
            fraction_for_samtools=$(echo "$fraction" | sed 's/^0*\.//')

            # Downsample the BAM file
            echo "Downsampling BAM file..." >&2
            samtools view -s 42."$fraction_for_samtools" -b "{input.bam}" > "$downsampled_bam"

            # Verify the downsampled BAM file is valid
            if ! samtools quickcheck "$downsampled_bam"; then
                echo "Downsampled BAM file is corrupted or invalid. Using original BAM file." >&2
                bam_file="{input.bam}"
                rm -f "$downsampled_bam"
            else
                # Index the downsampled BAM file if required by the R script
                samtools index "$downsampled_bam"

                # Verify the size of the downsampled BAM file
                downsampled_size=$(stat -c%s "$downsampled_bam")
                downsampled_size_gb=$(echo "scale=2; $downsampled_size / (1024^3)" | bc)
                echo "Downsampled BAM size: $downsampled_size_gb GB" >&2

                bam_file="$downsampled_bam"
            fi
        else
            echo "Original BAM size is less than or equal to {params.target_size_gb} GB. Skipping downsampling." >&2

            # Use the original BAM file
            bam_file="{input.bam}"
        fi

        # Run the R script with the appropriate BAM file
        echo "Running R script..." >&2
        Rscript "{params.script}" "{input.fasta}" "$bam_file" "{output.excel}" > "{log}" 2>&1

        # Clean up the downsampled BAM file if it exists
        if [ "$bam_file" = "$downsampled_bam" ]; then
            rm -f "$downsampled_bam" "$downsampled_bai"
        fi

        # Touch the output file to signal completion to Snakemake
        touch "{output.excel}"
        """

rule picb_combine_replicates:
    localrule: True
    input:
        xlsx=[
            config['resultDir'] + "/picb_result/{sid}/{sid}-{tp}.1.xlsx",
            config['resultDir'] + "/picb_result/{sid}/{sid}-{tp}.2.xlsx",
            config['resultDir'] + "/picb_result/{sid}/{sid}-{tp}.3.xlsx"
        ],
        bam=[
            config['resultDir'] + "/STAR_srna_strain_wise_picb/{sid}/{sid}-{tp}.1/Aligned.sortedByCoord.out.bam",
            config['resultDir'] + "/STAR_srna_strain_wise_picb/{sid}/{sid}-{tp}.2/Aligned.sortedByCoord.out.bam",
            config['resultDir'] + "/STAR_srna_strain_wise_picb/{sid}/{sid}-{tp}.3/Aligned.sortedByCoord.out.bam"
        ]
    output:
        config['resultDir'] + "/picb_result/{sid}/{sid}-{tp}.combined.xlsx"
    log:
        config['logDir'] + "/picb_combine/{sid}-{tp}.log"
    params:
        script = config['scriptDir'] + "/R/picb_combine_script.R",
        ref = config['resources'] + "/PICB/refFasta/{sid}_chromosomes_MT.fasta"
    shell:
        """
        Rscript "{params.script}" \
            --xlsxs {input.xlsx[0]} {input.xlsx[1]} {input.xlsx[2]} \
            --bams {input.bam[0]} {input.bam[1]} {input.bam[2]} \
            --ref "{params.ref}" \
            --out "{output}" \
            > "{log}" 2>&1
        """


rule run_picb_analysis_ex:
    input:
        fasta = config['resources'] + "/PICB/refFasta/{sid}_chromosomes_MT.fasta",
        bam = config['resultDir'] + "/STAR_srna_strain_wise_picb/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam",
        bai = config['resultDir'] + "/STAR_srna_strain_wise_picb/{sid}/{sid}-{tp}.{rep}/Aligned.sortedByCoord.out.bam.bai",
    output:
        excel = config['resultDir'] + "/picb_result/{sid}/{sid}-{tp}.{rep}.xlsx_ex",
    log:
        config['logDir'] + "/picb_cluster/{sid}-{tp}.{rep}.log"
    resources:
        mem = 800000,
        queue = '2004*',
    params:
        script = config['scriptDir'] + "/R/picb_script.R",
    priority: 12
    shell:
        """
        Rscript {params.script} {input.fasta} {input.bam} {output.excel} > {log} 2>&1
        """


rule picb_convert_to_bed:
    localrule: True
    input:
        excel=config['resultDir'] + r"/picb_result/{sid}/{sid}-{tp}.{rep}.xlsx"
    output:
        bed=config['resultDir'] + "/picb_result/{sid}/{sid}-{tp}.{rep}.bed"
    log:
        config['logDir'] + "/convert_picb_to_bed/{sid}-{tp}.{rep}.log"
    resources:
        mem=4000,
        queue='2004*'
    params:
        script=config['scriptDir'] + "/python/convert_picb_excel_to_bed6.py"
    priority: 13
    conda:
        config['envsDir'] +"/bedannotation.yaml"  # Define this in your config file
    shell:
        """
        python {params.script} {input.excel} {output.bed} > {log} 2>&1
        """

rule fix_chromosome_names:
    localrule: True
    input:
        bed=config['resultDir'] + "/picb_result/{sid}/{sid}-{tp}.{rep}.bed"
    output:
        fixed_bed=config['resultDir'] + "/picb_result/{sid}/{sid}-{tp}.{rep}_fixed.bed"
    log:
        config['logDir'] + "/fix_chromosome_names/{sid}-{tp}.{rep}.log"
    params:
        strain="{sid}"
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} {{
            if ($1 !~ /^{{params.strain}}#1#/) $1="{{params.strain}}#1#" $1;
            gsub(/^{{params.strain}}#1#chrM$/, "{{params.strain}}#1#MT", $1);
            print
        }}' {input.bed} > {output.fixed_bed} 2> {log}
        """

# Rule to liftOver the BED file
rule liftover_strains_to_GRCm39:
    input:
        bed=config['resultDir'] + "/picb_result/{sid}/{sid}-{tp}.{rep}_fixed.bed",
        chain=config['resultDir'] + "/minimap2_chain_alignments/{sid}_to_GRCm39.chain"
    output:
        bed=config['resultDir'] + "/picb_result/{sid}/{sid}-{tp}.{rep}_lifted.bed",
        unmapped=config['resultDir'] + "/picb_result/{sid}/{sid}-{tp}.{rep}_unmapped.bed"
    log:
        config['logDir'] + "/liftover/{sid}-{tp}.{rep}.log"
    params:
        extra="-fudgeThick"
    resources:
        mem_mb=4000,
        queue="normal"
    conda:
        config['envsDir'] + "/liftover.yaml"
    shell:
        """
        liftOver {params.extra} {input.bed} {input.chain} {output.bed} {output.unmapped} > {log} 2>&1
        """

rule convert_picb_to_gff:
    localrule: True
    input:
        excel=config['resultDir'] + "/picb_result/{sid}/{sid}-{tp}.{rep}.xlsx"
    output:
        gff=config['resultDir'] + "/picb_result/{sid}/{sid}-{tp}.{rep}.gff"
    log:
        config['logDir'] + "/convert_to_gff/{sid}-{tp}.{rep}.log"
    resources:
        mem=4000,
        queue='2004*'
    params:
        script=config['scriptDir'] + "/python/convert_picb_excel_to_gff.py"
    priority: 10
    conda:
        config['envsDir'] +"/bedannotation.yaml"  # Path to the Conda environment
    shell:
        """
        python {params.script} {input.excel} {output.gff} > {log} 2>&1
        """


rule liftoff_conversion:
    input:
        gff="path/to/input/{sid}-{tp}.{rep}.gff",  # Input GFF file
        reference_fasta="path/to/reference/{sid}_reference.fasta",  # Reference genome FASTA
        target_fasta="path/to/target/{sid}_target.fasta"  # Target genome FASTA
    output:
        lifted_gff="path/to/liftoff_result/{sid}-{tp}.{rep}_lifted.gff"
    log:
        "path/to/logs/liftoff/{sid}-{tp}.{rep}.log"
    params:
        output_prefix="path/to/liftoff_result/{sid}-{tp}.{rep}",  # Prefix for Liftoff output files
    conda:
         config['envsDir'] +"/iftoff_env.yml"  # Conda environment for Liftoff
    shell:
        """
        liftoff {input.target_fasta} {input.reference_fasta}  -g {input.gff}   -o {output.lifted_gff}   --chrom > {log} 2>&1
        """