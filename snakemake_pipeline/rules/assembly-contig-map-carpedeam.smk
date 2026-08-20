rule contig_mapping_carpedeam:
    input:
        contigs = f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.carpedeam.{{config}}.fasta',
        reads   = f'{config["rdir"]}/assembly-input-2/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.fastq.gz',
    output:
        bam_sorted       = f'{config["rdir"]}/assembly-contig-map/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mapping.carpedeam.{{config}}.sorted.bam',
    threads: config["threads_16"]
    params:
        assm_input_1 = "{assm_input_1}",
        assm_input_2 = "{assm_input_2}",
        # TSV file with filtered contigs (see your original TSV)
        #tsv          = f'{config["rdir"]}/assembly-mmseqs-filter/filter_uniref90/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.carpedeam.{{config}}_uniref90.tsv',
        seqkit_bin   = config["seqkit_bin"],
        bwa_bin      = config["bwa_bin"],        # new: path to BWA binary
        samtools_bin = config["samtools_bin"],
        wdir         = config["wdir"],
        rdir         = config["rdir"] + "/assembly-contig-map",
        # intermediate filenames for SAM and unsorted BAM
        sam = f'{config["rdir"]}/assembly-contig-map/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mapping.carpedeam.{{config}}.sam',
        bam = f'{config["rdir"]}/assembly-contig-map/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mapping.carpedeam.{{config}}.bam',
    log:
        mapping_log = f'{config["rdir"]}/assembly-contig-map/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mapping.carpedeam.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-contig-map/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mapping.carpedeam.{{config}}.bmk'
    message:
        """--- Mapping reads against filtered contigs using BWA."""
    shell:
        r"""
        set -euo pipefail
        
        # If TSV has only one line (only a header), just create dummy outputs and exit
        if [[ $(wc -l < {input.contigs}) -le 1 ]]; then
            touch {output.bam_sorted}
            exit 0
        fi

        # # Extract contig IDs (skip header) and filter the original contigs file
        # cut -f 1 params.tsv | tail -n +2 > contig_ids.txt
        # rm -f contig_ids.txt

        # Build a BWA index on the filtered contigs
        echo "Building BWA index on filtered contigs..." >> {log.mapping_log}
        {params.bwa_bin} index {input.contigs} >> {log.mapping_log} 2>&1

        # Run BWA aln to generate a .sai file (use a temporary file)
        echo "Running BWA aln..." >> {log.mapping_log}
        sai=$(mktemp)
        {params.bwa_bin} aln -l 1024 -o 1 -n 0.1 -t {threads} {input.contigs} {input.reads} > $sai 2>> {log.mapping_log}

        # Generate SAM file using BWA samse
        echo "Generating SAM file with BWA samse..." >> {log.mapping_log}
        {params.bwa_bin} samse {input.contigs} $sai {input.reads} > {params.sam} 2>> {log.mapping_log}

        # Convert SAM to BAM
        echo "Converting SAM to BAM..." >> {log.mapping_log}
        {params.samtools_bin} view -bS {params.sam} -o {params.bam}

        # Sort the BAM file
        echo "Sorting BAM file..." >> {log.mapping_log}
        tmp_sorted=$(mktemp)
        {params.samtools_bin} sort --threads {threads} {params.bam} -o $tmp_sorted

        # Add MD flag using samtools calmd (the final sorted BAM is our output)
        echo "Adding MD flag to BAM file..." >> {log.mapping_log}
        {params.samtools_bin} calmd --threads {threads} -b $tmp_sorted {input.contigs} > {output.bam_sorted}

        # Index the final BAM file
        echo "Indexing final BAM file..." >> {log.mapping_log}
        {params.samtools_bin} index {output.bam_sorted}

        # Clean up temporary files
        rm -f $sai {params.sam} {params.bam} $tmp_sorted
        """
