rule assembly_mmseqs_evaluation_megahit:
    input:
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.megahit.{{config}}.fasta',
    output:
        tsv=f'{config["rdir"]}/assembly-mmseqs-nr/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.megahit.{{config}}.tsv',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    params:
        seqkit_bin   = config["seqkit_bin"],
        ancientContigs="/net/mimer/mnt/tank/projects2/metagnm_asm/project3/data/results/ancient_headers_all.csv",
        ancientFasta=f'{config["rdir"]}/assembly-mmseqs-nr/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.megahit.{{config}}.ancient.fasta',
        tmp=f'{config["rdir"]}/assembly-mmseqs-nr/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.megahit.{{config}}_tmp',
        mmseqs_bin=config["mmseqs_bin"],
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-mmseqs-nr",
        reference=config["nr"],
    threads: config["threads_32"]
    log:
        mmseqs_log=f'{config["rdir"]}/assembly-mmseqs-nr/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.megahit.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-mmseqs-nr/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.megahit.{{config}}.bmk'
    message:
        """--- assemblies vs. reference (megahit) """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        N=$(grep -c '>' {input.contigs} || [[ $? == 1 ]])

        if [[ ${{N}} -eq 0 ]]; then
            touch {output.tsv}
            exit 0
        fi

        {params.seqkit_bin} grep \
            -f {params.ancientContigs} --quiet \
            {input.contigs} \
            > {params.ancientFasta}

        {params.mmseqs_bin} easy-search \
            {params.ancientFasta} \
            {params.reference} \
            {output.tsv} \
            {params.tmp} \
            --max-seq-len 1000000 --search-type 4 --cov-mode 1 --min-aln-len 300 \
            --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov,qorfstart,qorfend,torfstart,torfend \
            --threads {threads} --split-memory-limit 200G >> {log.mmseqs_log} 2>&1

        #format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov,qframe,tframe,qorfstart,qorfend,torfstart,torfend


        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
