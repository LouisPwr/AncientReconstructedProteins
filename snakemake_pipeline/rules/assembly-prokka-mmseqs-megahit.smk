rule assembly_prokka_mmseqs_evaluation_megahit:
    input:
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.megahit.{{config}}.fasta',
    output:
        tsv=f'{config["rdir"]}/assembly-prokka-mmseqs-uniref90/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.prokka.megahit.{{config}}.tsv',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    params:
        ancientFasta=f'{config["rdir"]}/assembly-mmseqs-nr/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.megahit.{{config}}.ancient.fasta',
        orfAA=f'{config["rdir"]}/assembly-prokka-mmseqs-uniref90/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.prokka.megahit.{{config}}.ancient.orf.faa',
        orfNuc=f'{config["rdir"]}/assembly-prokka-mmseqs-uniref90/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.prokka.megahit.{{config}}.ancient.orf.fasta',
        tmp=f'{config["rdir"]}/assembly-prokka-mmseqs-uniref90/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.prokka.megahit.{{config}}_tmp',
        mmseqs_bin=config["mmseqs_bin"],
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-prokka-mmseqs-uniref90",
        reference=config["uniref90"],
    threads: config["threads_28"]
    log:
        mmseqs_log=f'{config["rdir"]}/assembly-prokka-mmseqs-uniref90/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.prokka.megahit.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-prokka-mmseqs-uniref90/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.prokka.megahit.{{config}}.bmk'
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

        /home/ctools/conda/miniconda3/envs/prokka/bin/prodigal -c -i {params.ancientFasta} -a {params.orfAA} -d {params.orfNuc} >> {log.mmseqs_log} 2>&1

        {params.mmseqs_bin} easy-search \
            {params.orfAA} \
            {params.reference} \
            {output.tsv} \
            {params.tmp} \
            --max-seq-len 1000000 --search-type 1 --cov-mode 1 \
            --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov,qorfstart,qorfend,torfstart,torfend \
            --threads {threads} --split-memory-limit 200G >> {log.mmseqs_log} 2>&1

        #format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov,qframe,tframe,qorfstart,qorfend,torfstart,torfend


        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
