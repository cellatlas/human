import pandas as pd

# Read the observations.txt file and create a dictionary mapping observation IDs to tissues
observations_df = pd.read_csv("observations.txt", sep="\t", header=None, names=["observation", "tissue", "technology", "status"])
observations_df = observations_df.query("status == 'ok'")
observations_dict = dict(zip(observations_df["observation"], observations_df["tissue"]))

# Use the observation IDs from the first column of observations_df as the RUNS list
RUNS = observations_df["observation"].tolist()
RUNS = ["small_GSM3587010"]

rule all:
    input:
        expand("data/{run}/mx_out/assignments.txt", run=RUNS),
        "reference/index.idx",
        "reference/t2g.txt",
        "reference/transcriptome.fa",
	expand("data/{run}/kb_out/counts_unfiltered/cells_x_genes.mtx", run=RUNS)

rule download_references:
    output:
        json="reference/reference.json",
        fasta="reference/dna.fa.gz",
        gtf="reference/dna.gtf.gz"
    shell:
        """
        gget ref -o {output.json} -w dna,gtf homo_sapiens
        url=$(jq -r '.homo_sapiens.genome_dna.ftp' {output.json}) && curl -o {output.fasta} $url
        url=$(jq -r '.homo_sapiens.annotation_gtf.ftp' {output.json}) && curl -o {output.gtf} $url
        """

rule create_kb_reference:
    input:
        fasta="reference/dna.fa.gz",
        gtf="reference/dna.gtf.gz"
    output:
        index="reference/index.idx",
        t2g="reference/t2g.txt",
        transcriptome="reference/transcriptome.fa"
    shell:
        "kb ref -i {output.index} -g {output.t2g} -f1 {output.transcriptome} {input.fasta} {input.gtf}"

rule list_fastqs:
    output:
        fastq_list="data/{run}/fastqs/fastq_list.txt"
    shell:
        "ls data/{wildcards.run}/fastqs/*.fastq.gz > {output.fastq_list}"

rule kb_count:
    input:
        index="reference/index.idx",
        t2g="reference/t2g.txt",
        spec=lambda wildcards: f"data/{wildcards.run}/fastqs/spec.yaml",
        fastq_list=rules.list_fastqs.output.fastq_list
    output:
        onlist="data/{run}/kb_out/onlist.txt",
        inspect="data/{run}/kb_out/inspect.json",
        info="data/{run}/kb_out/kb_info.json",
        matrix="data/{run}/kb_out/matrix.ec",
        bus="data/{run}/kb_out/output.bus",
        unfiltered_bus="data/{run}/kb_out/output.unfiltered.bus",
        run_info="data/{run}/kb_out/run_info.json",
        transcripts="data/{run}/kb_out/transcripts.txt",
        # unfiltered_adata="data/{run}/kb_out/counts_unfiltered/adata.h5ad",
        unfiltered_barcodes="data/{run}/kb_out/counts_unfiltered/cells_x_genes.barcodes.txt",
        unfiltered_genes="data/{run}/kb_out/counts_unfiltered/cells_x_genes.genes.txt",
        unfiltered_gene_names="data/{run}/kb_out/counts_unfiltered/cells_x_genes.genes.names.txt",
        unfiltered_mtx="data/{run}/kb_out/counts_unfiltered/cells_x_genes.mtx"
    params:
        fastqs_csv=lambda wildcards, input: ",".join(open(input.fastq_list).read().strip().split("\n")),
        fastqs_space=lambda wildcards, input: " ".join(open(input.fastq_list).read().strip().split("\n")),
        outdir=lambda wildcards: f"data/{wildcards.run}/kb_out"
    shell:
        """
        kb count -t 64 -m 128G -i {input.index} -g {input.t2g} -o {params.outdir} \
            -w $(seqspec onlist -m rna -o {output.onlist} -s region-type -r barcode {input.spec}) \
            -x $(seqspec index -t kb -m rna -r {params.fastqs_csv} {input.spec}) \
            {params.fastqs_space}
        """

rule mx_filter:
    input:
        mtx="data/{run}/kb_out/counts_unfiltered/cells_x_genes.mtx",
        barcodes="data/{run}/kb_out/counts_unfiltered/cells_x_genes.barcodes.txt"
    output:
        barcodes="data/{run}/mx_out/barcodes.txt",
        mtx="data/{run}/mx_out/matrix.mtx"
    shell:
        "mx filter -c 2 2 -bi {input.barcodes} -bo {output.barcodes} -o {output.mtx} {input.mtx}"

rule mx_normalize:
    input:
        "data/{run}/mx_out/matrix.mtx"
    output:
        "data/{run}/mx_out/norm.mtx"
    shell:
        "mx normalize -o {output} {input}"

rule ec_index:
    input:
        markers=lambda wildcards: f"markers/{observations_dict[wildcards.run]}/markers.txt"
    output:
        groups="data/{run}/mx_out/groups.txt",
        targets="data/{run}/mx_out/targets.txt",
        ec="data/{run}/mx_out/markers.ec.txt"
    shell:
        "ec index -g {output.groups} -t {output.targets} -e {output.ec} {input.markers}"

rule ec_index_clean:
    input:
        markers="data/{run}/mx_out/clean.markers.txt"
    output:
        groups="data/{run}/mx_out/clean.groups.txt",
        targets="data/{run}/mx_out/clean.targets.txt",
        ec="data/{run}/mx_out/clean.markers.ec.txt"
    shell:
        "ec index -g {output.groups} -t {output.targets} -e {output.ec} {input.markers}"

rule mx_extract:
    input:
        mtx="data/{run}/mx_out/norm.mtx",
        targets="data/{run}/mx_out/targets.txt",
        genes_input="data/{run}/kb_out/counts_unfiltered/cells_x_genes.genes.names.txt"
    output:
        genes="data/{run}/mx_out/extract_genes.txt",
        mtx_extracted="data/{run}/mx_out/extract.mtx"
    shell:
        "mx extract -t {input.targets} -gi {input.genes_input} -go {output.genes} -o {output.mtx_extracted} {input.mtx}"

rule mx_extract_clean:
    input:
        mtx="data/{run}/mx_out/clean.mtx",
        genes="data/{run}/mx_out/clean_genes.txt",
        targets="data/{run}/mx_out/clean.targets.txt"
    output:
        genes="data/{run}/mx_out/extract.clean_genes.txt",
        mtx_extracted="data/{run}/mx_out/extract.clean.mtx"
    shell:
        "mx extract -t {input.targets} -gi {input.genes} -go {output.genes} -o {output.mtx_extracted} {input.mtx}"

rule mx_clean:
    input:
        genes="data/{run}/mx_out/extract_genes.txt",
        barcodes="data/{run}/mx_out/barcodes.txt",
        mtx="data/{run}/mx_out/extract.mtx"
    output:
        clean_genes="data/{run}/mx_out/clean_genes.txt",
        clean_mtx="data/{run}/mx_out/clean.mtx",
        clean_barcodes="data/{run}/mx_out/clean_barcodes.txt",
        bad_targets="data/{run}/mx_out/clean_genes.txt.bad"
    shell:
        "mx clean -gi {input.genes} -go {output.clean_genes} -bi {input.barcodes} -bo {output.clean_barcodes} -o {output.clean_mtx} --bad {input.mtx}"

rule ec_filter:
    input:
        markers=lambda wildcards: f"markers/{observations_dict[wildcards.run]}/markers.txt",
        bad_targets="data/{run}/mx_out/clean_genes.txt.bad"
    output:
        clean_markers="data/{run}/mx_out/clean.markers.txt"
    shell:
        "ec filter -bt {input.bad_targets} -o {output.clean_markers} {input.markers}"

rule mx_normalize_rank:
    input:
        "data/{run}/mx_out/extract.clean.mtx"
    output:
        "data/{run}/mx_out/rank.mtx"
    shell:
        "mx normalize -m rank -o {output} {input}"

rule assign:
    input:
        mtx="data/{run}/mx_out/rank.mtx",
        groups="data/{run}/mx_out/clean.groups.txt",
        genes="data/{run}/mx_out/extract.clean_genes.txt",
        barcodes="data/{run}/mx_out/clean_barcodes.txt",
        ec="data/{run}/mx_out/clean.markers.ec.txt"
    output:
        "data/{run}/mx_out/assignments.txt"
    shell:
        "mx assign -g {input.groups} -gi {input.genes} -bi {input.barcodes} -e {input.ec} -o {output} {input.mtx}"
