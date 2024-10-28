import pandas as pd
from scipy.io import mmread, mmwrite

# Read the observations.txt file and create a dictionary mapping observation IDs to tissues
RUNS = ["small_GSM3587010"]
RUNS = pd.read_csv("data.ids.txt", header=None, names=["observation"])["observation"].tolist()

rule all:
    input:
        expand("data/{run}/mx_out/assignments.txt", run=RUNS),
        "reference/nac/index.idx",
        "reference/nac/t2g.txt",
        expand("data/{run}/kb_out_nac/counts_unfiltered/cells_x_genes.cell.mtx", run=RUNS),
        # expand("data/{run}/mx_out/adata.h5ad", run=RUNS),
        expand("data/{run}/mx_out/assignments.matrix.mtx", run=RUNS),
        expand("data/{run}/mx_out/assignments.bcs.txt", run=RUNS),
        expand("data/{run}/mx_out/assignments.genes.txt", run=RUNS),
        expand("data/{run}/mx_out/assignments.json", run=RUNS),
        expand("data/{run}/mx_out/degs.txt", run=RUNS),
        expand("data/{run}/markers/data_markers.gid.txt", run=RUNS),
        expand("data/{run}/markers/data_markers.gid.json", run=RUNS),
        expand("data/{run}/markers/data_markers.txt", run=RUNS)

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

rule kb_ref:
    input:
        fasta="reference/dna.fa.gz",
        gtf="reference/dna.gtf.gz"
    output:
        index="reference/nac/index.idx",
        t2g="reference/nac/t2g.txt",
        spl_fa="reference/nac/spl.fa",
        unspl_fa="reference/nac/unspl.fa",
        spl_t2c="reference/nac/spl.t2c.txt",
        unspl_t2c="reference/nac/unspl.t2c.txt"
    shell:
        "kb ref --workflow nac -i {output.index} -g {output.t2g} -f1 {output.spl_fa} -f2 {output.unspl_fa} -c1 {output.spl_t2c} -c2 {output.unspl_t2c} {input.fasta} {input.gtf}"

rule list_fastqs:
    output:
        fastq_list="data/{run}/fastqs/fastq_list.txt"
    shell:
        "ls data/{wildcards.run}/fastqs/*.fastq.gz > {output.fastq_list}"

rule kb_count:
    input:
        index="reference/nac/index.idx",
        t2g="reference/nac/t2g.txt",
        spl_t2c="reference/nac/spl.t2c.txt",
        unspl_t2c="reference/nac/unspl.t2c.txt",
        spec=lambda wildcards: f"data/{wildcards.run}/spec.yaml"
    output:
        dir=directory("data/{run}/kb_out_nac"),
        adata_h5ad="data/{run}/kb_out_nac/counts_unfiltered/adata.h5ad",
        ambig_mtx="data/{run}/kb_out_nac/counts_unfiltered/cells_x_genes.ambiguous.mtx",
        barcodes="data/{run}/kb_out_nac/counts_unfiltered/cells_x_genes.barcodes.txt",
        cell_mtx="data/{run}/kb_out_nac/counts_unfiltered/cells_x_genes.cell.mtx",
        genes_names="data/{run}/kb_out_nac/counts_unfiltered/cells_x_genes.genes.names.txt",
        genes_txt="data/{run}/kb_out_nac/counts_unfiltered/cells_x_genes.genes.txt",
        mat_mtx="data/{run}/kb_out_nac/counts_unfiltered/cells_x_genes.mature.mtx",
        nascent_mtx="data/{run}/kb_out_nac/counts_unfiltered/cells_x_genes.nascent.mtx",
        matrix_ec="data/{run}/kb_out_nac/matrix.ec",
        inspect_json="data/{run}/kb_out_nac/inspect.json",
        kb_info="data/{run}/kb_out_nac/kb_info.json",
        output_bus="data/{run}/kb_out_nac/output.bus",
        output_unfiltered_bus="data/{run}/kb_out_nac/output.unfiltered.bus",
        run_info="data/{run}/kb_out_nac/run_info.json",
        transcripts="data/{run}/kb_out_nac/transcripts.txt"
    shell:
        """
        ol=$(seqspec onlist -o data/{wildcards.run}/onlist.txt -s region-type -i barcode -m rna {input.spec})
        x=$(seqspec index -t kb -s file -m rna {input.spec})
        f=$(seqspec file --fullpath -s read -f paired -k url -m rna {input.spec} | tr "\\t\\n" " ")
        cmd="kb count --strand unstranded --sum cell --filter bustools -t 8 -m 16G -x $x -w $ol -i {input.index} -g {input.t2g} -c1 {input.spl_t2c} -c2 {input.unspl_t2c} --h5ad --workflow=nac -o {output.dir} $f"
        echo $cmd
        bash -c "$cmd"
        """

rule mx_filter:
    input:
        umtx="data/{run}/kb_out_nac/counts_unfiltered/cells_x_genes.cell.mtx",
        ubarcodes="data/{run}/kb_out_nac/counts_unfiltered/cells_x_genes.barcodes.txt",
        mtx="data/{run}/kb_out_nac/counts_filtered/cells_x_genes.cell.mtx",
        barcodes="data/{run}/kb_out_nac/counts_filtered/cells_x_genes.barcodes.txt"
    output:
        barcodes="data/{run}/mx_out/barcodes.txt",
        mtx="data/{run}/mx_out/matrix.mtx"
    shell:
        "cp {input.barcodes} {output.barcodes} && cp {input.mtx} {output.mtx}"
        # "mx filter -bi {input.barcodes} -bo {output.barcodes} -o {output.mtx} {input.mtx}"

rule mx_normalize:
    input:
        "data/{run}/mx_out/matrix.mtx"
    output:
        "data/{run}/mx_out/norm.mtx"
    shell:
        "mx normalize -o {output} {input}"

rule ec_index:
    input:
        markers = "data/{run}/markers/markers.clean.txt",
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
        genes_input="data/{run}/kb_out_nac/counts_unfiltered/cells_x_genes.genes.names.txt"
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
        markers = "data/{run}/markers/markers.clean.txt",
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

rule mx_assign:
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

rule make_filtered_adata:
    input:
        mtx="data/{run}/mx_out/matrix.mtx",
        barcodes="data/{run}/mx_out/barcodes.txt",
        genes_txt="data/{run}/kb_out_nac/counts_unfiltered/cells_x_genes.genes.txt",
        gene_names="data/{run}/kb_out_nac/counts_unfiltered/cells_x_genes.genes.names.txt",
        assignments="data/{run}/mx_out/assignments.txt",
    output:
        # adata="data/{run}/mx_out/adata.h5ad",
        out_mtx="data/{run}/mx_out/assignments.matrix.mtx",
        out_bcs = "data/{run}/mx_out/assignments.bcs.txt",
        out_genes = "data/{run}/mx_out/assignments.genes.txt"
    run:
        from kb_python.utils import import_matrix_as_anndata
        import pandas as pd

        # Import the matrix as an AnnData object
        adata = import_matrix_as_anndata(
            input.mtx,
            input.barcodes,
            input.genes_txt
        )

        # Load gene names
        gene_names = pd.read_csv(input.gene_names, header=None, names=["gname"])
        adata.var["gname"] = gene_names["gname"].values

        # Load assignments and filter the AnnData object
        assignments = pd.read_csv(input.assignments, sep="\t", index_col=0)
        fadata = adata[assignments.index].copy()

        # Merge the assignments with the observations
        fadata.obs = pd.concat([fadata.obs, assignments], axis=1)

        # Save the filtered AnnData object
        fadata.write_h5ad(output.adata)
        mmwrite(output.out_mtx, fadata.X.tocsr())
        fadata.obs.index.to_series().to_csv(output.out_bcs, sep="\t", index=False, header=False)
        fadata.var.index.to_series().to_csv(output.out_genes, sep="\t", index=False, header=False)

rule format_assignments:
    input:
        assignments="data/{run}/mx_out/assignments.txt"
    output:
        json_file="data/{run}/mx_out/assignments.json"
    run:
        import pandas as pd
        import json
        
        # Read the assignments file
        a = pd.read_csv(input.assignments, sep="\t", index_col=0)
        
        # Create a DataFrame with counts and fractions for each label
        df = pd.DataFrame(a.groupby("label").size(), columns=["label_num"])
        df["label_frac"] = df["label_num"] / a.shape[0]
        df["label_id"] = df.index.map(a.reset_index()[["label", "label_id"]].drop_duplicates().set_index("label")["label_id"])
        df["label_ent_mean"] = df.index.map(a.groupby("label")["ent"].mean())
        df["label_ent_var"] = df.index.map(a.groupby("label")["ent"].var())

        
        # Convert the DataFrame to a list of dictionaries for JSON export
        json_list = df.reset_index().rename(columns={'index': 'celltype'}).to_dict(orient='records')
        
        # Write the JSON file
        with open(output.json_file, 'w') as f:
            json.dump(json_list, f, indent=4)

rule run_mx_norm_assignments:
    input:
        mtx = "data/{run}/mx_out/assignments.matrix.mtx"
    output:
        nmtx = "data/{run}/mx_out/assignments.norm.mtx"
    shell:
        "mx normalize -m log1pPF -o {output.nmtx} {input.mtx}"

rule run_mx_diff:
    input:
        assignments = "data/{run}/mx_out/assignments.txt",
        genes = "data/{run}/mx_out/assignments.genes.txt",
        barcodes = "data/{run}/mx_out/assignments.bcs.txt",
        norm_matrix = "data/{run}/mx_out/assignments.norm.mtx"
    output:
        degs = "data/{run}/mx_out/degs.txt"
    shell:
        """
        mx diff -a {input.assignments} -gi {input.genes} -b {input.barcodes} -o {output.degs} {input.norm_matrix}
        """

rule ec_mark:
    input:
        degs = "data/{run}/mx_out/degs.txt"
    output:
        markers_gid = "data/{run}/markers/data_markers.gid.txt",
        markers_gid_json = "data/{run}/markers/data_markers.gid.json"
    params:
        p = 0.05,
        f = 0.75,
        g = 10,
        m = 20
    shell:
        """
        ec mark -p {params.p} -f {params.f} -g {params.g} -m {params.m} -o {output.markers_gid} {input.degs}
        ec mark -p {params.p} -f {params.f} -g {params.g} -m {params.m} -o {output.markers_gid_json} -fmt json {input.degs}
        """

rule ec_convert_gid_gname:
    input:
        markers_gid = "data/{run}/markers/data_markers.gid.txt"
    output:
        bad_targets = "data/{run}/markers/data_bad_targets.txt",
        markers_gname = "data/{run}/markers/data_markers.txt"
    shell:
        """
        ec convert -m reference/gid_gname.map.txt -o {output.markers_gname} -b {output.bad_targets} {input.markers_gid}
        """
