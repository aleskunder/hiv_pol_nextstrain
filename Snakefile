configfile: "config.yaml"

RESULTS = config.get("results_dir", "results")
AUSPICE = config.get("auspice_dir", "auspice")

rule all:
    input:
        f"{AUSPICE}/hiv_pol.json",


if config.get("subset_n", 0):

    rule subset:
        input:
            sequences=config["sequences"],
        output:
            f"{RESULTS}/subset.fasta",
        params:
            n=config["subset_n"],
            seed=42, # for debugging purposes
        shell:
            """
            seqkit sample -n {params.n} {input.sequences} -s {params.seed} > {output}
            """
    seq_input = rules.subset.output

else:
    seq_input = config["sequences"]


rule download_dataset:
    output:
        dataset="/scicore/home/neher/kuznet0001/hiv_analysis/data/nextclade_hiv_ds"
    params:
        dataset_path=config["nextclade_dataset_path"],
    shell:
        """
        nextclade3 dataset get --name {params.dataset_path} --output-dir {output.dataset}
        """

rule align:
    input:
        sequences=seq_input,
        dataset=rules.download_dataset.output.dataset,
    output:
        alignment=f"{RESULTS}/aligned.fasta",
        tsv=f"{RESULTS}/nextclade.tsv",
    params:
        genes=config["cds"],
    threads: config["threads"]
    shell:
        """
        nextclade3 run -j {threads}\
            {input.sequences} \
            -D {input.dataset} \
            --output-fasta {output.alignment} \
            --output-tsv {output.tsv} \
            --min-seed-cover 0.27 \
            --include-reference
        """

# rule mask:
#     input:
#         aln=rules.align.output.alignment,
#     output:
#         masked=f"{RESULTS}/aligned_masked.fasta",
#     shell:
#         """
#         augur mask \
#           --sequences {input.aln} \
#           --mask-from-beginning 1798 \
#           --mask-from-end 4542 \
#           --output {output.masked}
#         """

rule extract_pol:
    input:
        aln=rules.align.output.alignment,
    output:
        pol=f"{RESULTS}/aligned_pol.fasta",
    shell:
        """
        seqkit subseq --region 1799:4639 {input.aln} > {output.pol}
        """



rule filter:
    input:
        aln=rules.extract_pol.output.pol,
        meta=rules.align.output.tsv,
    output:
        filtered=f"{RESULTS}/aligned_filtered.fasta",
    shell:
        """
        augur filter \
          --sequences {input.aln} \
          --metadata {input.meta} \
          --metadata-id-columns seqName \
          --exclude-where 'qc.overallStatus=bad' 'qc.overallStatus=mediocre'\
          --min-length 2500 \
          --output-sequences {output.filtered} \
        """


rule tree:
    input:
        aln=rules.filter.output.filtered,
    params:
        method=config["tree_method"],
        args=config["tree_args"],
    output:
        tree=f"{RESULTS}/tree_raw.nwk",
    threads: config["threads"]
    shell:
        """
        augur tree \
        --method {params.method} \
        --alignment {input.aln} \
        --nthreads {threads} \
        --output {output.tree}
        """


rule refine_tree:
    input:
        tree=rules.tree.output.tree,
        aln=rules.filter.output.filtered,
        meta=config["metadata"],
    output:
        tree=f"{RESULTS}/augur/hiv_refined.nwk",
        node_data=f"{RESULTS}/augur/hiv_node_data.json",
    shell:
        """
        augur refine \
          --tree {input.tree} \
          --alignment {input.aln} \
          --metadata {input.meta} \
          --root mid_point \
          --output-tree {output.tree} \
          --metadata-id-columns accession \
          --output-node-data {output.node_data}
        """

rule ancestral:
    input:
        tree=rules.refine_tree.output.tree,
        aln=rules.filter.output.filtered,
    output:
        f"{RESULTS}/augur/nt_muts.json",
    shell:
        """
        augur ancestral \
          --tree {input.tree} \
          --alignment {input.aln} \
          --output-node-data {output}
        """

rule translate:
    input:
        tree = rules.refine_tree.output.tree,
        node_data = rules.ancestral.output,
        reference = config["reference_ann"],
    output:
        f"{RESULTS}/augur/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output} \
            --alignment-output {RESULTS}/aligned_aa_%GENE.fasta
        """


rule export:
    input:
        tree=rules.refine_tree.output.tree,
        metadata=config["metadata"],
        branch_lengths=rules.refine_tree.output.node_data,
        nt_muts=rules.ancestral.output,
        aa_muts=rules.translate.output,
        auspice_config=config["auspice_config"],
    output:
        f"{AUSPICE}/hiv_pol.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} \
            --metadata-id-columns accession \
            --auspice-config {input.auspice_config} \
            --output {output}
        """


rule club:
    shell:
        f"rm -rf {RESULTS}/* {AUSPICE}/*"
