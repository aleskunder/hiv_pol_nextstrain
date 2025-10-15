configfile: "config.yaml"


rule all:
    input:
        "auspice/hiv_pol.json",


rule subset_sequences:
    input:
        sequences = config["sequences"]
    output:
        "results/subset.fasta"
    params:
        n = config["subset_n"]
    shell:
        """
        seqkit sample -n {params.n} {input.sequences} > {output}
        """


rule align:
    input:
        sequences = rules.subset_sequences.output,
        reference = config["reference_full"],
        annotation = config["reference_ann"]
    output:
        alignment = "results/aligned.fasta",
        tsv = "results/nextclade.tsv"
    params:
        cds = "pol"  # or a list, e.g. "pol,gag"
    shell:
        """
        nextclade3 run \
            {input.sequences} \
            --input-ref {input.reference} \
            --input-annotation {input.annotation} \
            --cds-selection {params.cds} \
            --output-fasta {output.alignment} \
            --output-tsv {output.tsv} \
            --include-reference
        """




# rule subset_alignment:
#     input:
#         aln=config["alignment"],
#     output:
#         "results/subset.fasta",
#     params:
#         n=config["subset_n"],
#     shell:
#         """
#         seqkit sample -n {params.n} {input.aln} > {output}
#         """


rule tree:
    input:
        aln=rules.align.output.alignment,
    params:
        method=config["tree_method"],
        args=config["tree_args"],
    output:
        tree="results/tree_raw.nwk",
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
        aln=rules.subset_sequences.output,
        meta=config["metadata"],
    output:
        tree="results/augur/hiv_refined.nwk",
        node_data="results/augur/hiv_node_data.json",
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
        aln=rules.subset_sequences.output,
    output:
        "results/augur/nt_muts.json",
    shell:
        """
        augur ancestral \
          --tree {input.tree} \
          --alignment {input.aln} \
          --output-node-data {output}
        """


rule translate:
    input:
        tree=rules.refine_tree.output.tree,
        node_data=rules.ancestral.output,
        reference=config["reference_ann"],
    output:
        "results/aa_muts.json",
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output} \
            --genes pol \
        """


rule export:
    input:
        tree=rules.refine_tree.output.tree,
        metadata=config["metadata"],
        branch_lengths=rules.refine_tree.output.node_data,
        nt_muts=rules.ancestral.output,
        aa_muts=rules.translate.output,
    output:
        "auspice/hiv_pol.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts}\
            --metadata-id-columns accession \
            --output {output}
        """


rule club:
    shell:
        "rm -rf results/* auspice/*"
