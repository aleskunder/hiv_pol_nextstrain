configfile: "config.yaml"


rule all:
    input:
        "auspice/hiv_pol.json",


rule subset:
    input:
        sequences=config["sequences"],
    output:
        "results/subset.fasta",
    params:
        n=config["subset_n"],
    shell:
        """
        seqkit sample -n {params.n} {input.sequences} > {output}
        """


rule align:
    input:
        sequences=rules.subset.output,
        reference=config["reference"],
        annotation=config["reference_ann"],
    output:
        alignment="results/aligned.fasta",
        tsv="results/nextclade.tsv",
        translation=expand("results/translations/gene.{gene}.fasta", gene="pol"),
    params:
        template_string=lambda w: "results/translations/gene.{gene}.fasta",
    shell:
        """
        nextclade3 run \
            {input.sequences} \
            --input-ref {input.reference} \
            --input-annotation {input.annotation} \
            --output-fasta {output.alignment} \
            --output-tsv {output.tsv} \
            --output-translations 'results/translations/gene.{{cds}}.fasta' \
            --min-seed-cover 0.2 \
            --include-reference
        """


# rule subset:
#     input:
#         aln=config["alignment"]
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
        aln=rules.align.output.alignment,
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
        aln=rules.align.output.alignment,
        translation=expand("results/translations/gene.{gene}.fasta", gene="pol"),
    params:
        annotation=config["reference_ann"],
        translation="results/translations/gene.%GENE.fasta",
    output:
        "results/augur/muts.json",
    shell:
        """
        augur ancestral \
          --tree {input.tree} \
          --alignment {input.aln} \
          --annotation {params.annotation} \
          --translations {params.translation} \
          --genes pol \
          --output-node-data {output}
        """


# rule translate:
#     input:
#         tree=rules.refine_tree.output.tree,
#         node_data=rules.ancestral.output,
#         reference=config["reference_ann"],
#     output:
#         "results/aa_muts.json",
#     shell:
#         """
#         augur translate \
#             --tree {input.tree} \
#             --ancestral-sequences {input.node_data} \
#             --reference-sequence {input.reference} \
#             --output-node-data {output} \
#         """


rule export:
    input:
        tree=rules.refine_tree.output.tree,
        metadata=config["metadata"],
        branch_lengths=rules.refine_tree.output.node_data,
        ancestral=rules.ancestral.output,
    output:
        "auspice/hiv_pol.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.ancestral}\
            --metadata-id-columns accession \
            --output {output}
        """


rule club:
    shell:
        "rm -rf results/* auspice/*"
