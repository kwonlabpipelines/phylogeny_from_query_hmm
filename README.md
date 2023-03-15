## Build a tree based on an alignment of query genes found in a set of genomes

*credit to [Phil Arevalo](https://github.com/philarevalo) for the initial pipeline, initially developed in the [Polz lab](https://vds-mes.univie.ac.at/people/martin-f-polz)*

![Example Phylogenetic Tree](readme_files/example_tree.png?raw=true "Phylogenetic Tree")

To use this pipeline, have snakemake, snakedeploy installed, then run:

snakedeploy deploy-workflow https://github.com/kwonlabpipelines/phylogeny_from_query_hmm <target_folder> --tag v0.1