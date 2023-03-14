## Configuration for ribosomal protein tree pipeline

`config.yaml` - main configuration for the pipeline. Here you supply a path to the csv file with genome_ids and genome_paths, as well as set parameters for the alignment and tree

```
genome_paths: "config/genome_paths.csv"
alignment_params:
  filter_proportion_missing: 0.1
tree_params:
  outgroup: "Escherichia coli"
  model: "GTR+G"
  num_parsimony_trees: 10
  num_bootstraps: 10
```

`genome_paths.csv` - a csv table mapping genome_ids to genome_paths. other columns can be added and are ignored by the pipeline.

```
genome_id,genome_path,population
Vibrio_tasmaniensis_1,test_data/10N.222.45.E7_contigs.fa.gz,Vibrio tasmaniensis
Vibrio_tasmaniensis_2,test_data/10N.222.45.A2_contigs.fa.gz,Vibrio tasmaniensis
Vibrio_tasmaniensis_3,test_data/10N.222.48.A2_contigs.fa.gz,Vibrio tasmaniensis
Vibrio_tasmaniensis_4,test_data/10N.222.48.B1_contigs.fa.gz,Vibrio tasmaniensis
...
```