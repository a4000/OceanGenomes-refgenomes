# Computational-Biology-OceanOmics/OceanGenomes-refgenomes: Parameters

## Introduction

All the parameters in the pipeline can be set in a config file, or they can be set at the command line.

### General parameters

- input: Input sample sheet .csv file (columns: sample, hifi_dir, hic_dir, version, tolid, taxid)
- outdir: Directory where output files will be published
- binddir: Directory to bind to the containers (e.g., /scratch)
- kvalue: kvalue to use with Meryl (default = 31).
- buscomode: Busco mode. Can be `genome`, `ogs`, or `trans` (default = `genome`).
- buscodb: Busco database.
- gxdb: FCS-GX database.
- rclonedest: Destination for rclone.
