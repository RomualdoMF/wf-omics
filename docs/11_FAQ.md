If your question is not answered here, please start a discussion on the [community](https://community.nanoporetech.com/). 

### How do I use Rerio models for SNV calling?

If you wish to use cutting-edge research release models from Rerio with Clair3, follow the instructions on the [Clair3 section of the Rerio repository](https://github.com/nanoporetech/rerio#clair3-models). The download script will download and extract the model to a directory which can then be provided to the workflow's `--clair3_model_path` option. Additionally, if the basecaller configuration you are using is not known to the workflow, you must set `--override_basecaller_cfg custom`.


### The number of SNVs and indels in the report do not sum up to the number of records, is that normal?

Yes; this can be due to some multiallelic sites carrying a mixture of SNV and indel alleles.