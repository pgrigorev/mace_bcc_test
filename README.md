This repository contains the data and the scripts to produce plots for section A.13 of [A foundation model for atomistic materials chemistry](https://arxiv.org/abs/2401.00096) preprint. The naming convention for pretrained models follows the [mace-foundations](https://github.com/ACEsuit/mace-foundations) repository. Follow the link to download the model files as well.

The notebooks contain the codes to produce the plots for two versions of the preprint. The folders are structured as follows:
- `data` folder contains the reference dft data with finetuning data and finetuned model in `A13_ft`  subfolder.
- `figures`: pdf versions of the figures produced by the notebooks
- `logs`: calculations logs to obtain the files in the results
- `results`: files with calculations results
- `scripts`: helper scripts to produce the  (are imported in the notebooks)

We are currently working on a public version of the workflow used to run the tests. If you would like to have an access to the workflow or test tour model please feel free to [contact Petr Grigorev](mailto:petr.y.grigorev@gmail.com).