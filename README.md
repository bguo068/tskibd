# `tskibd`

`tskibd` is a tool designed to calculate Identity by Descent (IBD) segments from genealogical trees represented in the [tree sequence](https://tskit.dev/tutorials/what_is.html#sec-what-is) format. 
It leverages the tskit C API library, which you can find at [tskit](https://github.com/tskit-dev/tskit) repository.


## System requirements

This software has been tested on MacOS and Linux operating systems. The
dependencies (including version numbers) for compilation and execution are
defined in the env.yml Conda recipe. Installation time is typically within a
minute.

## Conda environment for compilation

If you don't have `Conda` installed on your system, you can easily install it from [here](https://docs.conda.io/en/latest/miniconda.html).

Once `Conda` is installed, create the `tskit` environment using the recipe
provided in `env.yml` by running:
```sh
conda env create -f ./env.yml
conda activate tskibd
```

## Compilation

The `meson/ninja` build system is used to compile `tskibd`. 
Please refer to the previous section to prepare the compilation envrionment.

```sh
git clone https://github.com/bguo068/tskibd.git
cd tskibd
git submodule update --init --recursive # IMPORTANT!
mv tskit/c/VERSION tskit/c/VERSION.txt # Fix filename issue
meson setup build
ninja -C build tskibd
```
The compiled executable can be found at `build/tskibd`.



## Usage

```
Usage: tskibd <chromN> <<bp_per_cm/recombination_rate_map_path>> <sampling_window> <mincm> <treeseq_file> [<out>]

Positional parameters:

    <chr_no>          Chromosome number (integer). This parameter is used for
                      naming the output files, such as 1.ibd and 1.map.

    <bp_per_cm/recombination_rate_map_path>
                      This positional argument can be used in two ways:
                      (1) if it is a numeric type, it will be interpreted as the 
                          base pairs per centimorgan.
                          For example, use '1000000' for human or 
                         '15000' for p.f. (P. falciparum).
                      (2) If it is a string, it will be interpreted as the path
                      to a recombination rate map file. This file is a
                      space/tab-delimited two-column table without a header. The
                      first column contains integers indicating the coordinates
                      in base-pair, and the second column contains float numbers
                      indicating the coordinates in centimorgan.  For positions
                      that is not explicited included, we use linear 
                      interpolation/extrapolation for bp to cm mapping.

    <sampling_window> Sampling window size in base pairs. Use '1' to check all
                      trees or '1000' to check trees covering positions of 0,
                      1000, 2000 bp, and so on. We recommend a window size
                      equivalent to 0.01 centimorgan. For instance, for human
                      use 10000 bp, and for p.f. use 150.

    <mincm>           Minimum length in centimorgan of IBD (Identity by Descent)
                      segment to keep. IBD segments shorter than this value
                      will be ignored and not written to the output file.

    <treeseq_file>    Tree sequence file that has finished coalescent.

    <out>             Optional: output IBD filename.
                      If unspecified, it will be '{chrno}.ibd'
                      Otherwise, it will be '{out}.ibd'
```

Example running on the **test data**: 
```
# use default output prefix
build/tskibd 1 15000 150 2  example_data/chr1.trees
# or use specified output prefix
build/tskibd 1 15000 150 2  example_data/chr1.trees out/chr1

# from version 0.0.2, we allow specifying recombination rate map
# instead of a constant recombination rate
echo "750000 50.0" > rate_map.txt
echo "1500000 100.0" >> rate_map.txt
build/tskibd 1 rate_map.txt 150 2  example_data/chr1.trees out/chr1
```

## Input `*.trees` file:

The input genealogical trees are expected to be in the tree sequence format.
Please ensure that the input trees are fully coalesced. For trees simulated from
[`msprime`](https://github.com/tskit-dev/msprime), they are typically fully
coalesced. However, if the trees are from a forward simulator like
[`SLiM`](https://github.com/MesserLab/SLiM), additional coalescent simulation
might be required to ensure full coalescence (refer to [`PySLiM`
tutorial](https://tskit.dev/pyslim/docs/latest/tutorial.html) for more details).

Notes on mutational information of a tree sequence: the mutational information
is not used to infer IBD segments but is used to annotate IBD segments 
(the last column of the output `.ibd` file). 
If all mutations are included in the tree sequeunce, `tskibd` might run very slow.
We recommand directly use tree sequence before mutations are added (no mutations)
or keep only those for sites under selection (selected mutation).  

## Output `*.ibd` file:
- col1 - Id1: node id 1 (a sample node of the specified tree sequence).
- col2 - Id2: node id 2.
- col3 - Start: IBD segment start position in base pair.
- col4 - End: IBD segment end position in base pair.
- col5 - Ancestor: the node id of the ancestor from which the IBD segment is originated.
- col6 - Tmrca: the time (backward) of the ancestor node.
- col7 - HasMutation: whether the IBD segment carries a mutation. (1: yes, 0: no)


# Citation:

If you find this tool useful, please cite our preprint https://doi.org/10.1101/2023.07.14.549114
