# tskibd-filter

## Overview

When working with simulated data to infer effective population size (via IBDNe),
it's crucial to accurately filter out IBD (Identity-by-descent) segments from
close relatives. This is especially important when these segments originate from
ancestors with a Time to Most Recent Common Ancestor (TMRCA) less than 1.5
generations. Filtering out such segments is necessary because they can introduce
oscillations in the inferred Ne (Effective Population Size) trajectory,
particularly in the most recent 1-20 generations.

## Background

The IBDNe paper outlines a method for addressing this issue:

> When analyzing true IBD segments, we removed IBD segments with TMRCA less than
> 1.5 generations ago to align with the removal of half-siblings and closer
> relatives in real data. When analyzing inferred IBD segments, we removed
> segments that overlapped with a true IBD segment with TMRCA less than 1.5
> generations ago.

## Features

This program implements the IBD-filtering process described above. It takes both
inferred IBD and true IBD data from the **same** simulation as input. It is assumed that
the input IBD data follows the `tskibd` file format, and the TMRCA for true
IBD is assumed to be reliable. The output IBD file is the filterred version of the
inferred IBD. 

## Compililation

Run the following command:
```
cargo build --release
```

## Usage

You can use the tskibd-filter as follows:

```
$ target/release/tskibd-filter --help

Usage: tskibd-filter --out <OUT> <INFER_SET> <TRUE_SET>

Arguments:
  <INFER_SET>       Path to the inferred IBD file (in tskibd format)
  <TRUE_SET>        Path to the true IBD file (in tskibd format)

Options:
  -o, --out <OUT>   Path to the output IBD file
  -h, --help        Print help information
  -V, --version     Print version information

```


