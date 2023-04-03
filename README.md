# `r-pfbwt`

[![CMake](https://github.com/marco-oliva/r-pfbwt/actions/workflows/cmake.yml/badge.svg)](https://github.com/marco-oliva/r-pfbwt/actions/workflows/cmake.yml)

`r-pfbwt` is a tool to build the run-length encoded *BWT* and the *SA* values at the run heads from the prefix-free 
parsing of the input data.

## Docker
`r-pfbwt` is available on docker:

```bash
docker pull moliva3/r-pfbwt:latest
docker run moliva3/r-pfbwt:latest rpfbwt --help
```

If using singularity:
```bash
singularity pull rpfbwt_sif docker://moliva3/r-pfbwt:latest
./rpfbwt_sif rpfbwt --help
```

## Build
`r-pfbwt` can be built using: 

```shell
git clone https://github.com/marco-oliva/r-pfbwt.git
cd r-pfbwt
mkdir build && cd build
cmake ..
make 
```
 
## How to use `r-pfbwt`
`r-pfbwt` takes as input the prefix-free parse of the input data, namely a dictionary *D<sub>1</sub>* and a parse *P<sub>1</sub>*, and 
the dictionary *D<sub>2</sub>* and the parse *P<sub>2</sub>* obtained by prefix-free parsing *P<sub>1</sub>*. Note that `r-pfbwt` does not use *P<sub>1</sub>*,
it only uses *D<sub>1</sub>*, *D<sub>2</sub>* and *P<sub>2</sub>*. In order to compute the prefix-free parsing of the input data we can use `pfp++`
([link here](https://github.com/marco-oliva/pfp.git)). The following example computes the run-length encoded *BWT* of multiple
sequences of yeast.

```shell
wget https://gitlab.com/manzai/Big-BWT/-/blob/f67022fe74dae0234e516324103613a0fdd58a6e/yeast.fasta -O ./yeast.fasta
pfp++ -f yeast.fasta -w 10 -p 100 --output-occurrences 
pfp++ -i yeast.fasta.parse -w 5 -p 11 
rpfbwt  --l1-prefix yeast.fasta --w1 10 --w2 5 --threads 10
```


## Usage
We report here all the available parameters for `r-pfbwt`

```shell
rpfbwt
Usage: rpfbwt [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  --l1-prefix TEXT REQUIRED   Level 1 Prefix.
  --w1 UINT REQUIRED          Level 1 window length.
  --w2 UINT REQUIRED          Level 2 window length.
  -t,--threads UINT           Number of threads.
  --chunks UINT:INT in [1 - 1000]
                              Number of chunks.
  --integer-shift UINT        Integer shift used during parsing.
  --tmp-dir TEXT:DIR          Temporary files directory.
  --bwt-only                  Only compute the RLBWT. No SA values.
  --version                   Version number.
  --configure                 Read an ini file
```

## Citation

Please, if you use this tool in an academic setting cite the following papers:

```
@article{oliva2023recursiveBWT,
  title={Recursive Prefix-Free Parsing for Building Big BWTs},
  author={Oliva, Marco and Gagie, Travis and Boucher, Christina},
  journal={bioRxiv},
  pages={2023--01},
  year={2023},
  publisher={Cold Spring Harbor Laboratory},
  doi = {10.1101/2023.01.18.524557}
}
```

```
@article {oliva2023recursiveRI,
    title = {Building a Pangenome Alignment Index via Recursive Prefix-Free Parsing},
	author={Oliva, Marco and Gagie, Travis and Boucher, Christina},
	journal={bioRxiv},
    pages={2023--01},
    year={2023},
    publisher={Cold Spring Harbor Laboratory}
    doi = {10.1101/2023.01.26.525723},
}
```
