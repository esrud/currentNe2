# currentNe2
Estimation of current effective population using artificial neural networks
# System Requirements 
## Hardware requirements
`currentne2` requires only a standard computer with enough RAM to support the in-memory operations
## Software requirements
### OS Requirements
This program has been tested on the following operating sytems but should work on most linux distributions:
  - Linux: Ubuntu 20.04, Arch Linux (kernel 6.11.6), Debian Buster
### Software required for compiling
#### g++
It should be included in (almost) every linux distribution. To install it in debian-like distributions:
```
sudo apt install g++
```
#### openmp
Depending on your operating system, you may need to install openmp
```
sudo apt install libomp-dev
```
#### make
This is only needed if you want to use the make commands to compile the program. You could directly compile it using g++. To install it just run the following command (in debian-like distributions):
```
sudo apt install make
```
# Installation
## Pre-built version
It is recommended that you compile the program, as it will most likely run faster. However we provide a prebuilt version in the releases section of this repo.

## Compiling
Clone the github repo:
```
git clone https://github.com/esrud/currentNe2
```
Compile:
```
cd currentNe2
make
```
# Usage
Please note that the program makes extensive use of RAM memory. This is dependent on the amount of loci and individuals to be sampled, as well as the maximum number of chromosomes and the maximum distance between loci.
This can be increased (or reduced) by tweaking the constants *MAXLOCI*, *MAXIND*, *MAXCROMO* and *MAXDIST* respectively in the source code.
Note that increasing those values will increase the free RAM requirements.

If no output path is provided (option -o), running the program will produce a file with the same name as the input file (without the extension) with the added suffix of *_currentNe2_OUTPUT.txt* or *_currentNe2_mix_OUTPUT.txt* when running the program with "-x". This file is stored in the same path as the analyzed data.
```
currentNe2 - Current Ne estimator (v 2.0 - March 2024)
Authors: Enrique Santiago - Carlos Köpke

USAGE: ./currentne2 [OPTIONS] <filename_with_extension> [Genome_size_in_Morgans]
         where filename is the name of the data file in vcf, ped or tped
             format. The filename must include the .vcf, .ped or .tped
             extension, depending on the format.
         If Genome_size_in_Morgans (~chromosomes number) is specified,
         markers are assumed to be evenly distributed in the genome, and:
             If chromosome assignments are available in the input
               file (or in a accompanying map file when using ped format),
               two Ne estimates are calculated: one based on all the
               SNP pairs and another based only on SNP pairs located in
               different chromosomes. Chromosome sizes are assumed to be
               proportional to the number of markers.
             If chromosome assignments are not available, chromosomes are
               assumed to be one Morgan long and markers are assumed to
               be evenly distributed between chromosomes. Only one estimate
               based on all SNP pairs is made.
         If Genome_size_in_Morgans is not specified, then information on
         marker locations within chromosomes is obligatory, and:
             If a constant recombination rate in cM/Mb is specified with
               the option r, physical locations of markers in bp are
               converted to genetic distances in cM and the Ne estimate
               using all SNP pairs considers these distances.
             If the option r is not used, then the information about marker
               locations in a genetic map (i.e., in cM) is used (only
               applicable for ped/map and tped formats).

OPTIONS:
   -h    Print this manual.
   -s    Number of SNPs to use in the analysis (all by default).
   -k    -If a POSITIVE NUMBER is given, the number of full siblings that
         a random individual has IN THE POPULATION (the population is the
         set of reproducers). With full lifetime monogamy k=2, with 50%
         of monogamy k=1 and so on. With one litter per multiparous
         female k=2, with two litters per female sired by the same father
         k=2, but if sired by different fathers k=1, in general, k=2/Le
         where Le is the effective number of litters (Santiago et al. 2023).
         -If ZERO is specified (i.e., -k 0), each offspring is assumed to
         be from a new random pairing.
         -If a NEGATIVE NUMBER is specified, the average number of full
         siblings observed per individual IN THE SAMPLE. The number k of
         full siblings in the population is estimated together with Ne.
         -BY DEFAULT, i.e. if the modifier is not used, the average number
         of full siblings k is estimated from the input data.
   -r    Recombination rate to convert the physical locations in the
         input file to genetic locations in cM.
   -x    The sample consists of a random set of individuals from
         a metapopulation with two subpopulations of equal size.
         The k parameter cannot be used when this option is used.
         This option is only applicable when the assigments of markers
         to chromosomes are available. A single estimate of Ne is generated.
   -o    Specifies the output file name. If not specified, the output
         filename is built from the name of the input file.
   -t    Number of threads (default: 8)

EXAMPLES:
   - Random mating. Full sibs will be identified from the input data.
     Genetic distances between markers are obtained from the physical
     locations in the input file using a constant rate of 1cM/Mb:
         ./currentne2 -r 1.0 filename
   - Random mating and 20 Morgans (equivalent to a genome of 20 chromosomes),
     assuming that full siblings are no more common than expected under
     random pairing (each offspring from a new random pairing):
         ./currentne2 -k 0 filename 20
   - A subsample of 10000 SNPs will be analysed. Individuals are assumed to
     come from a metapopulation composed of two subpopulations. The length of
     the genome is 30 Morgans and the marker assigments to chromosomes are
     available in the input file:
         ./currentne2 -x -s 100000 filename 30
   - Full siblings will be identified from the genotyping data. The genetic
     locations (in cM) are available in the input file (only for ped/map and
     tped formats):
         ./currentne2 filename
   - Two full siblings per individual (k = 2) IN THE POPULATION:
         ./currentne2 -k 2 filename 20
   - An 80% of lifetime monogamy in the population. Output filename specified:
         ./currentne2 -k 1.6 -o output_filename filename 20
     (with a monogamy rate m = 0.80, the expected number of full
     siblings that a random individual has is k = 2*m = 1.6)
   - If 0.2 full siblings per individual are OBSERVED IN THE SAMPLE:
         ./currentne2 -k -0.2 filename 20
     (NOTE the MINUS SIGN before the number of full sibling 0.2)
```
# Demo
In this repo we provide some sample data for testing purposes, kindly provided by John Taggart (Institute of Aquaculture, University of Stirling). Here are two examples on how to test `currentne2` on this data:
- Random mating. Genetic distances between markers are obtained from the physical locations in the input file, using a constant rate of 1cM/Mb:

```
# Decompress the sample data
tar -xvf example/sample_data.tar.gz -C ./example

# Run currentne2 on the sample data. This should take 5-10 min
./currentne2 -r 1 example/baddoch_clean.ped

# Once it finishes running you can check the resulting file with your preferred text editor
less example/baddoch_clean_currentNe2_OUTPUT.txt
```
- Random mating. Individuals are assumed to come from a metapopulation composed of two subpopulations. Genetic distances between markers are obtained from the physical locations in the input file, using a constant rate of 1cM/Mb:

```
# Decompress the sample data (if it hasn't been decompressed previously)
tar -xvf example/sample_data.tar.gz -C ./example

# Run currentne2 on the sample data. This should take 5-10 min
./currentne2 -r 1 -x example/girnock_clean.ped

less example/girnock_clean_currentNe2_mix_OUTPUT.txt
```
# Acknowledgements
This study forms part of the Marine Science Programme (ThinkInAzul) supported by the Ministerio de Ciencia e Innovación and Xunta de Galicia with funding from the European Union NextGenerationEU (PRTR-C17.I1) and European Maritime and Fisheries Fund.
# How to cite
Santiago, E., Köpke, C. & Caballero, A. Accounting for population structure and data quality in demographic inference with linkage disequilibrium methods. Nat Commun 16, 6054 (2025). [https://doi.org/10.1038/s41467-025-61378-w](https://doi.org/10.1038/s41467-025-61378-w)
