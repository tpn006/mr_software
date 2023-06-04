# MR Software (CSE185 Spring 2023 Final Project)

CSE 185 Spring 2023 Final Project: MOTIF Recognition Software (MR Software)

This is the final project for group 18 for CSE 185 Spring 2023. The goal of the tool is to find which motifs are enriched in peaks. The tool takes a list of peak locations, a reference genome, and the motifs to check for enrichment. This tool is meant to replicate HOMERs known motif recognition but on a smaller scale using supplied known motifs. 

# Install instructions

As of now, our tool cannot be installed and must be ran as a python program. -This is a work in progress and will be able to be installed at a later date.-

Running requires the `pyfaidx` library to be installed. You can install `pyfaidx` with:

```
pip install pyfaidx 
```

Another option is to install `mrsoftware`. This is done by using `setuptools`. You can verify if `setuptools` is installed correctly by using `pip show setuptools`. Otherwise, you can install `setuptools` through `sudo apt-get install python3-setuptools` for python3 or `sudo apt-get install setuptools` for python2.  This allows you to use the command `python3 setup.py install` followed by `chmod +x mrsoftware` allowing the user to use `mrsoftware` through commandline. You can check for installation by running `mrsoftware --help`

# Basic usage

The basic usage of `mrsoftware` is:
```shell
python3 mrsoftware.py -m <Path to the .meme file> -g <Path to the .fa file> <Path to the .bed file>
```
or if mrsoftware is installed, you can run
```shell
mrsoftware -m <path to meme file> -g <path to genome> <path to bed>
``` 

To run `mrsoftware` on a small test example (using files in this repo):
* copy GRCm38.fa to /MR_SOFTWARE/tests/ (this file is very large ~2.6GB so not stored on github. We recommend: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/)
* go to the directory /MR_SOFTWARE/mrsoftware/
* run: 
```shell
python3 mrsoftware.py -m ../tests/short.meme -g ../tests/GRCm38.fa ../tests/short.bed -o outfile.txt
```

In this example, we choose to specify an output file using `-o outfile.txt` to store results in the file called `outfile.txt`.

The example command will output:
```
MOTIF: PO5F1_MOUSE.H11MO.0.A threshold score: 10.653805961203929 pval:1.1429535716149692e-27
MOTIF: ATOH1_MOUSE.H11MO.0.B threshold score: 7.361721218183671 pval:3.0205161724051136e-05
MOTIF: ATF2_MOUSE.H11MO.0.A threshold score: 7.143793953622376 pval:0.00309861101765903
MOTIF: AP2B_MOUSE.H11MO.0.B threshold score: 10.819923859254406 pval:0.013674086674685577
MOTIF: ASCL1_MOUSE.H11MO.0.A threshold score: 6.947104693031475 pval:0.028427706507898966
MOTIF: BACH1_MOUSE.H11MO.0.C threshold score: -18.51132984334289 pval:0.028427706507898966
MOTIF: AP2A_MOUSE.H11MO.0.A threshold score: 9.472461728667259 pval:0.05870136798384332
MOTIF: AP2C_MOUSE.H11MO.0.A threshold score: 13.188068968309581 pval:0.12041306253096071
MOTIF: ASCL2_MOUSE.H11MO.0.C threshold score: 10.114975640536887 pval:0.12041306253096071
MOTIF: ATF4_MOUSE.H11MO.0.A threshold score: 9.051958606423801 pval:0.12041306253096071
MOTIF: BACH2_MOUSE.H11MO.0.A threshold score: 11.625465502907595 pval:0.12041306253096071
MOTIF: BARX1_MOUSE.H11MO.0.C threshold score: 9.80302688149681 pval:0.12041306253096071
MOTIF: ARI5B_MOUSE.H11MO.0.C threshold score: 10.116102569642104 pval:0.245398773006135
MOTIF: ATF1_MOUSE.H11MO.0.B threshold score: 6.920941679586882 pval:0.4969325153374233
MOTIF: ATF3_MOUSE.H11MO.0.A threshold score: 10.78039215314772 pval:0.4969325153374233
MOTIF: AHR_MOUSE.H11MO.0.B threshold score: 12.499396890474504 pval:1.0
MOTIF: AIRE_MOUSE.H11MO.0.C threshold score: 8.743721746541322 pval:1.0
MOTIF: ALX1_MOUSE.H11MO.0.B threshold score: 9.668304152254336 pval:1.0
MOTIF: ANDR_MOUSE.H11MO.0.A threshold score: 8.928136456214904 pval:1.0
MOTIF: ARNT_MOUSE.H11MO.0.B threshold score: 10.501204725085982 pval:1.0
```
# Comparison to HOMER's findMotifsGenome.pl
Our tool recognizes motifs from peaks and determines if there is enrichment. We were inspired to make our tool by [HOMER](http://homer.ucsd.edu/homer/) which was also created at UCSD. Although the tools have very different levels of capabilities, HOMER has a `findMotifsGenome.pl` that we tried to replicate in our `mrsoftware`.
To compare to the output of `findMotifsGenome.pl`, run:
* `findMotifsGenome.pl <peaks.bed> <reference genome.fa> <output directory> -mask -size 100`
 * Note: `findMotifsGenome.pl` takes over 20 minutes to run so we recommend using `nohup`. 
What follows is an example of how to run `findMotifsGenome.pl` from the `/mr_software/` directory:
```
mkdir HOMER

findMotifsGenome.pl /tests/short.bed /tests/GRCm38.fa HOMER 
```
The outputs between HOMER's `findMotifsGenome.pl`
# mrsoftware options

`mrsoftware.py` requires 3 inputs: 
* `-g` followed by a reference genome in .fa format. 
* `-m` followed by a file containing motif information in a .meme format.
* and lastly a file containing peak information in a .bed format.

There is one optional input as of now:
* `-o` an output file to write output to. If the option is not used, the tool uses the default which is sysout
The plan is to have the output be a required directory to create the .png of the motif logos and the .html to easily view all results.

# File format

The output file format is currently a .txt but the plan is to create a directory with .html and .png that allows one to view the outputs easily.

# Contributors

This repository was generated by Tasha Nguyen and Jonathan Narita, following the example project for CSE 185 Spring 2023 made by Melissa Gymerrek. 