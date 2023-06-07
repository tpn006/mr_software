# MR Software (CSE185 Spring 2023 Final Project)

CSE 185 Spring 2023 Final Project: MOTIF Recognition Software (MR Software)

This is the final project for group 18 for CSE 185 Spring 2023. The goal of the tool is to find which motifs are enriched in peaks. The tool takes a list of peak locations, a reference genome, and the motifs to check for enrichment. This tool is meant to replicate HOMERs known motif recognition but on a smaller scale using supplied known motifs. 

# Install instructions

Our tool can be ran as a unix command line tool.

### Required Dependencies
Either way, there are a two required dependencies. Running our tool requires the `pyfaidx` and `scipy` libraries to be installed. You can install those libraries using `pip`:

```shell 
pip install scipy
pip install pyfaidx 
```

### Installing tool
Another option is to install `mrsoftware`. This is done by using `setuptools`. 
You can verify if `setuptools` is installed correctly by using `pip show setuptools`. If `setuptools` is not installed, you can install `setuptools` through `sudo apt-get install python3-setuptools` for python3. You can check if the package is installed properly by running `dpkg-query -l python3-setuptools`; if the package is installed properly, you will see information about the package installation, if not installed you will see a "no packages found" message.

Using Setuptools allows you to install the tool using the command `python3 mr_software/setup.py install` followed by `chmod +x mrsoftware`. You can now check for installation by running `mrsoftware --help`. Congratualtions, you have now installed `mrsoftware`. 

#### Troubleshooting Installation
* You may have to add the installation directory to the path (or install to a directory on the path) using `--install-dir` option in `easy_install` (such as `easy_install . --install-dir /usr/bin`). 
* On DataHub, you may not have write permissions, so run `python3 setup.py install --user` then look for the directory where the tool got installed to. In my case, it was `/home/{me}/.local/bin`, so now I can run `/home/{me}/.local/bin/mrsoftware` to use our tool.

# Basic usage

### MR Software options

`mrsoftware` has 4 options. 
##### 3 required inputs: 
* `-g` followed by a reference genome in `.fa` format. 
* `-m` followed by a file containing motif information in a `.meme` format.
* and lastly a file containing peak information in a `.bed` format.

##### 1 optional input:
* `-o` followed by an output `.txt` file to write output to. If not specified, the tool prints results to the terminal.

## File format

The output file format is a `.txt` with the name of the motif and the p-value for enrichment.

The basic usage of `mrsoftware` is:
```shell
mrsoftware -m < .meme file > -g < .fa file > < .bed file > -o < output file >
``` 

### Example

To run `mrsoftware` on a small test example (using files in this repository):
* download the [GRCm38.fa reference genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/) (this file is very large ~2.6GB)
* copy GRCm38.fa to /MR_SOFTWARE/tests/
* run: 
* Go to `/mr_software/` directory
```shell
mrsoftware -m ./tests/short.meme -g ./tests/GRCm38.fa ./tests/short.bed -o outfile.txt
```

In this example, we choose to specify an output file using `-o outfile.txt` to store results in the file called `outfile.txt`.

Our program uses randomly collected background sequences so your output may not be exactly the same. When we ran the example command our output was:
```
MOTIF: PO5F1_MOUSE.H11MO.0.A	p-value:3.6991452376732615e-17
MOTIF: ATF4_MOUSE.H11MO.0.A	p-value:0.0001458113625425742
MOTIF: AP2B_MOUSE.H11MO.0.B	p-value:0.006532207010200117
MOTIF: AP2A_MOUSE.H11MO.0.A	p-value:0.05870136798384332
MOTIF: AP2C_MOUSE.H11MO.0.A	p-value:0.12041306253096071
MOTIF: BACH2_MOUSE.H11MO.0.A	p-value:0.12041306253096071
MOTIF: BARX1_MOUSE.H11MO.0.C	p-value:0.12041306253096071
MOTIF: ARI5B_MOUSE.H11MO.0.C	p-value:0.4969325153374233
MOTIF: BACH1_MOUSE.H11MO.0.C	p-value:0.4969325153374233
MOTIF: AHR_MOUSE.H11MO.0.B  p-value:1.0
MOTIF: AIRE_MOUSE.H11MO.0.C p-value:1.0
MOTIF: ALX1_MOUSE.H11MO.0.B p-value:1.0
MOTIF: ANDR_MOUSE.H11MO.0.A	p-value:1.0
MOTIF: ARNT_MOUSE.H11MO.0.B	p-value:1.0
MOTIF: ASCL1_MOUSE.H11MO.0.A	p-value:1.0
MOTIF: ASCL2_MOUSE.H11MO.0.C	p-value:1.0
MOTIF: ATF1_MOUSE.H11MO.0.B	p-value:1.0
MOTIF: ATF2_MOUSE.H11MO.0.A	p-value:1.0
MOTIF: ATF3_MOUSE.H11MO.0.A	p-value:1.0
MOTIF: ATOH1_MOUSE.H11MO.0.B	p-value:1.0
```

# Comparison to HOMER's findMotifsGenome.pl
Our tool recognizes motifs from peaks and determines if there is enrichment. We were inspired to make our tool by [HOMER](http://homer.ucsd.edu/homer/) which was also created at UCSD. Although the tools have very different levels of capabilities, HOMER has a `findMotifsGenome.pl` that we tried to replicate in our `mrsoftware`. Unlike Homer, out tool does not perform *de novo* motif finding.
To compare to the output of `findMotifsGenome.pl`, run:
* `findMotifsGenome.pl <peaks.bed> <reference genome.fa> <output directory> -mask -size 100`
    * Note: `findMotifsGenome.pl` takes over 20 minutes on our computer to run so we recommend using `nohup`. 
What follows is an example of how to run `findMotifsGenome.pl` from the `/mr_software/` directory on the same data:
```shell
mkdir HOMER

findMotifsGenome.pl ./tests/short.bed ./tests/GRCm38.fa HOMER -mask -size 100
```
One can compare the outputs of HOMER's `findMotifsGenome.pl` and our `mrsoftware` by comparing the p-value for enrichment for a given motif.


# Contributors

This repository was generated by Tasha Nguyen and Jonathan Narita, following the [example project](https://github.com/gymreklab/cse185-demo-project) for CSE 185 Spring 2023 made by Melissa Gymrek. Inspired by HOMER.