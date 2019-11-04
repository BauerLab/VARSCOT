# Variant-aware scoring of off-targets for CRISPR Genome Engineering (VARSCOT)

This repository contains the stand-alone VARSCOT program, a tool that takes sequence variants of individuals for CRISPR/Cas9 off-target search and scoring into account.
VARSCOT can be installed via the included Dockerfile.

Tools or libraries used by VARSCOT are the TUSCAN on-target activity predictor and the SeqAn library, can be installed and run via the Dockerfile. See _Building and Running VARSCOT with Docker_. The Dockerfile can also be used as a reference for dependencies and build commands.

## Installation

The simplest approach is to create a Docker image and run VARSCOT as a Docker container as described in _Building and Running VARSCOT with Docker_ below.

However, to install VARSCOT manually, the ```Dockerfile``` can be consulted for detailed steps, but the essential steps follow. Linux would be the most appropriate operating system; at least a bash shell, gcc (5.4.0 or higher), cmake, wget and git are assumed.

* Install Python 2.7 and dependencies:
  * numpy, scipy, sklearn 0.19.0, pybedtools
* Install R 3.4.2 and dependencies:
  * randomForest
* git clone this repository
	* the cloned directory is referred to hereafter as $APP
* mkdir $APP/lib
* cd $APP/lib
* Install Seqan
	* git clone --branch seqan-v2.4.0rc2 --depth 1 https://github.com/seqan/seqan.git
* Install TUSCAN
	* cd $APP/lib
    * wget --no-verbose https://github.com/BauerLab/TUSCAN/archive/master.zip -O tuscan-master.zip
	* unzip tuscan-master.zip
	* rm tuscan-master.zip
	* mv TUSCAN-master TUSCAN
* Install VARSCOT
    * mkdir -p $APP/build/read_mapping_build
    * cd $APP/build/read_mapping_build
	* cmake ../../read_mapping -DCMAKE_BUILD_TYPE=Release
	* make
    * mkdir -p $APP/build/variant_processing_build
	* cd $APP/build/variant_processing_build
	* cmake ../../variant_processing -DCMAKE_BUILD_TYPE=Release
	* make

The next section describes how to use VARSCOT having installed it.

## Usage

Please note that for the user-defined reference genome already a corresponding index must be present.
To build the index run:

```
$ ./build/read_mapping_build/bidir_index -G /path/to/ref/genome.fa -I /path/to/index/prefix
```

The index is built using secondary memory.
If you're getting a runtime error, you're most likely runing out of disk space or quota.
You can change the TMPRDIR environment variable TMPRDIR on UNIX systems (and TEMP on Windows).

```
$ export TMPDIR=/somewhere/else/with/more/space
```

To use VARSCOT run:

```
Example with variant information
$ VARSCOT -f /path/to/variants.vcf -b /path/to/ontargets.bed -o /path/to/output.txt -g /path/to/genome.fa -i /path/to/index/prefix

Example without variant information
$ VARSCOT -f -b /path/to/ontargets.bed -o /path/to/output.txt -g /path/to/genome.fa -i /path/to/index/prefix


Options: -f,  --vcf               Path to input variant file (.vcf), optional argument.
                                  The VCF must contain at least one sample.
         -s,  --sample            Column index of the sample in VCF file to use for off-target search, optional argument.
                                  The index must be 0-based (e.g. first sample has index 0).
                                  If a VCF file is given as input and -s option is missing than first sample in file is chosen as default.
         -b,  --bed               Path to on-target file (.bed).
                                  The on-target information must be present in BED6 format (chr - start - end - name - score - strand) where start/end is 0-based with half open intervals.
                                  Only 23 bp sequences are accepted.
         -o,  --output            Path to output file (.txt).
         -g,  --genome            Path to reference genome (.fa/.fasta)
         -i,  --index             Path to index (see above instructions for building the index).
         -m,  --mismatch          Maximum number of mismatches (default 8).
                                  At maximum 8 mismatches can be detected which includes mismatches at the 21. position (N-position of the PAM).
         -t,  --threads           Number of threads to use (default 8).
         -p,  --pam               Additional non-canonical PAM that should be allowed for off-target search besides (N)GG and (N)GA (default), provide only two letters (upper case, e.g. AG), optional argument
         -T,  --temp-dir          temporary directory (default temp_files)
         -e,  --evaluation        Type of activity scoring, either 'mit', 'class' or 'prob' (default 'mit').
                                  The default evaluation returns the MIT score based on 'http://crispr.mit.edu/'.
                                  In addition, a classification by a Random Forest ML model can be chosen.
                                  Either the exact classes ('class') or the probability of being active ('prob') can be included into the output.
         -v,  --verbose           Keep all temporary files.
         -h,  --help              Help message.
```

Please note: VARSCOT finds off-targets for Cas9 with canonical PAM 5'-NGG-3' or non-canonical PAM 5'-NGA-3'.
Off-targets that include one or more variants are tagged accordingly in the output.
You will need at least as much RAM as the index is large, so please consider running this program on an adequate machine or cluster.

## Building and Running VARSCOT with Docker

### Building
On a host with Docker installed, in a shell, change to the directory in which the Dockerfile exists and type:

    docker build -t varscot .

This will create a Docker image that can be used to run VARSCOT in a container. Only if the Dockerfile or the repositories it refers to change does this build step need to be repeated.

### Running
To run VARSCOT as a Docker container, a command of the following form is required:

    docker run -v $GENOME_DIR:/genome \
               -v $BED_DIR:/bed \
               -v $VCF_DIR:/vcf \
               -v $INDEX_DIR:/index \
               -v $OUTPUT_DIR:/output \
               varscot <VARSCOT-parameters-as-described-above>

The -v options map a directory (e.g. VCF file directory) on your machine to a named volume on the running Docker container. The environment variables may be defined as shown below (the export command assumes you are using a bash shell) or replaced with suitable paths.

    export GENOME_DIR=<location-of-genome-directory>
    export BED_DIR=<location-of-bed-directory>
    export VCF_DIR=<location-of-vcf-directory>
    export INDEX_DIR=<location-of-index-directory>
    export OUTPUT_DIR=<location-of-output-directory>

Note that ```varscot``` must be lower case in the ```docker run``` command above. VARSCOT parameters must refer to the container's volume names when referring to files, as shown in the examples below.

Othe than what has already been noted in this section that is specific to Docker, information provided in _Usage_ section still applies.

### Examples
#### VARSCOT Usage
	docker run varscot -h
	VARSCOT [ARGUMENTS]

	    Arguments:
	    -f,  --vcf               path to input variant file (.vcf)
	    -s,  --sample            column index of the sample in VCF file to use (0-based, default 0)
	    -b,  --bed               path to ontarget file (.bed, must be a bed6 file)
	    -o,  --output            path to output file (.txt)
	    -g,  --genome            path to reference genome (.fa/.fasta)
	    -i,  --index             path to bidirectional index with index prefix
	    -m,  --mismatch          maximum number of mismatches (default 6)
	    -t,  --threads           number of threads to use (default 8)
        -p,  --pam               additional non-canonical PAM that should be allowed for off-target search besides (N)GG and (N)GA (default), provide only two letters (upper case)
	    -e,  --evaluation        either 'mit', 'class' or 'prob' (default 'mit')
	    -v,  --verbose           keep all temp files

#### Running pipeline with VCF
    docker run -v $VCF_DIR:/vcf \
               -v $BED_DIR:/bed \
               -v $GENOME_DIR:/genome \
               -v $INDEX_DIR:/index \
               -v $PWD:/output \
               varscot \
               -f /vcf/example.vcf \
               -b /bed/example.bed \
               -g /genome/hg19.fa \
               -i /index/index \
               -m 2 \
               -o /output/output.txt

#### Running pipeline without VCF
    docker run -v $BED_DIR:/bed \
               -v $GENOME_DIR:/genome \
               -v $INDEX_DIR:/index \
               -v $PWD:/output \
               varscot \
               -b /bed/example.bed \
               -g /genome/hg19.fa \
               -i /index/index \
               -m 2 \
               -o /output/output.txt

#### Running commands within the VARSCOT Docker container
To run commands (e.g. ```bidir_index```) within the VARSCOT Docker container, it is possible to start an interactive shell in a container. A command like this isn't quite right:

    docker run -it varscot /bin/bash

since running the ```varscot``` Docker image corresponds to running the ```VARSCOT``` script. This command would yield a ```VARSCOT``` help message.

What is needed instead is to specify a environment in which the VARSCOT script would be run (within the container). So, given the tail end of ```docker build``` output such as this:

	...
	Removing intermediate container 7a6dcc43930b
	Step 53/56 : VOLUME /output
	---> Running in c299829c679b
	---> 693278a1f70f
	Removing intermediate container c299829c679b
	Step 54/56 : WORKDIR /app
	---> 7f61d3c61fbc
	Removing intermediate container 597d23c886d8
	Step 55/56 : COPY ./entrypoint.sh .
	---> be019493a15d
	Removing intermediate container 9aed2ff91226
	Step 56/56 : ENTRYPOINT /app/entrypoint.sh
	---> Running in a34cbdcd1e26
	---> a995e772b72b
	Removing intermediate container a34cbdcd1e26
	Successfully built a995e772b72b

step 54 of 56 would be suitable since this immediately precedes the point at which VARSCOT (via entrypoint.sh) would be run. The numeric ID ```7f61d3c61fbc``` on the line after "Step 54/56 ..." can then be used instead of the image name ```varscot```:

	docker run -it 7f61d3c61fbc /bin/bash

The result will be an interactive bash shell prompt such as this:

	root@a51dd1386554:/app#

at which commands may be issued. Note that this numeric ID will be different for each ```docker build``` run.

Volumes may also need to be mapped for particular commands, e.g.

	docker run -v $GENOME_DIR:/genome -v $INDEX_DIR:/index -it 7f61d3c61fbc /bin/bash

A ```bidir_index``` command could then be run, e.g.

	./build/read_mapping_build/bidir_index -G /genome/genome.fa -I /index/prefix

where ```prefix``` is the prefix of each file in the index directory.

As per the _Usage_ section setting TMPDIR may be necessary. The ```docker run``` command would need to include a temporary directory volume, e.g.

        docker run -v $ROOT:/bed -v $ROOT:/genome -v $ROOT/bidir-index:/index -v /somewhere/else/with/more/space:/tmpdir -it 7f61d3c61fbc /bin/bash

then in the container before the ```bidir_index``` command:

        export TMPDIR=/tmpdir

## Output
The output of VARSCOT consists of off-targets detected for each target site.
Off-targets are numbered based on the corresponding on-target, but numbers don't imply a ranking of any sort.
Positional information is reported as well as sequence, score, mismatch number and positions.
Positions in general refer to the reference genome, especially if they were not detected within variant regions (marked by tag REF).
However, please note that positions of off-targets that are detected within variant regions might be slightly shifted.
This holds for off-targets containing variants (marked by VAR_chr_start) which start within a variant.
In those cases it was not possible to find a valid position relative to the reference as those off-targets only exist due to a variant and thus never would match any reference position.
