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
	* git clone https://github.com/seqan/seqan.git
	* cd $APP/lib/seqan
	* git checkout develop
    * git checkout d7a4805454dd2183e34a7e2c776485ce4ada8732
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

The next section describes usage how to use VARSCOT having installed it.

## Usage

Please note that for the user-defined reference genome already a corresponding index must be present.
To build the classifier run:

```
$ build/read_mapping_build/bidir_index -G /path/to/ref/genome.fa -I /path/to/index/prefix
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
         -i,  --index             Path to index with index prefix (see above instructions for building the index).
         -m,  --mismatch          Maximum number of mismatches (default 8).
                                  At maximum 8 mismatches can be detected which includes mismatches at the 21. position (N-position of the PAM).
         -t,  --threads           Number of threads to use (default 8).
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

This will create a Docker image that can be used to run VARSCOT in a container.

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
#### With VCF
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

#### Without VCF
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

#### Getting Usage
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
	    -e,  --evaluation        either 'mit', 'class' or 'prob' (default 'mit')
	    -v,  --verbose           keep all temp files
