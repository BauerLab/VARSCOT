VARSCOT
=======

How to build the project
------------------------

```sh
   $ git clone git@github.com:cpockrandt/VARSCOT.git --recursive
   $ mkdir VARSCOT-build && cd VARSCOT-build
   $ cmake ../VARSCOT -DCMAKE_BUILD_TYPE=Release
   $ make -j
```

Documentation
-------------

All binaries document possible flags that you can get with --help, e.g. './search --help'

How to run it
-------------

First you need to create an index of a FASTA-file (e.g. the genome). Currently only DNA5 is supported. The fasta file can hold up to 65536 sequences but the total length cannot exceed 4 giga bases.

```sh
   ./build_index -G file_to_genome.fasta -I path/to/index/indexname
```

The index is built using secondary memory. If you're getting a runtime error, you're most likely runing out of disk space or quota. You can change the TMPRDIR environment variable TMPRDIR on UNIX systems (and TEMP on Windows).

```sh
   $ export TMPDIR=/somewhere/else/with/more/space
```

After building the index, you can run the search executable:

```sh
   $ ./search -G path/to/index/indexname -R reads.fasta .....
```

For a complete list of all arguments, please check the help function './search --help'.
