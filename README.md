# fastkmers

_Author: Koki_

Python script for reading and counting the frequency of N sized kmers

The following script was made on python3.7.7
and requires the numpy package (version 1.18.1)

File overview

    .
    └── src
       ├── dev.ipynb  # development notebook
       ├── fastkmers.py  # main script
       └── utils
          ├── __init__.py
          └── decorator.py

The script can be used interactively in a notebook by importing the file

```Python
import fastkmers
```

Or ran in terminal

```bash
fastkmers.py --fastq $fastqfile \
--N $kmersize \
--output $outputfname \
--cores $ncores
```

| Parameters | Description                                                                                                                                                                                |
| :--------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| --fastq    | Full path to fastq file of interest. Currently supporting reading zipped files                                                                                                             |
| --N        | Size of kmers to generate. Note that computationl time can drastically increase depending on the specified size. Please refer to the `dev.ipynb` notebook for an estimation of performance |
| --output   | Name of output file, which will be a tab-delimited file with the kmer and corresponding frequency                                                                                          |
| --cores    | Number of cpu cores to parallelize the task on. Can lead to drastic performance improvements                                                                                               |
