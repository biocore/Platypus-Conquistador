Platypus Conquistador Tutorial
==============================

The files for this tutorial are from [this](http://metagenomics.anl.gov/metagenomics.cgi?page=MetagenomeOverview&metagenome=4502816.3) study that was published [here - in review](empty). These files have been randomly subsample to create the tutorial and have them run in any computer; same for the [Integrated Microbial Genomes](http://img.jgi.doe.gov/) reference database.

Download the [tutorial](ftp://thebeast.colorado.edu/pub/platypus/tutorial_example.tgz) files.

Untar the files and move to the new created folder:
`tar zxvf tutorial_example.tgz; cd tutorial`

Now we need to split our database in the interest/rest databases. For this example, we will look for any occurrences of the word salmonella:
`platypus_split_db.py -t bacteria.contigs.txt -f bacteria.contigs.fna -q 'salmonella' -o databases`

The next step is to create a basic [QIIME](http://qiime.org) mapping file for each input file to QC them. Note that this is not a necessary step but will help us remove some problematic and noisy sequences. Additionally, that in this example we are using the already QC-ed sequences from [MG-RAST](http://metagenomics.anl.gov/) but this step still remove some low quality seqs. 
``` bash
mkdir maps
for i in `ls input_samples/`; do printf "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\n%s\t\t\t%s\n" $i $i > maps/${i/.fna}.txt; done
for i in `ls input_samples/`; do split_libraries.py -f input_samples/${i} -l 75 -b 0 -m maps/${i/.fna}.txt -o splitLib/${i/.fna} -c -p; done
cat splitLib/*/seqs.fna > splitLib/seqs.fna
```

After this we will have a single seqs.fna file with all the QC-ed files so we are ready to BLAST with relaxed parameters agains the interest and rest database. In this example, we are going to use 2 jobs for the processing.
``` bash
mkdir blastOut
parallel_blast.py -O 2 -n 10 -w 11 -e .1 -r databases/rest.fna -i splitLib/seqs.fna -o blastOut/rest
parallel_blast.py -O 2 -n 10 -w 11 -e .1 -r databases/interest.fna -i splitLib/seqs.fna -o blastOut/interest
```

Now we should have 2 files with the results of BLASTing against the interest and the rest. Consequently, we can use these resulting files to compare the results and confirm the presence of Salmonella in our files. In the example command, we are using different levels of similarity: 80,85,90,95,100, and different lengths for the alignment: 80,100.
`platypus_compare.py -i blastOut/rest/seqs_blast_out.txt -j blastOut/interest/seqs_blast_out.txt -p 80,85,90,95,100 -a 80,100 -o compare `

This last step will create several files with the summary of each combination of options (% identity, and alignment length), each summary will have the sequences that matched at any given level and two compiled results: compile_output_no_nohits.txt, results without those sequences that weren't found in the first DB, and compile_output.txt, all results.