Scripts for augmenting a Qiime database
=====

This is a small collection of scripts to ease the augmentation of a Qiime reference database.  BioPerl is required, as is `sort` and a local copy of the NCBI taxonomy database (see below).  Two of these scripts started life written by others and are extensively modified here.

Blast for GenBank hits 
------

The workflow requires a few steps and is quite simple.  We start with a set of ITS sequences to blast, in Fasta format, say in `seqs.fa`.  We will have blast results returned in a table format that includes taxonomic information.

```bash
blastn -db nt -query seqs.fa -outfmt "6 std staxids sscinames sskingdoms sblastnames" > seqs.bl6
```

Only `"6 std staxids"` is strictly required by the next step but the other columns are useful for checking the taxonomic content of results.

Extract hit information from Blast results
------

Run these blast results through a small filter to extract the target gi ID, GenBank ID, taxon IDs, and start and end of the HSP in the target in a tab-delimited format.  The `qiime_get_blast_ids_for_genbank.pl` script depends upon the `-outfmt` string for the blast being as you see it above.  The `seqs.ids` file produced by this step is used in the two following steps.

```bash
qiime_get_blast_ids_for_genbank.pl seqs.bl6 | sort -k1,2 -u > seqs.ids
```

Note that the `sort -k1,2 -u` command removes redundant blast subject sequences for which the first two columns (the merged gi and taxon IDs and the GenBank accession ID) are identical.  There is no provision for choosing the most appropriate HSP, if you feel there might be some issues with the HSP chosen as a reference for a particular taxon, first check the filtering of the blast results here.

Fetch GenBank sequences for the hits
------

Fetch the GenBank sequences corresponding to these hits.  The `qiime_get_genbank_seqs.pl` script was originally the BioPerl script [`bp_download_query_genbank.pl`](https://github.com/bioperl/bioperl-live/blob/master/scripts/utilities/bp_download_query_genbank.pl), and the initial commit of the script to this repository was with a copy of that script so modifications can be tracked.  The output of this script is to STDOUT.

```bash
qiime_get_genbank_seqs.pl --gifile seqs.ids > gb_seqs.fa
```

Create local copy of NCBI taxonomy database
------

We require a local copy of the NCBI taxonomy database in an `ncbi/` directory within your working directory.  This is available from the [NCBI Taxonomy Database FTP site](ftp://ftp.ncbi.nih.gov/pub/taxonomy), as file `taxdump.tar.gz`, `taxdump.tar.Z`, or `taxdmp.zip`; download whichever is convenient for your operating system and file decompressor.  For `taxdump.tar.gz`, you can fetch and unpack the latest version into the `ncbi/` directory with:

```bash
mkdir ncbi
cd ncbi
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xvzf taxdump.tar.gz
cd ..
```

After you unpack the database, create an index for it.  By default this is created into `ncbi-indices/`.  This will speed up usage of the script considerably.  These indices must be recreated each time a new version of the NCBI taxonomy database is downloaded.

```bash
mkdir ncbi-indices
perl -MBio::DB::Taxonomy -e 'Bio::DB::Taxonomy->new(-source=>"flatfile", -nodesfile=>"ncbi/nodes.dmp", -namesfile=>"ncbi/names.dmp", -directory=>"ncbi-indices", -force);'
```

You might see an error while running this, especially on a Mac:

~~~~
$ perl -MBio::DB::Taxonomy -e 'Bio::DB::Taxonomy->new(-source=>"flatfile", -nodesfile=>"ncbi/nodes.dmp", -namesfile=>"ncbi/names.dmp", -directory=>"ncbi-indices", -force);'
HASH: Out of overflow pages.  Increase page size
perl(3790,0x7fff789a1310) malloc: *** mach_vm_map(size=18446744073704988672) failed (error code=3)
*** error: can't allocate region
*** set a breakpoint in malloc_error_break to debug
~~~~

If you check the `ncbi-indices/` directory, you should still see the indices there (sizes will not be exact):

~~~~
$ ls -l ncbi-indices/
total 844928
-rw-r--r--  1 Douglas  staff   50885040 Dec  8 13:06 id2names
-rw-r--r--  1 Douglas  staff  296517632 Dec  8 13:06 names2id
-rw-r--r--  1 Douglas  staff   37530524 Dec  8 13:05 nodes
-rw-r--r--  1 Douglas  staff   47665152 Dec  8 13:05 parents
~~~~

As far as I have been able to tell, the indices are still correct at this point despite this error.

If you are keeping the database or indices in directories other than `ncbi/` and `ncbi-indices/` and thus will be using the `--db-directory` and/or `--db-index-directory` options to the following script, then directory names will need to be modified when performing the steps above.


Assemble taxonomic hierarchies for GenBank hits
------

The previous scripts formatted sequence names and IDs to be helpful when searching for taxonomic information in this step.  The `qiime_get_taxonomy_from_seqs.pl` script originally started as the [`taxonomy.pl`](https://github.com/hyphaltip/mobedac-fungi/blob/master/scripts/taxonomy.pl) script in the [MOBeDAC Fungi Database](https://github.com/hyphaltip/mobedac-fungi) repository, and the initial commit of the script to this repository was with a copy of that script so modifications can be tracked.

**Note**: If you choose not to create the indices above, this step may have a long initiation time and may generate errors, see above.

If the taxonomic databases or their indices are in directories other than `ncbi/` and `ncbi-indices/`, respectively, their locations may be specified with `--db-directory` and `--db-index-directory`.

We now find full taxonomic information associated with the retrieved sequences and put that into a format useful for Qiime.  The `--exclude` option may be used to exclude specific taxa from the final results.  This option may be used one or more times to specify regular expressions against which the taxonomic hierarchy is compared.  Information about excluded taxa is written to a file.

Two other options of note are `--min-to-truncate` (default 1000), which sets a minimum length (bp) of GenBank target sequence which will be subject to target truncation (shortening around the HSP), and `--min-after-truncate` (default 300), which is the minimum length of the truncation result, expanded equally up- and downstream of the HSP.  Target truncation shortens the sequence to just the region indicated by the start and end positions of the HSP returned by the Blast results above, and then expands the site.  This can be very useful for producing consistently-sized sequences if the ITS hit is within a large (perhaps multi-Mbp) GenBank sequence. 

There are also several other options that might be useful, including a facility for replacing taxonomic hierarchies that are incomplete.  Find out more by using the `--help` option.

```bash
qiime_get_taxonomy_from_seqs.pl --exclude "uncultured fungus" --accessionfile seqs.ids gb_seqs.fa
```

Output files
------

The output of this script is in several files.  For a set of sequences called `gb_seqs.fa` the files are prefixed with `gb_seqs.`.

`gb_seqs.taxonomy.fa` and `gb_seqs.taxonomy` are the final files to present to Qiime.

`gb_seqs.taxonomy.fa` contains sequences derived from `gb_seqs.fa`, potentially shortened to surround the blast HSP and with a modified identifier and a description that contains full taxonomic information in a Qiime-compatible format.  **Unidentified, redundant and excluded sequences do not appear here.**

`gb_seqs.taxonomy` contains taxonomic information for each sequence, in the same order as in `gb_seqs.taxonomy.fa`.  **Unidentified, redundant and excluded sequences do not appear here.**

`gb_seqs.unidentified` contains one line per unidentified sequence.  If taxonomic information was not found in the database or was poorly formatted in the blast results, the identifiers of the sequences will appear here.

`gb_seqs.incomplete` contains one line per sequence with an incomplete taxonomic hierarchy, and one line per sequence for which the hierarchy was replaced following the contents of the `--incompletefile` file.

`gb_seqs.redundant` contains one line per sequence that was excluded from final results for duplicating a taxonomic hierarchy also matched by another sequence with a longer HSP.

`gb_seqs.excluded` contains one line per sequence that was excluded from output because its taxonomic hierarchy matched a regular expression provided with the `--exclude` option.

`gb_seqs.truncated` contains one line per sequence that was truncated from its original size to the HSP following `--min-to-truncate`, and also contains information if the sequence was then expanded to a larger size following `--min-after-truncate`.


Completing incomplete taxonomic hierarchies
------

Within the `.incomplete` file, lines beginning with `HIERARCHY_INCOMPLETE` contain accessions with incomplete taxonomic information, missing one of the major levels of classification.  The `qiime_get_taxonomy_from_seqs.pl` script has an option `--incompletefile` which takes an argument which is a file containing incomplete hierarchies in the first column, and their completed counterparts in the second column.  The `.incomplete` file generated by a run of this pipeline can be used to create such a file with a bit of handwork and taxonomic information.  If the `--incompletefile` option is used, then the `.incomplete` file will additionally contain lines beginning with `HIERARCHY_REPLACED` listing replacements enabled by the file.

For our project, we have produced a set of replacement taxonomic hierarchies for over 300 fungal taxonids, and provide these here in the file `completed-taxonomies-with-ids.20141213-fungi.txt`.  Columns 1 and 2 of this file are the GenBank accession number and the taxonid and columns 3 and 4 are the incomplete and complete hierarchies, respectively.  This file can be provided to the script by using a named pipe to extract columns 3 and 4:

    qiime_get_taxonomy_from_seqs.pl --incompletefile <(cut -f3-4 completed-taxonomies-with-ids.20141213-fungi.txt) ...


