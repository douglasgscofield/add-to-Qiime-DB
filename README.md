Qiime-DB-enhance
================

A small collection of scripts to ease the enhancement of a Qiime reference database.  BioPerl is required, as is a local copy of the NCBI taxonomy database (see below).  Two of these scripts started life written by others and are extensively modified here.

The workflow requires a few steps and is quite simple.  We start with a set of ITS sequences to blast, in Fasta format, say in `seqs.fa`.  We will have blast results returned in a table format that includes taxonomic information.

```bash
blastn -db nt -query seqs.fa -outfmt "6 std staxids sscinames sskingdoms sblastnames" > seqs.bl6
```

Run these blast results through a small filter to extract the gi, GenBank, and taxon IDs in a simple format.  The `qiime_get_blast_ids_for_genbank.pl` script depends upon the `-outfmt` string for the blast being as you see it above.

```bash
qiime_get_blast_ids_for_genbank.pl seqs.bl6 > seqs.ids
```

Fetch the GenBank sequences corresponding to these taxon IDs.  The `qiime_get_genbank_seqs.pl` script was originally the BioPerl script [`bp_download_query_genbank.pl`](https://github.com/bioperl/bioperl-live/blob/master/scripts/utilities/bp_download_query_genbank.pl), and the initial commit of the script to this repository was with a copy of that script so modifications can be tracked.

The output of this script is two files, one of sequences to STDOUT and one of IDs and GenBank accession numbers  to STDERR.  Both are used in the next step.  The content of the accession number file should be the same as the `seqs.ids` file, but the entries are likely to be in a different order that matches the order of sequences in the sequence file.

```bash
qiime_get_genbank_seqs.pl --gifile seqs.ids > gb_seqs.fa 2> gb_ids_acc.txt
```

The sequence names and IDs are formatted to be helpful when searching for taxonomic information in the following step.  The `qiime_get_taxonomy_from_seqs.pl` script originally started as the [`taxonomy.pl`](https://github.com/hyphaltip/mobedac-fungi/blob/master/scripts/taxonomy.pl) script in the [MOBeDAC Fungi Database](https://github.com/hyphaltip/mobedac-fungi) repository, and the initial commit of the script to this repository was with a copy of that script so modifications can be tracked.

This script also requires a local copy of the NCBI taxonomy database in an `ncbi/` directory within your working directory.  This is available from the [NCBI Taxonomy Database FTP site](ftp://ftp.ncbi.nih.gov/pub/taxonomy), as file `taxdump.tar.gz`, `taxdump.tar.Z`, or `taxdmp.zip`; download whichever is convenient for your operating system and file decompressor.  For `taxdump.tar.gz`, you can fetch and unpack the latest version into the `ncbi/` directory with:

```bash
mkdir ncbi
cd ncbi
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xvzf taxdump.tar.gz
cd ..
```

**Note**: This next step might challenge the memory capacity of desktop computers because of the reading of the taxonomic database.

We then use this database to find full taxonomic information associated with the retrieved sequences and put that into a format useful for Qiime.  The `--exclude` is optional, and can contain any regular expression that would match taxonomic information for taxa to be excluded from the final results.

```bash
qiime_get_taxonomy_from_seqs.pl --exclude "uncultured fungus" --accessionfile gb_ids_acc.txt gb_seqs.fa
```

The output of this script is in several files.  For a set of sequences called `gb_seqs.fa` the files called `gb_seqs.taxonomy.fa` and `gb_seqs.taxonomy` are the final files to present to Qiime.

`gb_seqs.taxonomy.fa`: 
sequences from `gb.seqs.fa` but with a modified identifier and a description that contains full taxonomic information in a Qiime-compatible format.  **Unidentified and excluded sequences do not appear here.**

`gb_seqs.taxonomy`:
one line per sequence, containing identifiers and full taxonomic information in the same order as in `gb_seqs.taxonomy.fa`.  **Unidentified and excluded sequences do not appear here.**

`gb_seqs.unidentified`:
one line per unidentified sequence.  If taxonomic information was not found in the database or was poorly formatted, the identifiers of the sequences appear here.

`gb_seqs.excluded`:
one line per excluded sequence.  Sequences appear here because their taxonomic information matched the regular expression provided with the `--exclude` option.
