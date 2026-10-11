plastweaver

Plastome assembly, ordering, annotation and GenBank submission from paired-end reads.

For every library, plastweaver chooses a seed and a reference (SeedRef), assembles the plastome from the reads (AssemblerComplete), orders the circle as LSC-IRa-SSC-IRb, annotates it (PlastomeAnnotator) and writes gene tables, a circular genome map, an IR junction plot and the files for a GenBank submission (PlastomeExport). The sequence always comes from the reads; the references only guide the assembly and the annotation.

Requirements
Julia 1.10 or newer (https://julialang.org/downloads/)
Linux or macOS, bash, tar
Internet access to NCBI while running: SeedRef downloads the reference plastome chosen for each library (and, at the first run, the marker panel)
Free disk space for the reference dataset (about 3.2 GB unpacked, plus the downloaded archive while unpacking), the reads and the results
Installation
Get the code:
git clone https://github.com/kubek78/plastweaver.git
cd plastweaver
(or Code -> Download ZIP on this page, unpack it and go into the folder)
Install the Julia dependencies (once):
julia --project=. -e 'using Pkg; Pkg.instantiate()'
Download the reference dataset (one file, plastweaver_dataset_refseq237.tar.gz, 1.4 GB; NCBI RefSeq release 237) from Zenodo: https://doi.org/10.5281/zenodo.23303797 Check the file (optional) and unpack it into the refseq_plastid folder of plastweaver:
sha256sum plastweaver_dataset_refseq237.tar.gz
# 5f634f31a9150cb2c42a6d44b91ab4b0725dd666cd3b8d961c5ea67bdb5bdf69
tar -xzf plastweaver_dataset_refseq237.tar.gz -C refseq_plastid
Afterwards refseq_plastid/dataset/ must contain genomes.fasta.gz, records.tsv, records_final.tsv, labels.tsv, features.tsv, ir_table.tsv, ref_sketches.bin, refindex_all/ and models_final/ (about 3.2 GB). The downloaded .tar.gz file can then be deleted.
Running

Put the paired reads of each library in one folder, two files per library:

reads/LIBRARY1_R1.fastq.gz   reads/LIBRARY1_R2.fastq.gz
reads/LIBRARY2_R1.fastq.gz   reads/LIBRARY2_R2.fastq.gz

Then, from the plastweaver folder:

julia plastweaver.jl run --input reads

Options:

--samples LIB1,LIB2     only these libraries
--results DIR           output folder (default: results)
--threads N             threads per step (default: auto)
--resume                keep the assemblies that already exist in the results folder
Results

One folder per library, results/LIBRARY/:

LIBRARY_plastome.fasta        ordered plastome(s), LSC-IRa-SSC-IRb
LIBRARY_plastome.gff3         annotation (coordinates of the FASTA above)
LIBRARY_plastome_regions.tsv  LSC, IR1 (= IRa), SSC, IR2 (= IRb) coordinates
genes/                        gene tables by function, intron table, gene summary
figures/                      circular genome map with GC content, IR junction plot
genbank/                      metadata.tsv to fill in, README.txt (submission files after
                              the genbank step, see below)
seed/                         seed, reference and the assembler settings used
assembly/                     all assembler outputs (assembly graph, reports)
plastweaver.log

results/plastweaver_summary.tsv has one row per library: the primary record, how its circle was closed (closure) and whether it can be submitted to GenBank as a complete circular genome (genbank_ready).

Orientation of the ordered plastome: the IR that follows the LSC is IRa; the SSC starts with ndhF next to IRa; the molecule is oriented so that rpoB lies on the + strand (when rpoB cannot be found, the orientation of the nearest RefSeq plastome is used).

Closure of the circle: validated, ir_collapse, read_bridge and graph_path are conclusive closures of the assembler; other values (ref_scaffold, inferred_join, topology_inferred, estimated_repeat, gap_unknown) mean that the circle was closed with less evidence and should be checked. Such records get genbank_ready = no.

GenBank submission
Fill in results/LIBRARY/genbank/metadata.tsv (organism and source information; instructions in results/LIBRARY/genbank/README.txt).
Run
julia plastweaver.jl genbank results/LIBRARY
Check genbank/report_records.tsv and genbank/report.tsv, then submit submission.fsa, submission.tbl and source_modifiers.tsv through the GenBank Submission Portal (or validate them first with NCBI table2asn).

Gene tables, figures and the metadata template can be written again from an existing assembly with

julia plastweaver.jl export results/LIBRARY
Settings

Assembler settings (number of read pairs, k-mer size, thresholds) are in config_template.txt. Gene categories and map colours: gene_categories.tsv; CDS product names for GenBank: gene_products.tsv.

License

GPL-3.0 (see LICENSE). The reference dataset (https://doi.org/10.5281/zenodo.23303797) is derived from NCBI RefSeq plastid genomes (release 237).
