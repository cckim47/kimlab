## Synopsis

Wrapper scripts for running different RNA-Seq read alignment and expression estimation tools.

There are 3 general wrapper types: read aligers only (align_@@.pl), expression modelers only (model_@@.pl), and wrappers that both align and model expression (alignAndModel_@@.pl).

## Code Example

./alignAndModel_KaKa.pl

Or, as these workflows can require long runtimes, you can use nohup or screen.

nohup ./alignAndModel_KaKa.pl > log_KaKa 2>&1 &

## Expected directory input

Within the parent directory in which the wrappers are executed, the wrapper will search for subdirectories containing a fastq file of the same name as the subdirectory. So, for example, a directory "classical01" will be processed if the file "classical01.fastq" exists. There is also a check for whether the sub-subdirectory for the aligner/modeler's output exists -- this permits launching of parallel processes (although a word of caution: this is not intended for highly parallelized use, as there is no job control -- i.e., crashes are possible in very rare instances).

## Other parameters

A number of parameters are hard-coded into the header of each wrapper, and should be modified to suit your environment. These include paths to relevant indices, genome files, and annotation (gtf) files; output directory names; input directory names; and CPU threads, in some cases.

## Contributors

Claire Williams, Alyssa Baccarella, Charlie Kim

## License

See LICENSE in the top level directory.
