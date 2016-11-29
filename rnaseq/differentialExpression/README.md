## Synopsis

Wrapper scripts for running different RNA-Seq differential expression analysis tools.

## Code Example

Rscript differentialExpression_Bl.R

Or, as some of the tools require long runtimes, you can use nohup or screen.

nohup Rscript differentialExpression_Bl.R > log_Bl 2>&1 &

## Expected directory input

Most scripts expect an input directory called 'em' that contains expression matrices, organized into subdirectories named by unit type -- e.g. 'countsGn', 'fpkmTx', etc.

Each unit type subdirectory should contain expression matrix files, typically in the format 'RaEm_expressionMatrix.txt', where Ra = Read aligner and Em = Expression modeler.

A subset of the scripts (for Bs, Cd, and Su) require alignment or expression data in special formats. These scripts will search for the needed files within an input directory called 'mapAndModel' that contains subdirectories for each sample such as 'classical01', each with an Ra/Em subdirectory such as 'HsBs'.

For all, differential expression analysis outputs are directed to the directory 'de'.

## Other parameters

Most of these scripts require an annotation conversion to append gene symbols to the outputs. It should be a tab-delimited text file containing 3 columns: ENST (Transcript ID), ENSG (Gene ID), and geneName (Gene Symbol). These must match those used in the indices and in the reference files.

## Contributors

Claire Williams, Alyssa Baccarella, Charlie Kim

## License

See LICENSE in the top level directory.
