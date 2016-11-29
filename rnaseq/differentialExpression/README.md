## Synopsis

Wrapper scripts for running different RNA-Seq differential expression analysis tools.

## Code Example

Rscript differentialExpression_Bl.R

Or, as some of the tools require long runtimes, you can use nohup or screen.

nohup Rscript differentialExpression_Bl.R > log_Bl 2>&1 &

## Expected directory input

Expects an input directory called 'em' that contains expression matrices, organized into subdirectories named by unit type -- e.g. 'countsGn', 'fpkmTx', etc.

Each unit type subdirectory should contain expression matrix files, typically in the format 'RaEm_expressionMatrix.txt', where Ra = Read aligner and Em = Expression modeler.

Differential expression analysis outputs are directed to the directory 'de'.

## Other parameters

Most of these scripts require an annotation conversion to append gene symbols to the outputs. It should be a tab-delimited text file containing 3 columns: Transcript ID, Gene ID, and Gene Symbol. These must match those used in the indices and in the reference files.

## Contributors

Claire Williams, Alyssa Baccarella, Charlie Kim

## License

See LICENSE in the top level directory.
