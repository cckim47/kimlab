## Synopsis

A script that calculates performance, as estimated by recall and precision against reference datasets.

## Code Example

./calculatePerformance.pl

Runtime is generally short (minutes), and doesn't require nohup or screen.

## Expected directory input

Relative to the parent directory, 2 sets of directories are needed.
  * At the same level as the parent directory, a directory called 'de' that also contains the specified subdirectories in the array '@unitTypes' -- e.g. 'fpkmGn', 'countsTx'
  * Within the parent directory, a subdirectory called 'references' that contains the symbol lists. The filename prefixes and suffixes are hard-coded, but can be edited along with the reference names in the array '@references'.

## Contributors

Claire Williams, Alyssa Baccarella, Charlie Kim

## License

See LICENSE in the top level directory.
