#pyatac

The pyatac command line function contains several subcommands.  To run a subcommand, run `pyatac subcommand ...`.  The available subcommands are:

* ins
* bias
* vplot
* bias_vplot
* sizes
* pwm

To see available subcommands, you can also run `pyatac --help`.  For each subcommand, options can be viewed as `pyatac subcommand --help`.


###ins
`pyatac ins` makes a bedgraph file with insertion frequencies for every base.

If a bedfile is provided, insertion track will only span regions in bedfile. Otherwise track will span whole genome.

Output:  A tabix-indexed bedgraph file with insertion densities.

###bias
`pyatac bias` computes a bias track.

If a bedfile is provided, bias track will only span regions in bedfile. Otherwise track will span whole genome.

Output:  A tabix-indexed bedgraph file with relative Tn5 preferences.  These preferences are in log space.

###vplot
`pyatac vplot` generates a vplot around sites in input bed file.

Output:  An eps file with vplot image as well as a text file with vplot data.

###bias_vplot
`pyatac bias_vplot` computes the vplot based on the sequence bias and the global fragment size distribution

Output: An eps file with vplot image as well as a text file with vplot data.

###sizes
`pyatac sizes` computes the fragment size distribution

Output: text file with sizes information

###pwm
`pyatac pwm` computes the relative nucleotide frequencies around sites

Output: text file with PWM information

###nucleotide
`pyatac nucleotide` computes the nucleotide frequencies around sites in the input bedfile

Output: text file with frequencies










