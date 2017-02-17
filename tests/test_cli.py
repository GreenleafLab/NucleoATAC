from unittest import TestCase
from nucleoatac.cli import nucleoatac_parser, nucleoatac_main
from pyatac.cli import pyatac_parser, pyatac_main


class NucleoATACTestCase(TestCase):
    """
    Base TestCase class for nucleoatac commands
    """
    def setUp(self):
        self.parser = nucleoatac_parser()
    def test_run(self):
        cmd = "nucleoatac run --bam example/example.bam --bed example/example.bed --fasta example/sacCer3.fa --out example/test_results/test --cores 2" 
        args = self.parser.parse_args(cmd.split()[1:])
        nucleoatac_main(args)
    def test_occ(self):
        cmd = "nucleoatac occ --bam example/example.bam --bed example/example.bed --fasta example/sacCer3.fa --out example/test_results/test --cores 2" 
        args = self.parser.parse_args(cmd.split()[1:])
        nucleoatac_main(args)
    def test_vprocess(self):
        cmd = "nucleoatac vprocess --out example/test_results/test --sizes example/example_results/example.nuc_dist.txt" 
        args = self.parser.parse_args(cmd.split()[1:])
        nucleoatac_main(args)
    def test_nuc(self):
        cmd = "nucleoatac nuc --bam example/example.bam --bed example/example.bed --fasta example/sacCer3.fa --out example/test_results/test --vmat example/example_results/example.VMat --cores 2" 
        args = self.parser.parse_args(cmd.split()[1:])
        nucleoatac_main(args)
    def test_merge(self):
        cmd = "nucleoatac merge --occpeaks example/example_results/example.occpeaks.bed.gz --nucpos example/example_results/example.nucpos.bed.gz --out example/test_results/test" 
        args = self.parser.parse_args(cmd.split()[1:])
        nucleoatac_main(args)
    def test_nfr(self):
        cmd = "nucleoatac nfr --bed example/example.bed --occ_track example/example_results/example.occ.bedgraph.gz --calls example/example_results/example.nucmap_combined.bed.gz --out example/test_results/test --bam example/example.bam --fasta example/sacCer3.fa"
        args = self.parser.parse_args(cmd.split()[1:])
        nucleoatac_main(args)
 
class PyATACTestCase(TestCase):
    """
    Base TestCase class for pyatac commands
    """
    def setUp(self):
        self.parser = pyatac_parser()
    def test_sizes(self):
        cmd = "pyatac sizes --bam example/example.bam --out example/test_results/test_sizes" 
        args = self.parser.parse_args(cmd.split()[1:])
        pyatac_main(args)
 











