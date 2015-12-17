from unittest import TestCase
from nucleoatac.cli import nucleoatac_parser, main
import os
from pkg_resources import resource_filename

class NucleoATACTestCase(TestCase):
    """
    Base TestCase class for nucleoatac commands
    """
    def setUp(self):
        self.parser = nucleoatac_parser()
    def test_run(self):
        cmd = "nucleoatac run --bam example/example.bam --bed example/example.bed --fasta example/sacCer3.fa --out example/test_results/test --cores 2" 
        args = self.parser.parse_args(cmd.split()[1:])
        main(args)
    def test_occ(self):
        cmd = "nucleoatac occ --bam example/example.bam --bed example/example.bed --fasta example/sacCer3.fa --out example/test_results/test --cores 2" 
        args = self.parser.parse_args(cmd.split()[1:])
        main(args)
    def test_vprocess(self):
        cmd = "nucleoatac vprocess --out example/test_results/test --sizes example/example_results/example.nuc_dist.txt" 
        args = self.parser.parse_args(cmd.split()[1:])
        main(args)
    def test_nuc(self):
        cmd = "nucleoatac nuc --bam example/example.bam --bed example/example.bed --fasta example/sacCer3.fa --out example/test_results/test --vmat example/example_results/example.VMat --cores 2" 
        args = self.parser.parse_args(cmd.split()[1:])
        main(args)
    def test_merge(self):
        cmd = "nucleoatac merge --occpeaks example/example_results/example.occpeaks.bed.gz --nucpos example/example_results/example.nucpos.bed.gz --out example/test_results/test" 
        args = self.parser.parse_args(cmd.split()[1:])
        main(args)
    def test_nfr(self):
        cmd = "nucleoatac nfr --bed example/example.bed --occ_track example/example_results/example.occ.bedgraph.gz --calls example/example_results/example.nucmap_combined.bed.gz --out example/test_results/test --bam example/example.bam --fasta example/sacCer3.fa"
        args = self.parser.parse_args(cmd.split()[1:])
        main(args)
 
 











