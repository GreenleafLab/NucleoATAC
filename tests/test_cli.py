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
        cmd = "nucleoatac run --bam example/example.bam --bed example/example.bed --fasta example/sacCer3.fa --out example/test_example --cores 2" 
        args = self.parser.parse_args(cmd.split()[1:])
        main(args)
    def test_occ(self):
        cmd = "nucleoatac occ --bam example/example.bam --bed example/example.bed --fasta example/sacCer3.fa --out example/test_example --cores 2" 
        args = self.parser.parse_args(cmd.split()[1:])
        main(args)
    def test_vprocess(self):
        cmd = "nucleoatac vprocess --out example/test_example --sizes example/test_example.nuc_dist.txt" 
        args = self.parser.parse_args(cmd.split()[1:])
        main(args)
    def test_nuc(self):
        cmd = "nucleoatac nuc --bam example/example.bam --bed example/example.bed --fasta example/sacCer3.fa --out example/test_example --vmat example/test_example.VMat --cores 2" 
        args = self.parser.parse_args(cmd.split()[1:])
        main(args)
    def test_merge(self):
        cmd = "nucleoatac merge --occpeaks example/test_example.occpeaks.bed.gz --nucpos example/test_example.nucpos.bed.gz --out example/test_example" 
        args = self.parser.parse_args(cmd.split()[1:])
        main(args)
    def test_nfr(self):
        cmd = "nucleoatac nfr --bed example/example.bed --occ_track example/test_example.occ.bedgraph.gz --calls example/test_example.nucmap_combined.bed.gz --out example/test_example --bam example/example.bam --fasta example/sacCer3.fa"
        args = self.parser.parse_args(cmd.split()[1:])
        main(args)
 
 











