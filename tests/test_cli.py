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
        cmd = "nucleoatac run --bam example/example.bam --bed example/example.bed --fasta example/sacCer3.fa --out example/test_example" 
        args = self.parser.parse_args(cmd.split()[1:])
        main(args)


