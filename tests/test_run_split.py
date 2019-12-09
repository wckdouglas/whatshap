"""
Integration tests that use the command-line entry points run_split
"""
from tempfile import TemporaryDirectory
import os
from collections import namedtuple
from collections import defaultdict

from pytest import raises, fixture, mark
import pysam
from whatshap.cli.haplotag import run_haplotag
from whatshap.cli.split import run_split

def test_split():
	with TemporaryDirectory() as tempdir:
		outlist = tempdir + '/outlist.txt'
		outbam = tempdir + '/outbam.bam'
		# produce a list of read assignments using haplotag
		run_haplotag(variant_file='tests/data/haplotag_poly.vcf.gz', alignment_file='tests/data/haplotag_poly.bam', ploidy=4, output=outbam, haplotag_list=outlist)
		# use list as input for split
		outbam_h1 = tempdir + '/outbamh1.bam'
		outbam_h2 = tempdir + '/outbamh2.bam'
		outbam_h3 = tempdir + '/outbamh3.bam'
		outbam_h4 = tempdir + '/outbamh4.bam'
		run_split(reads_file=outbam, list_file=outlist, output_h1=outbam_h1, output_h2=outbam_h2, output_h3=outbam_h3, output_h4=outbam_h4)

		expected_splits = {3 :'S1_31286_NA19240_HAP2', 0 : 'S1_248595_HG00514_HAP1', 2 : 'S1_284251_NA19240_HAP1', 1 : 'S1_103518_HG00514_HAP2'}
		for hap, outfile in enumerate([outbam_h1, outbam_h2, outbam_h3, outbam_h4]):
			nr_reads = 0
			for alignment in pysam.AlignmentFile(outfile):
				assert expected_splits[hap] == alignment.query_name
				nr_reads += 1
			assert nr_reads == 1

	

