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


def test_split_bam():
	with TemporaryDirectory() as tempdir:
		outlist = tempdir + '/outlist.tsv'
		outbam1 = tempdir + '/outbam.bam'
		# produce list of read assignments to haplotypes
		run_haplotag(variant_file='tests/data/haplotag_1.vcf.gz', alignment_file='tests/data/haplotag.bam', haplotag_list=outlist, output=outbam1)
		outbam_h1 = tempdir + '/outbamh1.bam'
		outbam_h2 = tempdir + '/outbamh2.bam'
		outbam_untagged = tempdir + '/outbamuntagged.bam'
		run_split(reads_file='tests/data/haplotag.bam', list_file=outlist, output_h1=outbam_h1, output_h2=outbam_h2, output_untagged=outbam_untagged)

		expected_splits = {}
		for alignment in pysam.AlignmentFile(outbam1):
			if alignment.has_tag('HP'):
				expected_splits[alignment.query_name] = alignment.get_tag('HP')
			else:
				expected_splits[alignment.query_name] = 0
		assert len(expected_splits) > 0

		total_reads = 0
		for hap, outfile in enumerate([outbam_untagged, outbam_h1, outbam_h2]):
			for alignment in pysam.AlignmentFile(outfile):
				hap_key = hap
				assert expected_splits[alignment.query_name] == hap_key
				total_reads += 1
		assert total_reads == len(expected_splits)
				

def test_split_bam_large():
	with TemporaryDirectory() as tempdir:
		outlist = tempdir + '/outlist.tsv'
		outbam1 = tempdir + '/outbam.bam'
		# produce list of read assignments to haplotypes
		run_haplotag(variant_file='tests/data/haplotag.large.vcf.gz', alignment_file='tests/data/haplotag.large.bam', haplotag_list=outlist, output=outbam1)
		outbam_h1 = tempdir + '/outbamh1.bam'
		outbam_h2 = tempdir + '/outbamh2.bam'
		outbam_untagged = tempdir + '/outbamuntagged.bam'
		run_split(reads_file='tests/data/haplotag.large.bam', list_file=outlist, output_h1=outbam_h1, output_h2=outbam_h2, output_untagged=outbam_untagged)

		expected_splits = {}
		for alignment in pysam.AlignmentFile(outbam1):
			if alignment.has_tag('HP'):
				expected_splits[alignment.query_name] = alignment.get_tag('HP')
			else:
				expected_splits[alignment.query_name] = 0
		assert len(expected_splits) > 0

		total_reads = 0
		for hap, outfile in enumerate([outbam_untagged, outbam_h1, outbam_h2]):
			for alignment in pysam.AlignmentFile(outfile):
				hap_key = hap
				assert expected_splits[alignment.query_name] == hap_key
				total_reads += 1
		assert total_reads == len(expected_splits)


def test_split_poly():
	with TemporaryDirectory() as tempdir:
		outlist = tempdir + '/outlist.tsv'
		outbam = tempdir + '/outbam.bam'
		# produce list of read assignments to haplotypes
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

	

