#require 'helper'
require "test/unit"
gempath = File.dirname(File.dirname(__FILE__))
path = gempath + '/lib/bio-pangenome.rb'
require path

 
class Test::Pangenome < Test::Unit::TestCase
	def setup
		@genes = Bio::Pangenome.load_genes("./test/data/genes.txt")
		@lines = Bio::Pangenome.load_lines("./test/data/lines.txt")
		@projected_genes = Bio::Pangenome.load_projected_genes("./test/data/mapping.csv.gz", genes: @genes)
	end

	def teardown
		## Nothing really
	end

	def test_load_genes
		genes = Bio::Pangenome.load_genes("./test/data/genes.txt")
		assert_equal(50, genes.size)
		genes = Bio::Pangenome.load_genes("./test/data/genes.txt", window:0, no_windows: 2)
		assert_equal(25, genes.size)
		genes = Bio::Pangenome.load_genes("./test/data/genes.txt", window:1, no_windows: 2)
		assert_equal(25, genes.size)
	 end

	def test_load_lines
		lines = Bio::Pangenome.load_lines("./test/data/lines.txt")
		assert(15, lines.inspect)
	end

	def test_load_sequences
		sequences = Bio::Pangenome.load_sequences(varieties: @lines, genes: @projected_genes, prefix: "./test/data/2000bp", suffix: "_2000bp_RefSeqv1.1.fa.gz", set_id: "2000bp")

		assert_equal( 11, sequences["TraesCS6A02G002800"].length )
		assert_equal( 13, sequences["TraesCS6A02G006600"].length )
		assert_equal(  5, sequences["TraesCS6A02G008400"].length )
		assert_equal( 12, sequences["TraesCS6A02G017400"].length )
		assert_equal( 15, sequences["TraesCS6A02G027000"].length )
		assert_equal( 11, sequences["TraesCS6A02G032300"].length )
		assert_equal( 14, sequences["TraesCS6A02G052900"].length )
		assert_equal( 15, sequences["TraesCS6A02G077500"].length )
		assert_equal(  9, sequences["TraesCS6A02G081800"].length )
		assert_equal( 15, sequences["TraesCS6A02G084500"].length )
		assert_equal( 12, sequences["TraesCS6A02G096900"].length )
		assert_equal( 14, sequences["TraesCS6A02G097400"].length )
		assert_equal( 12, sequences["TraesCS6A02G103000"].length )
		assert_equal( 11, sequences["TraesCS6A02G120500"].length )
		assert_equal( 15, sequences["TraesCS6A02G134500"].length )
		assert_equal( 15, sequences["TraesCS6A02G137400"].length )
		assert_equal( 15, sequences["TraesCS6A02G143300"].length )
		assert_equal( 10, sequences["TraesCS6A02G163100"].length )
		assert_equal( 13, sequences["TraesCS6A02G177100"].length )
		assert_equal( 15, sequences["TraesCS6A02G177700"].length )
		assert_equal( 14, sequences["TraesCS6A02G183800"].length )
		assert_equal( 12, sequences["TraesCS6A02G197000"].length )
		assert_equal( 15, sequences["TraesCS6A02G212100"].length )
		assert_equal( 15, sequences["TraesCS6A02G224000"].length )
		assert_equal( 15, sequences["TraesCS6A02G235700"].length )
		assert_equal( 15, sequences["TraesCS6A02G239100"].length )
		assert_equal( 10, sequences["TraesCS6A02G243900"].length )
		assert_equal( 15, sequences["TraesCS6A02G256700"].length )
		assert_equal( 15, sequences["TraesCS6A02G261200"].length )
		assert_equal( 11, sequences["TraesCS6A02G264300"].length )
		assert_equal( 12, sequences["TraesCS6A02G418100"].length )
		assert_equal( 10, sequences["TraesCS6A02G282100"].length )
		assert_equal( 15, sequences["TraesCS6A02G295200"].length )
		assert_equal( 15, sequences["TraesCS6A02G307900"].length )
		assert_equal( 11, sequences["TraesCS6A02G310500"].length )
		assert_equal( 14, sequences["TraesCS6A02G312900"].length )
		assert_equal( 13, sequences["TraesCS6A02G329200"].length )
		assert_equal( 15, sequences["TraesCS6A02G329600"].length )
		assert_equal( 9, sequences["TraesCS6A02G339000"].length )
		assert_equal( 15, sequences["TraesCS6A02G389900"].length )
		assert_equal( 15, sequences["TraesCS6A02G399100"].length )
		assert_equal( 11, sequences["TraesCS6A02G403300"].length )
		assert_equal(  10, sequences["TraesCS6A02G012900"].length )
		assert_equal(  8, sequences["TraesCS6A02G087800"].length )
		assert_equal(  7, sequences["TraesCS6A02G020700"].length )
		assert_equal(  8, sequences["TraesCS6A02G025000"].length )
		assert_equal(  4, sequences["TraesCS6A02G090800"].length )
		assert_equal(  5, sequences["TraesCS6A02G206300"].length )
		assert_equal(  6, sequences["TraesCS6A02G378300"].length )
		assert(sequences["TraesCS6A02G008400"]["arinalrfor"].sequence.start_with? "AACAGTAAATCCAAAAATTAGAAAAGAAATTCAAAAGAATCTGGAATTTTTGGACCAGACCTATCCGGTTTTGAGACTTGAAGGACCAAATGTCTACCAAAAAATCTTGAGGGACCAAGTCTGTAATTTTGGGAACTTGAGAGACCAAAAACATAAATTTCTCTTATAAGATC", "Sequence dont start as expected")

	end
end
