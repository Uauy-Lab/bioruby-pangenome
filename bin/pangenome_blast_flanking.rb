#!/usr/bin/env ruby
#
# BioRuby bio-pangenome Plugin BioPangenome
# Author:: homonecloco
# Copyright:: 2019

USAGE = "panggenome_blast_flanking.rb [options]"

gempath = File.dirname(File.dirname(__FILE__))
$: << File.join(gempath,'lib')

VERSION_FILENAME=File.join(gempath,'VERSION')
version = File.new(VERSION_FILENAME).read.chomp

# print banner
print "panggenome_blast_flanking #{version} by Ricardo H. Ramirez-Gonzalez 2019\n"

if ARGV.size == 0
  print USAGE
end

path = gempath + '/lib/bio-pangenome.rb'
require path
#require 'bio-pangenome'
require 'optparse'
require 'tmpdir'

options = { 
  :transcript_mapping => "sorted_filtered_mapping.csv.gz",
  :lines => "lines.txt",
  :genes => "genes.txt",
  :no_windows => 0, 
  :window => 0,
  :distance => 2000
}
opts = OptionParser.new do |o|
  o.banner = "Usage: #{File.basename($0)} [options]"

  o.on('-t', '--transcript_mapping [sorted_filtered_mapping.csv.gz]', 'File with the  mappings across  transcriptomes') do |arg|
    options[:transcript_mapping] = arg
  end

  o.on('-g','--genes [genes.txt]', 'File with the list of genes') do |arg|
    options[:genes] = arg
  end

  o.on('-n','--no_windows INT', "Number of chunks to divide the genes list. 0 to not split") do |arg|
    options[:no_windows] = arg.to_i
  end

  o.on('-w', "--window INT", "Current window to run") do |arg|
    options[:window] = arg.to_i 
  end

  o.on('-d', "--distance DIST", "Name of the distance set. Ues 'cds' to align cds. default 2000") do |arg|
    options[:distance] = arg
  end

  o.on("-b", "--basepath PATH", "Folder with the sequences and mapping positions across genomes") do |arg|
    options[:basepath] = arg
  end

  o.on("-o", "--output PATH", "Folder with the output. If  there are chunks, they will be used") do |arg|
    options[:output] = arg
  end

  o.on("-l", "--lines PATH", "File containing the lines to be analysed") do |arg|
    options[:lines] = arg
  end

  o.separator ""
  o.on_tail('-h', '--help', 'display this help and exit') do
    options[:show_help] = true
  end
end

opts.parse!(ARGV)

genes = BioPangenome.load_genes(options[:genes], window: options[:window], no_windows: options[:no_windows] )
puts "Genes count: #{genes.size}"

lines = BioPangenome.load_lines(options[:lines])

projected_genes = BioPangenome.load_projected_genes options[:transcript_mapping], genes: genes

variety_coordinates = BioPangenome.load_mapping_hash(
  varieties:lines, 
  genes: projected_genes ,
  prefix: options[:basepath],
  distance: options[:distance]
  )

seqs = BioPangenome.load_sequences_from_hash(
  coordinates:variety_coordinates,
  prefix: options[:basepath],
  distance: options[:distance],
  projected_genes: projected_genes
  )

output = options[:output].to_s
output = output + "_" + options[:window].to_s if options[:no_windows] > 0

Dir.mktmpdir do |temp_dir|
  puts "Aligning in  #{temp_dir}"
  BioPangenome.align_gene_groups(
    seqs: seqs, 
    distance:options[:distance],
    output: output,
    tmp_folder: temp_dir)
end




