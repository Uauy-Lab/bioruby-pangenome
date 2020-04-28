#!/usr/bin/env ruby
#
# BioRuby bio-pangenome Plugin BioPangenome
# Author:: Rciardo H. Ramirez Gonzalez
# Copyright:: 2020

#!/usr/bin/env ruby
#
# BioRuby bio-pangenome Plugin BioPangenome
# Author:: Rciardo H. Ramirez Gonzalez
# Copyright:: 2019, 2020

USAGE = " #{File.basename($0)} [options]"

gempath = File.dirname(File.dirname(__FILE__))
$: << File.join(gempath,'lib')

VERSION_FILENAME=File.join(gempath,'VERSION')
version = File.new(VERSION_FILENAME).read.chomp

# print banner
print "panggenome_blast_flanking #{version} by Ricardo H. Ramirez-Gonzalez 2020\n"

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

  o.on("-p", "--gff_path DIR", "The directory where the gff files are") do |arg|
  	options[:path ] = arg
  end

  o.separator ""
  o.on_tail('-h', '--help', 'display this help and exit') do
    options[:show_help] = true
    puts o
    exit
  end
end