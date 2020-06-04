#!/usr/bin/env ruby
#
# BioRuby bio-pangenome Plugin Bio::Pangenome
# Author:: Ricardo H Ramirez Gonzalez
# Copyright:: 2019,2020

USAGE = "pangenome_gene_bed_files.rb [options]"

gempath = File.dirname(File.dirname(__FILE__))
$: << File.join(gempath,'lib')

VERSION_FILENAME=File.join(gempath,'VERSION')
version = File.new(VERSION_FILENAME).read.chomp

# print banner
print "pangenome_gene_bed_files #{version} by Ricardo H. Ramirez-Gonzalez 2020\n"

if ARGV.size == 0
  print USAGE
  exit
end

path = gempath + '/lib/bio-pangenome.rb'
require path
#require 'bio-pangenome'
require 'optparse'
require 'tmpdir'

require 'csv'
require 'zlib'
require 'bio'
require 'bio-svgenes'


options = {
    lines:     "lines.txt",
    distances: [0,1000,2000,5000],
    path: "./gff/",
    is_gz: true, 
    gff_suffix: "gff.gz"
}

opts = OptionParser.new do |o|

    o.on("-l", "--lines PATH", "File containing the lines to be analysed") do |arg|
        options[:lines] = arg
    end

    o.on("-d", "--distances 0,1000,2000", "File containing the distances to be analysed") do |arg|
        options[:distances] = arg.split(",").map { |e|  e.to_i }
    end

    o.on("-p", "--gff_path DIR", "The directory where the gff files are") do |arg|
        options[:path ] = arg
    end

    o.on("-s", "--gff-suffix gff.gz", "The extension of the gff.") do |arg|
        options[:gff_suffix] = arg 
    end

    o.separator ""
    o.on_tail('-h', '--help', 'display this help and exit') do
        options[:show_help] = true
    end

end

opts.parse!(ARGV)

puts options.inspect

lines = File.foreach(options[:lines]).map { |line| line.chomp }
distances = options[:distances]
puts lines



gffs = MultipleGFFs.new(folder: options[:path], lines:lines, suffix:  options[:gff_suffix], 
    is_gz:options[:is_gz] )

distances.each do |d|
    gffs.bedAround(distance: d, prefix: options[:path], suffix: ".bed" )
end


 def bedAroundToRegions(lines:[], distance: 1000, prefix: "../flanking/filtered/", suffix: ".RefSeqv1.1.reg" , suffix_in: ".RefSeqv1.1.bed" )
     lines.each do |k|
         path="#{prefix}#{k}_#{distance}bp_#{suffix_in}"
         path2="#{prefix}#{k}_#{distance}bp_#{suffix}"
         path3="#{prefix}#{k}_#{distance}bp_#{suffix}.map"
         puts path
         out=File.open(path2, "w")
         out2=File.open(path3, "w")
         File.foreach(path) do |line|
            # puts line
             arr = line.chomp!.split "\t"
             first=arr[1]
             last=arr[2]
             name=arr[0]
             #if(arr[5] == "-")
             #    first=arr[2]
             #    last=arr[1]
             #end
           
             reg =  "#{name}:#{first}-#{last}"
             out.puts reg
             out2.puts [reg,arr[3]].join "\t"
             #puts reg
             #v.bedAroundGene(distance:distance, out:out)
             #break
         end
         out.close
         out2.close
     end
 end



 distances.each do |d|
     bedAroundToRegions(lines:lines,
         distance: d,
         prefix: options[:path],
         suffix_in: ".bed", 
         suffix: ".reg")
 end






