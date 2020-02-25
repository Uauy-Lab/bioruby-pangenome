#!/usr/bin/env ruby
#
# BioRuby bio-pangenome Plugin BioPangenome
# Author:: homonecloco
# Copyright:: 2019

USAGE = "pangenome_gene_bed_files.rb [options]"

gempath = File.dirname(File.dirname(__FILE__))
$: << File.join(gempath,'lib')

VERSION_FILENAME=File.join(gempath,'VERSION')
version = File.new(VERSION_FILENAME).read.chomp

# print banner
print "pangenome_gene_bed_files #{version} by Ricardo H. Ramirez-Gonzalez 2020\n"

if ARGV.size == 0
  print USAGE
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


lines=["cadenza", "claire", "paragon", "robigus","weebil"]
distances=[0,1000,2000,5000]
gffs = MultipleGFFs.new(folder: "/Volumes/PanGenome/References/releaseV3/gff", lines:lines, suffix:".gff.gz", 
    is_gz:true )

distances.each do |d|
    gffs.bedAround(distance: d, prefix: "/Volumes/PanGenome/GeneRegions/201910_v2_v3/", suffix: ".RefSeqv1.1.bed" )
end

lines2=["arinalrfor", "ashsyn", "chinese", "jagger", "julius", "lancer", "landmark", "mace", "norin61", "spelta", "stanley", "sy_mattis"]
distances=[0,1000,2000,5000]
gffs = MultipleGFFs.new(folder: "/Volumes/PanGenome/References/releasePGSBv2.0/gff", lines:lines2, suffix:".gff.gz", 
    is_gz:true )

distances.each do |d|
    gffs.bedAround(distance: d, prefix: "/Volumes/PanGenome/GeneRegions/201910_v2_v3/", suffix: ".RefSeqv1.1.bed" )
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
    bedAroundToRegions(lines:lines+lines2,
        distance: d,
        prefix: "/Volumes/PanGenome/GeneRegions/201910_v2_v3/",
        suffix_in: ".RefSeqv1.1.bed", 
        suffix: ".RefSeqv1.1.reg")
end






