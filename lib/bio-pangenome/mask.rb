require 'bio'


module Bio::Pangenome::Mask 
  def output_fasta(size: 60)
    str = alignment.keys.map { |name|  ">#{name}\n#{alignment[name]}"}.join("\n")
    str += "\n>Mask\n#{mask}"
    str
  end

end

class Bio::Pangenome::HaplotypeMask
  attr_reader :alignment ,:match_line_char, :gap_char, :mismatch_char, :missing_char, :flanking_bases

  include Bio::Pangenome::Mask
  def initialize( alignment, match_line_char: ".", gap_char: "-", mismatch_char: "|", flanking_bases: 50, missing_char: "N")
    @match_line_char = match_line_char.freeze
    @gap_char = gap_char.freeze
    @mismatch_char = mismatch_char.freeze
    @alignment = alignment.freeze
    @missing_char = missing_char.freeze
    @flanking_bases = flanking_bases.freeze
  end

  def site_mask(site)
    a = site.collect { |c| c.upcase }.sort.uniq
    chr = mismatch_char
    chr = gap_char if site.has_gap?
    chr = match_line_char if a.size == 1 
    chr = missing_char if a.include? "N"
    chr
  end

  def mask
    return @mask if @mask
    str = "." * alignment.alignment_length
    i = 0
    alignment.each_site  do |s|
      str[i] = site_mask(s)
      i += 1
    end 
    @mask = str.freeze
  end

  def stats
    return @stats if @stats
    @stats = Hash.new
    @stats[:gaps] = mask.count(gap_char)
    @stats[:mismatches] = mask.count(mismatch_char) 
    @stats[:matches] = mask.count(match_line_char)
    @stats[:missing] = mask.count(missing_char)
    @stats[:length] = mask.length
    @stats[:mismatches_per_kbp] = 1000.0 * @stats[:mismatches] / ( @stats[:length]  - @stats[:missing] - @stats[:gaps] )
    @stats
  end

  def snp_positions
    return @snp_positions if @snp_positions 
    @snp_positions = []
    mask.each_char.each_with_index do |c, i|  
      @snp_positions << i if c == mismatch_char
    end
    @snp_positions.freeze
  end
 
  def valid_snps 
    return @valid_snps if @valid_snps
    @valid_snps = []
    snp_positions.each do |pos|
      next if pos < flanking_bases
      start = pos - flanking_bases
      tmp_mask = mask[start , flanking_bases*2 + 1]
      dots = tmp_mask.count(match_line_char)
      @valid_snps << pos if dots == flanking_bases * 2
    end
    @valid_snps
  end

end

