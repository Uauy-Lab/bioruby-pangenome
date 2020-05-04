require 'bio'

class Array
  def sum
    inject(0.0) { |result, el| result + el }
  end

  def mean
    sum / size
  end
end


module Bio::Pangenome::Mask 
  def output_fasta(size: 60)
    str = alignment.keys.map { |name|  ">#{name}\n#{alignment[name]}"}.join("\n")
    str += "\n>Mask\n#{mask}"
    str
  end

end

class Bio::Pangenome::HaplotypeMask
  attr_accessor :alignment, :match_line_char, :gap_char, :mismatch_char
  include Bio::Pangenome::Mask
  def initialize( alignment, match_line_char: ".", gap_char: "-", mismatch_char: "|")
    self.match_line_char = match_line_char
    self.gap_char = gap_char
    self.mismatch_char = mismatch_char
    self.alignment = alignment
  end


  def mask
    return @mask if @mask
    str = "." * alignment.alignment_length
    i = 0
    alignment.each_site  do |s|
      a = s.collect { |c| c.upcase }.sort.uniq
      chr = mismatch_char
      chr = gap_char if s.has_gap?
      chr = match_line_char if a.size == 1 
      str[i] = chr
      i += 1
    end 
    @mask = str.freeze
  end

  def stats
    return @stats if @stats
    @stats = Hash.new
    @stats[:gaps] = mask.count(gap_char)
    @stats[:mismatchs] = mask.count(mismatch_char) 
    @stats[:matches] = mask.count(match_line_char)
    @stats[:length] = mask.length
    @stats
  end
 

end

