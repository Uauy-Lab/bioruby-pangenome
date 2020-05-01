 class Bio::Alignment::SequenceHash
 	def to_fasta( size: 60)
 		self.keys.map { |name, val|  self[name].to_fasta(name, size) }.join("")
 	end

 	def masks
 		
 	end

 end

 class Bio::PanGenome::MaskSite
 	attr_accessor :gap, :ambiguous, :snp 
 end