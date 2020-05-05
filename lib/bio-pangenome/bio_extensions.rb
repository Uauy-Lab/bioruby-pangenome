require 'bio'
module Bio::Extensions
	def self.fasta_gz_iterator(path)
		infile = open(path)
		io = Zlib::GzipReader.new(infile) 
		Bio::FlatFile.open(Bio::FastaFormat, io) do |fasta_file|
			fasta_file.each do |entry|
				yield entry
			end
		end
		io.close
	end
end