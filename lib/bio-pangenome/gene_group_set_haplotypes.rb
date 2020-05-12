module Bio::Pangenome
	class GeneGroupSetHaplotypes 
		attr_reader :varieties, :gene_order, :chunks, :chunk_no, :blocks
		def initialize(chunks: 1, chunk_no: 0)
			@varieties = Set.new 
			@gene_order = []
			@chunks = chunks
			@chunk_no = chunk_no
			@blocks = Hash.new do|h,k| 
				tmp = GeneGroupSet.new { |h, k| h[k] = GeneGroup.new }
				tmp.variety_haplotype = Hash.new  
				tmp.genes     = []
				tmp.varieties = Set.new
				h[k]=tmp
			end
		end

		def groupset_for_gene
			return @groupset_for_gene if @groupset_for_gene
			blocks.freeze
			@groupset_for_gene = Hash.new
			@blocks.each_pair do |k,v|
				v.genes.each do |gene|
					@groupset_for_gene[gene] = k
				end
			end
			return @groupset_for_gene
		end

		def load_sequences(genes:{}, prefix: "../flanking/filtered/",  suffix: ".cds.fa.gz", set_id: "cds" )
	
			varieties.each do |variety|
				path = "#{prefix}/#{variety}#{suffix}"
				Bio::Extensions.fasta_gz_iterator(path) do |entry|
					seq_name, definition, region, arr = Bio::Pangenome::parse_entry_name(entry)
					next unless genes.include? seq_name.gene
					row = genes[seq_name.gene]
					row = genes[arr[0]] unless region
					seq_name.gene = row["gene"]
					seq = entry.seq
					seq_name.sequence = seq
					base_gene = row["gene"]
					block = groupset_for_gene[base_gene]
					@blocks[block][base_gene][variety] = seq_name  unless @blocks[block][base_gene][variety]
					@blocks[block][base_gene].gene = base_gene
					@blocks[block][base_gene].expected_varieties = varieties
				end
			end
		end
	end


	def self.parse_entry_name(entry)
		definition, region = entry.definition.split("::")
		seq_name = parseSequenceName(region, definition) if region
		seq_name = GeneFlankingRegion.new(entry.definition, nil, "",
			"", entry.definition, set_id, nil, variety ) unless region
		arr = definition.split(".")

		[seq_name, definition, region, arr]
	end

	def self.load_haplotye_groups(genes_file, haplotypes_file, chunk_no: 0, chunks: 1)
		ret = GeneGroupSetHaplotypes.new(chunk_no: chunk_no, chunks: chunks)
		CSV.foreach(haplotypes_file, col_sep:"\t", headers: true) do |row|
			ret.varieties << row["line"]
			ret.blocks[row["block"]].varieties << row["line"]
			ret.blocks[row["block"]].variety_haplotype[row["line"]] = row["haplotype"]
		end

		ret.blocks.keys.each_with_index do |e, i|
			ret.blocks.delete(e) unless i % chunks == chunk_no
		end

		CSV.foreach(genes_file, col_sep:"\t", headers: true) do |row|
			next unless ret.blocks.keys.include? row["block"]
			ret.gene_order << row["gene"]
			ret.blocks[row["block"]].genes << row["gene"]
		end
		ret
	end
end