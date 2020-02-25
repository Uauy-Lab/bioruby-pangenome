require 'bio-gff3'

module Bio::GFFbrowser::FastLineParser
	module_function :parse_line_fast
end

MrnaStats = Struct.new(:cds_count, :cds_max_gap)

class GFF3
	CDS_feature = Struct.new(:start, :end, :color, :ref_chr,:ref_start, :ref_end, :offset_start)

	def initialize(file: "", is_gz: true)
		@file = file
		@is_gz = is_gz
	end

	def each
		return enum_for(:each) unless block_given? 
		io = nil
		if @is_gz
			infile = open(@file)
			io = Zlib::GzipReader.new(infile) 
		else
			io =  File.open(@file)
		end
		parser = Bio::GFFbrowser::FastLineParser
		io.each_line do |line|  
			line.encode!('UTF-8', 'UTF-8', :invalid => :replace)
			line.strip!
			break if line == '##FASTA'
			next if line.length == 0 or line =~ /^#/
			begin
				record = Bio::GFFbrowser::FastLineRecord.new(parser.parse_line_fast(line))
				yield record
			rescue Exception => e
				$stderr.puts "Unable to parse '#{line}'\n#{e}" 
				throw e
			end
		end
	end

	def each_gene
		return enum_for(:each_gene) unless block_given? 
		self.each do |record|
			next unless record.feature == "gene"
			yield record
		end 
	end

	def each_mrna
		return enum_for(:each_mrna) unless block_given? 
		self.each do |record|
			next unless record.feature == "mRNA"
			yield record
		end 
	end

	def each_cds
		return enum_for(:each_mrna) unless block_given? 
		self.each do |record|
			next unless record.feature == "CDS"
			yield record
		end 
	end

	def calculate_mrna_stats
		return if @mrna_stats
		@mrna_stats = Hash.new {|h,k| h[k] = MrnaStats.new(0,0) }
		last_mrna = ""
		last_record = nil
		each_cds do |record|
			parent = record.get_attribute "Parent"
			mrna = @mrna_stats[parent]
			mrna.cds_count += 1
			if last_mrna == parent
				distance =  record.start - last_record.end 
				mrna.cds_max_gap = distance if distance > mrna.cds_max_gap
			end
			last_record = record
			last_mrna   = parent
		end
		return
	end

	def mrna_info(id)
		calculate_mrna_stats
		@mrna_stats[id] 
	end

	def bedAroundGene(distance:1000, out:$stdout)
		each_mrna do |record|
			start = record.start-distance
			start = 1 if start < 1
			reg_end=record.end + distance
			out.puts [record.seqid, start, reg_end, "#{record.id}_#{record.source}_#{distance}bp", ".", record.strand].join "\t"
		end
	end


	def cds_to_print(mrna,cannonical_exons:[], colors:["#a6cee3", "#1f78b4", "#b2df8a" , "#33a02c", "#fb9a99",  "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"])

		cds_features = [] 
		i = 0
		offset=0
		offset_start=0
		each_cds do |record|
			target = record.get_attribute "Target"
			arr = target.split(" ")
			col = colors[i % colors.size ]
			start = arr[1].to_i + offset
			ends = arr[2].to_i + offset
			offset_start = record.start  if offset_start == 0
			tmp = CDS_feature.new(start, ends, col, 
				record.seqid, record.start,record.end, record.start - offset_start )
			cds_features << tmp
			i += 1
		end
		cds_features
	end

end
