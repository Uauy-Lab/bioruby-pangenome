
require 'zlib'
require 'bio'
require 'csv'
require 'set'
require 'bio-blastxmlparser'

module Bio



module Pangenome
	Transcript = Struct.new(:id, :gene, :chromosome,:version,:count,:transcript,:confidence, :count_int, :isoform)
	GeneFlankingRegion = Struct.new(:transcript, :gene, :ann, :region, :id, :flank_length, :sequence, :line)

	def self.parseTranscript name
		arr=name.split(".")
		match = /TraesCS(?<chr>[[:alnum:]]{1,2})(?<ver>[[:digit:]]{2})G(?<count>[[:digit:]]+)(?<conf>[[:upper:]]*)/.match arr[0]
		raise "Unable to parse: #{name}" unless match
		Transcript.new(name, arr[0],match[:chr],match[:ver],match[:count],arr[1],match[:conf], match[:count].to_i, arr[1])
	end
	def self.parseEITranscript name
		arr=name.split(".")
		match = /Traes(?<chr>[[:upper:]]{3}_scaffold_[[:digit:]]*)_(?<ver>[[:digit:]]{2})G(?<count>[[:digit:]]+)(?<conf>[[:upper:]]*)/.match arr[0]
		raise "Unable to parse: #{name}" unless match
		Transcript.new(name, arr[0],match[:chr].downcase,match[:ver],match[:count],arr[1],match[:conf], match[:count].to_i, arr[1])
	end

	def self.parsePGSBTranscript name
		arr=name.split(".")
		match = /Traes(?<variety>[[:upper:]]{3})(?<chr>[[:alnum:]]{1,2})(?<ver>[[:digit:]]{2})G(?<count>[[:digit:]]+)(?<conf>[[:upper:]]*)/.match arr[0]

		raise "Unable to parse: #{name}" unless match
		Transcript.new(name, arr[0],match[:chr],match[:ver],match[:count],match[:variety], match[:conf],match[:count].to_i, arr[1])
	end
	def self.parseSequenceName region, name
		match = /(?<transcript>[[:alnum:]].+)_(?<ann>.+)_(?<flank_length>[[:digit:]]+bp)/.match name
		arr2=match[:transcript].split "."
		GeneFlankingRegion.new(match[:transcript],arr2[0],match[:ann], region, name, match[:flank_length] , nil, nil)
	end

	def self.load_mapping_hash(varieties:[],  transcripts:[], genes:[], distance: 1000, prefix: "../flanking/releasePGSBV1/",  suffix: ".RefSeqv1.1")
		ret = Hash.new { |h, k| h[k] = Hash.new }
		varieties.each do |v|
			path = "#{prefix}#{distance}bp/#{v}_#{distance}bp_#{suffix}.reg.map"
			$stderr.puts path
			File.foreach(path) do |line|
				line.chomp!
				arr = line.split("\t")
				begin
					parsed = parseSequenceName(arr[0], arr[1])
				rescue Exception => e
					throw "Unable to parse #{line} (#{v}) [#{e.to_s}]" 
				end
				next unless transcripts.include? parsed.transcript or genes.include? parsed.gene
				ret[v][parsed.region] = parsed
			end
		end
		ret
	end

	def self.blast_pair_fast(path_a, path_b, out_path, program: "blastn")
		cmd = "#{program} -query #{path_a} -subject #{path_b} -task #{program} -out #{out_path} -outfmt '5' "
		system cmd
		n = Bio::BlastXMLParser::XmlIterator.new(out_path).to_enum
		max_length = 0
		max_pident = 0.0
		n.each do | iter |
			iter.each do | hit |
				hit.each do | hsp |
					if hsp.align_len > max_length
						max_length = hsp.align_len
						max_pident = 100 * hsp.identity.to_f / hsp.align_len.to_f
					end
				end
			end
		end
		[max_length, max_pident]
	end



	def self.load_sequences_from_hash(coordinates:{},  prefix: "../flanking/filtered/",  suffix: "RefSeqv1.1", distance: 1000, projected_genes: {})
		ret = Hash.new { |h, k| h[k] = Hash.new }
		coordinates.each_pair do |variety, coords|

			path = "#{prefix}/#{distance}bp/#{variety}_#{distance}bp_#{suffix}.fa.gz"
			puts "Loading: #{path}"
			infile = open(path)
			io = Zlib::GzipReader.new(infile) 
			Bio::FlatFile.open(Bio::FastaFormat, io) do |fasta_file|
				fasta_file.each do |entry|
					next unless coords[entry.definition]
					seq_name = coords[entry.definition]
					seq = entry.seq
					seq.gsub!(/N*$/, '')
					seq.gsub!(/^N*/, '')
					seq_name.sequence = seq
					base_gene = projected_genes[seq_name.gene]["gene"]
					ret[base_gene][variety] = seq_name unless ret[base_gene][variety]
				end
			end
			io.close
		end
		ret
	end

	def self.align_gene_groups( seqs:{}, tmp_folder:"/Volumes/PanGenome/GeneRegions/201910_v2_v3/tmp", output:"../pairwise_blast_oct_2019/varieties_6A_identites", distance: 0 )
		out_tmp="#{tmp_folder}/out.blast"
		FileUtils.mkdir_p(tmp_folder)
		out = File.open("#{output}_#{distance}bp.tab", "w")
		out.puts [ "transcript" , "query", "subject" ,  "var_query", "var_subject", "aln_type", "length" , "pident" , "Ns_query", "Ns_subject", "Ns_total", "Flanking"   ].join("\t")
		seqs.each_pair do |transcript, transcript_seqs|
			vars = transcript_seqs.keys
			vars_done = []
			ns = {}
			vars.each do |v1|
				tmp =  tmp_folder  + "/" + v1 + ".fa"
				s = transcript_seqs[v1]
				seq = ">#{s.id}\n#{s.sequence}"
				File.open(tmp, 'w') {|f| f.write(seq) }
				ns[v1] = s.sequence.count('Nn')
			end
			vars.each do |v1|
				tmp1 =  tmp_folder  + "/" + v1 + ".fa"
				s1 = transcript_seqs[v1]
				next unless s1.sequence.length > 0
				vars.each do |v2|
					next if v1 == v2
					next if vars_done.include? v2
					s2 = transcript_seqs[v2]
					next unless s2.sequence.length > 0
					tmp2 =  tmp_folder  + "/" + v2 + ".fa"
					to_print = [transcript, s1.id , s2.id , v1,v2,"#{v1}->#{v2}"]
					to_print << blast_pair_fast(tmp1, tmp2, out_tmp) 
					to_print << ns[v1] 
					to_print << ns[v2]
					to_print << ns[v1] + ns[v2]
					to_print << distance
					out.puts to_print.join("\t")
				end
				vars_done << v1
			end
		end
		out.close
	end

	def self.load_cds_sequencess( varieties:[], genes:{}, prefix: "../flanking/filtered/",  suffix: ".cds.fa.gz", set_id: "cds" )
		return load_sequences(varieties: varieties, genes: genes, prefix: prefix, suffix: suffix, set_id: set_id)
	end

	def self.load_sequences( varieties:[], genes:{}, prefix: "../flanking/filtered/",  suffix: ".cds.fa.gz", set_id: "cds" )
		return GeneGroup.load_sequences(varieties: varieties, genes: genes, prefix: prefix, suffix: suffix, set_id: set_id)
	end

	def self.load_projected_genes(transcript_mapping, genes:[])
		projected_genes = {}
		Zlib::GzipReader.open(transcript_mapping) do |gzip|
			csv = CSV.new(gzip, headers: true)
			csv.each do |row|
				next unless genes.include? row["gene"]
				projected_genes[row["projected_gene"]] = row
			end
		end
		projected_genes
	end

	def self.load_genes(filename, window: 0, no_windows: 0)
		genes = File.readlines(filename).map do |t|
			next if t.nil?
			t.chomp!.split(".")[0]
		end
		if no_windows > 0
			puts "'loading window #{window} of #{no_windows}'"
			window_size = genes.size/no_windows
			start = window * window_size
			genes = genes[start, window_size]
		end
		genes
	end

	def self.load_lines(filename)
		File.readlines(filename).map do |t|
			t.chomp!.rstrip
		end
	end
end
end
