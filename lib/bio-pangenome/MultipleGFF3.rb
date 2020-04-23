class MultipleGFFs
    attr_reader :lines_gffs

    def initialize(folder: "../mapping/", lines:[], suffix:".SM1.cds.sorted.gff", is_gz:false )
        @folder = folder
        @lines = lines
        @suffix = suffix
        @lines_gffs = Hash.new
        @lines.each do |l|
            path ="#{folder}/#{l}#{suffix}"
            @lines_gffs[l] = GFF3.new(file: path, is_gz: is_gz)
        end
    end
    
    def each_gff
        @lines_gffs.each_pair{|k,v| yield k, v }
    end

    def bedAround(distance: 1000, prefix: "../flanking/releasePGSBV1/", suffix: ".RefSeqv1.1.bed" )
        each_gff do |k, v|
            path="#{prefix}#{k}_#{distance}bp_#{suffix}"
            puts path
            out=File.open(path, "w")
            v.bedAroundGene(distance:distance, out:out)
            out.close
        end
    end

    def summary
        ret = []
        each_gff do |k,v|
            v.each_mrna do |record|
                tmp = {}
                tmp[:line] = k
                tmp[:id] = record.get_attribute "Name"
                tmp[:chr] = record.seqid
                tmp[:start] = record.start
                tmp[:end] = record.end
                tmp[:strand] = record.strand
                tmp[:genomic_length] = record.end - record.start
                tmp[:coverage] = record.get_attribute "coverage"
                tmp[:identity] = record.get_attribute "identity"
                tmp[:matches]  = record.get_attribute "matches"
                tmp[:mismatches]  = record.get_attribute "mismatches"
                tmp[:indels] = record.get_attribute "indels"
                tmp[:unknowns] = record.get_attribute "unknowns"
                mrna_stats = @lines_gffs[k].mrna_info(record.id)
                tmp[:cds_count]   = mrna_stats.cds_count
                tmp[:cds_max_gap] = mrna_stats.cds_max_gap
                ret << tmp
            end
        end
        ret
    end

    def to_svg(mrna: "Sm1_CDS.mrna1", positions: false, out: nil)
        p = Bio::Graphics::Page.new(width: 800,
         height: 1000, 
         number_of_intervals:10,
         background_color: "white"
         )
        each_gff do |k,v|
            generic_track = p.add_track(:glyph => :generic, 
                :name => k, 
                :label => true  )
            v.cds_to_print(mrna).each do |cds|

                f_id = positions ? cds.offset_start : nil
                feature = Bio::Graphics::MiniFeature.new(start: cds.start, 
            end: cds.end,
            fill_color: cds.color, 
            id: f_id)
                generic_track.add(feature) 
            end 
        end
    end
end
