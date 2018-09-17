#
#	String class mods for parsing sequences
#
class String

	def indices_of_chars(chars)
		start, indices = -1, []
		indices << start while start = (self.index /[#{chars}]/, start + 1)
		indices
	end

	def indices_of_mods()
		tmp=self
		indices = []
		while (i = tmp.index /\(/)
			indices << i - 1
			tmp.sub! /\(.*?\)/,''
		end
		indices
	end

end


class Evidence

	attr_accessor :sequence, :protein, :cleaned_sequence, :alignments, :modified_positions,
		:intensity, :experiment, :msmscount, :matched_amino_acids

	@@dir = "PeptideMatchOutput"

	def initialize(h={})
		@sequence = h['Modified sequence']     #	1
		@protein = h['Leading razor protein']  #	15
		@intensity = h['Intensity']            #	54
		@experiment = h['Experiment']          #	19
		@msmscount = h['MS/MS Count']          #	49
		@cleaned_sequence = sequence.gsub(/_/,'').gsub(/\([[:alpha:]]*\)/,'')
		@alignments = []
		@modified_positions = sequence.gsub(/_/,'').indices_of_mods
		@matched_amino_acids = false
	end

	def query( database )
		Dir.mkdir(@@dir) unless Dir.exists?(@@dir)
		outfile="#{@@dir}/#{database}.#{cleaned_sequence}.txt"

		#	puts outfile

		if !File.exists?( outfile )
			puts "Querying."
			system("java -jar PeptideMatchCMD_1.0.jar -a query -i #{database} -q #{cleaned_sequence} -o #{outfile}")
		else
			puts "Using existing query."
		end

		File.open(outfile,'rb').each do |line|
			next if line =~ /^#/
			l = line.split
			if l[0] == cleaned_sequence && l[1] == protein
				##Query	Subject	SubjectLength	MatchStart	MatchEnd
				alignments.push({ match_start: l[3].to_i, match_end: l[4].to_i,
					absolute_positions: {}, modified_positions: [] })
			end
		end
	end

end
