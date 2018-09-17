#!/usr/bin/env ruby

require 'csv'
require 'optparse'

require 'parallel'

def usage(options={})
	puts
	puts "Usage"
	puts
	puts "Produces MatchedModification.txt, ProteinModification.txt, "
	puts "MatchedModificationNONE.txt, MatchedModificationMULTIPLE.txt"
	puts "and ModificationNormalizerPeptides.txt"
	puts
	puts "#{$0} [options] "#<formatted moda output file> <observed spectra mgf file>"
	puts
	puts "Options:"
	puts "	--amino_acids COMMA SEPARATED TEXT"
	puts "	--evidence TEXT TSV FILE"
	puts "	--protein FASTA FILE"
	puts
	puts "Defaults/current values:"
	puts
	puts "Examples:"
	puts "#{$0} --amino_acids STY,M --evidence evidence.txt --protein uniprot-organism+homo+sapiens.fasta"
	puts
	exit
end


######################################################################


#	Setup future redirect of STDOUT if desired
#stdout_redirected = false
#original_stdout = $stdout.clone
#	In the future, could redirect STDOUT to /dev/null
#	if !stdout_redirected && options.has_key?(:log_count) && record_number > options[:log_count].to_i
#		puts "Record number(#{record_number}) has exceeded requested log count(#{options[:log_count]})."
#		puts "Redirecting the rest of the output to /dev/null."
#		stdout_redirected = true
#		$stdout.reopen("/dev/null", "a")
#	end


#	Must be called before option parsing as they remove the items.
puts "Command: #{$0} #{$*.join(' ')}"


options = {}
OptionParser.new do |opt|
	opt.on("--amino_acids COMMA SEPARATED TEXT"){|o| options[:amino_acids] = o }
	opt.on("--evidence TEXT TSV FILE"){|o| options[:evidence_file] = o }
	opt.on("--protein FASTA FILE"){|o| options[:protein_file] = o }
end.parse!

puts options

######################################################################

if options[:amino_acids].to_s.empty?
	puts "\nAmino Acids Required"
	usage(options)
end
if options[:evidence_file].to_s.empty? or !File.exists?(options[:evidence_file])
	puts "\nEvidence file Required"
	usage(options)
end
if options[:protein_file].to_s.empty? or !File.exists?(options[:protein_file])
	puts "\nProtein fasta file Required"
	usage(options)
end

#usage(options) #if ARGV.length < 2

puts "Starting"
puts Time.now


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

##################################################
#
#		Create protein database from protein file
#
##################################################

protein_base_name=options[:protein_file].gsub(/.[^.]*$/,'').gsub(/.*\/(.*)$/,'\1')

if !Dir.exists?( protein_base_name )
	puts "Creating #{protein_base_name} protein database."
	system("java -jar PeptideMatchCMD_1.0.jar -a index -d #{options[:protein_file]} -i #{protein_base_name}")
else
	puts "#{protein_base_name} protein database exists."
end



##################################################
#
#		Read evidence file.
#		Stop creating select_evidence_file.
#		Create MatchedModification*.txt files.
#		Create ModificationNormalizerPeptides.txt
#		Create experiments array.
#
##################################################


puts "Producing ModificationNormalizerPeptides.txt ..."
modification_normalizer=CSV.open("ModificationNormalizerPeptides.txt",'w', {col_sep: "\t" })
modification_normalizer.puts ['Leading razor protein','Peptides without specified AAs']
peptides_without_specified_aas={}

no_matches = CSV.open("MatchedModificationNONE.txt",'w', {col_sep: "\t" })
no_matches.puts ['Modified sequence','Leading razor protein']

multiple_matches = CSV.open("MatchedModificationMULTIPLE.txt",'w', {col_sep: "\t" })
multiple_matches.puts ['Modified sequence','Leading razor protein']

matched_mod = CSV.open("MatchedModification.txt",'w', {col_sep: "\t" })
matched_mod_header = ['Modified sequence','Leading razor protein','Start position','End position']

options[:amino_acids].split(",").each do |acid|
	matched_mod_header.push "#{acid} position"
	matched_mod_header.push "#{acid} modification position"
end
matched_mod.puts matched_mod_header

evidences = []

acids = options[:amino_acids].gsub(/[^[:alpha:]]/,'')
acids_regex=Regexp.new("[#{acids}]")

c = CSV.open(options[:evidence_file],'rb')
record_count = c.readlines.size
c.close

(c=CSV.open(options[:evidence_file],'rb',
		{ headers: true, col_sep: "\t" })).each do |line|

	puts "#{Time.now} : Processing record #{c.lineno}/#{record_count} : #{line['Modified sequence']} : #{line['Leading razor protein']}"

	e = Evidence.new(line.to_hash)
	evidences.push e

	if line['Modified sequence'].match(acids_regex)
		#	MatchedModification*.txt

		e.query( protein_base_name )
		e.matched_amino_acids = true

		if e.alignments.length < 1
			puts "None of the matches match the expected protein."
			no_matches.puts line
		else
			if e.alignments.length > 1
				puts "More than one matches match the expected protein. Reporting all #{e.alignments.length}."
				multiple_matches.puts line
			end

			e.alignments.each do |alignment|

				matched_mod_line = [e.sequence,e.protein,alignment[:match_start],alignment[:match_end]]

				acids.split(//).each do |acid|
					positions = e.cleaned_sequence.indices_of_chars(acid).collect{|i| alignment[:match_start] + i }
					positions.each { |position| alignment[:absolute_positions][position] = acid }
				end

				options[:amino_acids].split(",").each do |acid_group|

					absolute_positions=alignment[:absolute_positions].select{|k,v| acid_group.match(v) }.collect{|k,v| k }

					modified_positions=e.modified_positions.collect{|i| alignment[:match_start] + i } & absolute_positions
					alignment[:modified_positions] += modified_positions

					matched_mod_line << absolute_positions.join(";")
					matched_mod_line << modified_positions.join(";")
				end

				matched_mod.puts matched_mod_line

			end	#	e.alignments.each do |alignment|

		end	#	e.alignments.length >= 1

	else
		#	ModificationNormalizerPeptides.txt

		peptides_without_specified_aas[ line['Leading razor protein'] ] ||= []
#		unless peptides_without_specified_aas[ line['Leading razor protein'] ].include? line['Modified sequence']
			peptides_without_specified_aas[ line['Leading razor protein'] ] << line['Modified sequence']
#		end

	end


end	#	(c=CSV.open(options[:evidence_file],'rb',

puts "Writing ModificationNormalizerPeptides.txt ..."

peptides_without_specified_aas.keys.sort.each do |protein|
	#	mnp
	peptides_without_specified_aas[protein] = peptides_without_specified_aas[protein].sort.uniq
	modification_normalizer.puts [ protein, peptides_without_specified_aas[protein].join(';') ]
end
modification_normalizer.close

puts "Generating experiment list."
experiments = evidences.collect{|e| e.experiment }.sort.uniq

puts "Selecting evidences"
select_evidences = evidences.select{|e| e.matched_amino_acids == true }

no_matches.close
multiple_matches.close
matched_mod.close


##################################################
#
#		Loop over each protein and create ProteinModification.txt
#
##################################################

protein_mod=CSV.open("ProteinModification.txt",'w', {col_sep: "\t" })
protein_mod_line = ['Leading razor protein',
	'AA',
	'AA location',
	'Sequences with Modified AA',
	'Sequences with Unmodified AA']
experiments.each do |k|
	protein_mod_line.push(
		"Sequences with Modified AA intensity.#{k}",
		"Sequences with Unmodified AA intensity.#{k}",
		"ModificationNormalizerPeptides intensity.#{k}",
		"Sequences with Modified AA sum.log2.intensity.#{k}",
		"Sequences with Unmodified AA sum.log2.intensity.#{k}",
		"ModificationNormalizerPeptides sum.log2.intensity.#{k}",
		"Sequences with Modified AA MSMSCount.#{k}",
		"Sequences with Unmodified AA MSMSCount.#{k}",
		"ModificationNormalizerPeptides MSMSCount.#{k}"
	)
end
protein_mod.puts protein_mod_line


#	Leading razor protein	AA	AA location	Sequences with Modified AA	Sequences with Unmodified AA
#
#	PROTEIN1	S	4		_ASSST[80]TY_;_SSS[80]T[80]_
#	PROTEIN1	S	5		_ASSST[80]TY_;_SSS[80]T[80]_
#	PROTEIN1	S	6	_SSS[80]T[80]_	_ASSST[80]TY_
#	PROTEIN1	T	7	_ASSST[80]TY_;_SSS[80]T[80]_


puts "Producing ProteinModification.txt ..."

puts "Selecting proteins."
proteins = select_evidences.select{|e|e.alignments.length > 0}.collect{|e|e.protein}.uniq.sort

puts "Looping over all proteins."

#proteins.each_with_index do |protein, protein_i|

#	Use parallel as this process would take weeks without for large 5,000,000 line evidence file
#	Takes days (week) on 64 CPU Azure machine anyway.

Parallel.each_with_index(proteins) do |protein, protein_i|

	puts "Processing protein #{protein_i+1}/#{proteins.length} : #{protein}"

	puts "Selecting evidences for this protein."
	select_evidences_for_this_protein = select_evidences.select{|e|e.protein == protein}

	positions = {}

	alignments = select_evidences_for_this_protein.collect do |e|
		e.alignments.each do |alignment|
			positions.update alignment[:absolute_positions]

			#	copy sequence to alignment
			alignment[:sequence] = e.sequence

			#	copy intensity to alignment
			alignment[:intensity] = e.intensity

			#	copy experiment to alignment
			alignment[:experiment] = e.experiment

			#	copy msmscount to alignment
			alignment[:msmscount] = e.msmscount

		end
		e.alignments
	end.flatten

	positions_length = positions.keys.length
	positions.keys.sort.each_with_index do |position,position_i|

		puts "Processing position #{position_i+1}/#{positions_length} : #{position}."

		acid=positions[position]

		alignments_at_this_position = alignments.select{|a| a[:match_start] <= position && a[:match_end] >= position }

		with_modified_positions = alignments_at_this_position.select{|a|
			a[:modified_positions].include? position }
		sequences_with_modified_positions = with_modified_positions.collect{|a|
			a[:sequence] }.sort.uniq.join(';')
		sequences_with_modified_positions = "NA" if sequences_with_modified_positions.empty?

		without_modified_positions = alignments_at_this_position.select{|a|
			!a[:modified_positions].include? position }
		sequences_without_modified_positions = without_modified_positions.collect{|a|
			a[:sequence] }.sort.uniq.join(';')
		sequences_without_modified_positions = "NA" if sequences_without_modified_positions.empty?

		protein_mod_line = [ protein, acid, position,
			sequences_with_modified_positions, sequences_without_modified_positions ]

		experiments.each_with_index do |experiment,experiment_i|

			puts "Processing experiment #{experiment_i+1}/#{experiments.length} : #{experiment}."

			exp_with_modified_positions = with_modified_positions.select{|a|
				a[:experiment] == experiment }
			exp_without_modified_positions = without_modified_positions.select{|a|
				a[:experiment] == experiment }

			#	protein_mod_line.push("Sequences with Modified AA intensity.#{experiment}",
			#		"Sequences with Unmodified AA intensity.#{experiment}",
			#		"ModificationNormalizerPeptides intensity.#{experiment}")


			#	<Extract intensities of "Sequences with Modified AA" for each Experiment and
			#	attach them to ProteinModification.txt>
			#
			#	We will have columns "Sequences with Modified AA intensity.[ExperimentName]".
			#	For example, Sequences with Modified AA intensity.C1, Sequences with Modified AA intensity.C2, and so on.
			#
			#	"experimentalDesignTemplate_FINAL.txt" can be used to create all the column names.
			#
			#	If there is no Sequences with Modified AA, then please list NA.
			#
			#	If we have any Sequences with Modified AA, please sum their intensities from
			#	"evidence.txt" for the corresponding Experiment.
			#
			#	For example, let say we have _MMMS(ph)M_. And we have 3 experiments, Y1, Y2, and Y3.
			#	Then look for intensities of _MMMS(ph)M_ for Y1. They may be more than one, so please add
			#	them up and listed under the column name – "Modified AA intensity.Y1".
			#	If there is no intensity, then put NA. You can do the same for Y2, and Y3.

			protein_mod_line.push( exp_with_modified_positions.collect{|a|
				a[:intensity].to_i }.inject(:+) || "NA" )



			#	<Extract intensities of "Sequences with Unmodified AA" for each Experiment and
			#	attach them to ProteinModification.txt>
			#
			#	Same thing as above. The column name can be "Sequences with Unmodified AA intensity.[ExperimentName]"

			protein_mod_line.push( exp_without_modified_positions.collect{|a|
				a[:intensity].to_i }.inject(:+) || "NA" )



			#	<Extract intensities of "ModificationNormalizerPeptides.txt" for each Experiment and
			#	attached them to ProteinModificaiton.txt>
			#
			#	The column name can be "ModificationNormalizerPeptides intensity.[ExperimentName]"
			#
			#	I want to sum all intensities of "peptides without specified AAs" for each protein
			#	listed in "ProteinModificaiton.txt".
			#
			#	Let say that "ProteinModification.txt" has ProteinABC. Find "Peptide without specified AAs"
			#	of ProteinABC in "ModificationNormalizerPeptides.txt".
			#
			#	Then, sum intensities of all "Peptide without specified AAs" in "evidence.txt" for the
			#	corresponding Experiment and record it in ProteinModification.txt.
			#
			#
			#	For example,
			#
			#	"PorteinModification.txt"
			#
			#	Leading razor protein	AA	AA location	Sequences with Modified AA	Sequences with Unmodified AA
			#	ProteinABC	M	22	_	EEEEESEEK_
			#
			#
			#	"ModificationNormalizerPeptides.txt"
			#
			#	Leading razor protein	Peptides without specified AAs
			#	ProteinABC	_PLLVEPEGLEK_;_VLLVDGK_
			#
			#
			#	"evidence.txt"
			#
			#	Modified.sequence	Experiment	Intensity	MS.MS.Count
			#	_PLLVEPEGLEK_	Y1	10	1
			#	_PLLVEPEGLEK_	Y1	12	1
			#	_VLLVDGK_	Y1	100	3
			#
			#
			#	Then, for the column "ModificationNormalizerPeptides intensity.Y1", ProteinABC has 10+12+100=122.
			#
			#	"ProteinModificaiton.txt " has sp|A0A0B4J2F0|PIOS1_HUMAN in its first line.
			#	But I think "ModificationNormalizerPeptides.txt" does not have sp|A0A0B4J2F0|PIOS1_HUMAN.
			#
			#	So, NA will be listed.


			#puts "---"
			#puts "Protein:#{protein}:"

			#	mnp = peptides_without_specified_aas[protein] || []

			#	m1 = evidences.select{|e|
			#		e.protein == protein && e.experiment == experiment && mnp.include?( e.sequence ) }

			#	m2 = evidences.select{|e|
			#		e.protein == protein && e.experiment == experiment && e.matched_amino_acids == false }

			#	if !m1.empty? || !m2.empty?
			#		puts "Comparing"
			#		puts m1.inspect
			#		puts m2.inspect
			#		if m1 == m2
			#			puts "SAME"
			#		else
			#			puts "DIFFERENT"		#	never happened during test
			#		end
			#	end


			#puts protein
			#	tr|V9HW72|V9HW72_HUMAN
			#puts mnp.inspect
			#	["_AAALEFLNR_", "_AAALEFLNRFEEAK_"]
			#
			#	There is no position here so ???

#			mnp_evidences_for_these_sequences = evidences.select{|e|
#				e.protein == protein }.select{|a|
#				a.experiment == experiment }.select{|e|
#				mnp.include? e.sequence } #	this takes about 3 seconds
#	3 separate selects?

#	faster or slower to join into 1? neither it seems

#			mnp_evidences_for_these_sequences = evidences.select{|e|
#				e.protein == protein && e.experiment == experiment && mnp.include?( e.sequence ) }

#	why mnp and not just?
#	e.matched_amino_acids == false

#	all three seem to take the same amount of time

			mnp_evidences_for_these_sequences = evidences.select{|e|
				e.protein == protein && e.experiment == experiment && e.matched_amino_acids == false }


			#puts "experiment:#{experiment.to_s}"
			#puts "MNP:#{mnp.to_s}"
			#puts "AFTS:#{mnp_evidences_for_these_sequences.collect{|e| e.sequence }.sort.uniq.inspect}"

			protein_mod_line.push( mnp_evidences_for_these_sequences.collect{|a|
				a.intensity.to_i }.inject(:+) || "NA" )


			#	I wonder whether it is better to sum up the intensity like above or log2-transform & sum them up.
			#
			#	So please log2-transform peptide intensities before summing them up and sum
			#	log2-transformed intensities and make additional column.
			#	For the above example for ProteinABC, log2(10)+log2(12)+log2(100)=?
			#
			#	The column names will be something like "Sequences with Modified AA sum.log2.intensity.Y1", ...

			protein_mod_line.push( exp_with_modified_positions.collect{|a|
				Math.log2(a[:intensity].to_i) }.inject(:+) || "NA" )

			protein_mod_line.push( exp_without_modified_positions.collect{|a|
				Math.log2(a[:intensity].to_i) }.inject(:+) || "NA" )

			protein_mod_line.push( mnp_evidences_for_these_sequences.collect{|a|
				Math.log2(a.intensity.to_i) }.inject(:+) || "NA" )


			#	Now, MS.MS.Count in "evidence.txt" will be listed in the same way as intensity.
			#	Summing up MS.MS.Count for the corresponding column.
			#
			#	For example, we will have "Sequences with Modified AA MSMSCount.Y1",
			#	"Sequences with Modified AA MSMSCount.Y2", "Sequences with Modified AA MSMSCount.Y3",
			#
			#	"Sequences with Unmodified AA MSMSCount.Y1", "Sequences with Unmodified AA MSMSCount.Y2",
			#	"Sequences with Unmodified AA MSMSCount.Y3",
			#
			#	"ModificationNormalizerPeptides MSMSCount.Y1", "ModificationNormalizerPeptides MSMSCount.Y2",
			#	and "ModificationNormalizerPeptides MSMSCount.Y3".

			protein_mod_line.push( exp_with_modified_positions.collect{|a|
				a[:msmscount].to_i }.inject(:+) || "NA" )

			protein_mod_line.push( exp_without_modified_positions.collect{|a|
				a[:msmscount].to_i }.inject(:+) || "NA" )

			protein_mod_line.push( mnp_evidences_for_these_sequences.collect{|a|
				a.msmscount.to_i }.inject(:+) || "NA" )

		end

		protein_mod.puts protein_mod_line

	end

end	#	proteins.each_with_index do |protein, protein_i|

protein_mod.close

##################################################

puts "Done"
puts Time.now




#	awk -F"\t" '{print $(NF-180),$(NF-171),$(NF-162),$(NF-153),$(NF-144),$(NF-135),$(NF-126),$(NF-117),$(NF-108),$(NF-99),$(NF-90),$(NF-81),$(NF-72),$(NF-63),$(NF-54),$(NF-45),$(NF-36),$(NF-27),$(NF-18),$(NF-9),$NF}' ProteinModification.txt | uniq

#	
#	
#	 head -1 modification_locator/evidence.1000.txt | awk -F"\t" '{ for(i=1;i<=NF;i++)print i,$i}'
#	1 Sequence
#	2 Length
#	3 Modifications
#	4 Modified sequence
#	5 Oxidation (M) Probabilities
#	6 Phospho (STY) Probabilities
#	7 Oxidation (M) Score Diffs
#	8 Phospho (STY) Score Diffs
#	9 Acetyl (Protein N-term)
#	10 Oxidation (M)
#	11 Phospho (STY)
#	12 Missed cleavages
#	13 Proteins
#	14 Leading proteins
#	15 Leading razor protein
#	16 Type
#	17 Raw file
#	18 Fraction
#	19 Experiment
#	20 MS/MS m/z
#	21 Charge
#	22 m/z
#	23 Mass
#	24 Resolution
#	25 Uncalibrated - Calibrated m/z [ppm]
#	26 Uncalibrated - Calibrated m/z [Da]
#	27 Mass Error [ppm]
#	28 Mass Error [Da]
#	29 Uncalibrated Mass Error [ppm]
#	30 Uncalibrated Mass Error [Da]
#	31 Max intensity m/z 0
#	32 Retention time
#	33 Retention length
#	34 Calibrated retention time
#	35 Calibrated retention time start
#	36 Calibrated retention time finish
#	37 Retention time calibration
#	38 Match time difference
#	39 Match m/z difference
#	40 Match q-value
#	41 Match score
#	42 Number of data points
#	43 Number of scans
#	44 Number of isotopic peaks
#	45 PIF
#	46 Fraction of total spectrum
#	47 Base peak fraction
#	48 PEP
#	49 MS/MS Count
#	50 MS/MS Scan Number
#	51 Score
#	52 Delta score
#	53 Combinatorics
#	54 Intensity
#	55 Reverse
#	56 Potential contaminant
#	57 id
#	58 Protein group IDs
#	59 Peptide ID
#	60 Mod. peptide ID
#	61 MS/MS IDs
#	62 Best MS/MS
#	63 AIF MS/MS IDs
#	64 Oxidation (M) site IDs
#	65 Phospho (STY) site IDs
#	
