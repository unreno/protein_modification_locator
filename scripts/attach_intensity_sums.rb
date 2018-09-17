#!/usr/bin/env ruby

require 'csv'
require 'optparse'
#require 'parallel'
#require File.join(File.dirname(__FILE__),'evidence.rb')

def usage(options={})
	puts
	puts "Usage"
	puts
	puts "Produces MatchedModificationWithIntensitySums.txt"
	puts
	puts "#{$0} [options]"
	puts
	puts "Options:"
	puts "	--evidence TEXT TSV FILE"
	puts
	puts "Defaults/current values:"
	puts
	puts "Expects the presence of MatchedModification.txt"
	puts
	puts "Examples:"
	puts "#{$0} --evidence evidence.txt"
	puts
	exit
end


######################################################################

#	Must be called before option parsing as they remove the items.
puts "Command: #{$0} #{$*.join(' ')}"


options = {}
OptionParser.new do |opt|
	opt.on("--evidence TEXT TSV FILE"){|o| options[:evidence_file] = o }
end.parse!

puts options

######################################################################

if options[:evidence_file].to_s.empty? or !File.exists?(options[:evidence_file])
	puts "\nEvidence file Required"
	usage(options)
end

#usage(options) #if ARGV.length < 2

puts "Starting"
puts Time.now


##################################################
#
#		Read evidence file.
#		Read MatchedModification.txt file.
#		Create MatchedModificationWithIntensitySums.txt file.
#
##################################################

#	Given that I don't use all of the features of evidences, just use a hash.
#	Its much faster to load and select anyway.

#evidences = []
evidences = {}

#c = CSV.open(options[:evidence_file],'rb')
#record_count = c.readlines.size		#	this is slow when file is huge ( 30 seconds / 5 million )
#c.close
record_count = `wc -l #{options[:evidence_file]}`.split()[0].to_i

(c=CSV.open(options[:evidence_file],'rb',
		{ headers: true, col_sep: "\t" })).each do |line|

	puts "#{Time.now} : Processing record #{c.lineno}/#{record_count} : #{line['Modified sequence']} : #{line['Leading razor protein']}"

	#puts ['Modified sequence','Experiment','Intensity'].collect{|f| line[f] }.join(",")

#	evidences.push Evidence.new(line.to_hash)

	evidences[line['Modified sequence']] ||={}
	evidences[line['Modified sequence']][line['Experiment']] ||= []
	evidences[line['Modified sequence']][line['Experiment']].push(line['Intensity'].to_i)

end	#	(c=CSV.open(options[:evidence_file],'rb',

puts "Generating experiment list."
#experiments = evidences.collect{|e| e.experiment }.sort.uniq

experiments = evidences.keys.collect{|s|evidences[s].keys}.flatten.sort.uniq
#puts experiments



puts "Writing MatchedModificationWithIntensitySums.txt ..."


c = CSV.open("MatchedModification.txt",'rb',{col_sep: "\t" })
matched_modification_header = c.readline
c.close

#c = CSV.open("MatchedModification.txt",'rb',{col_sep: "\t" })
#record_count = c.readlines.size
#c.close
record_count = `wc -l MatchedModification.txt`.split()[0].to_i

matched_mod = CSV.open("MatchedModificationWithIntensitySums.txt",'w', {col_sep: "\t" })

matched_mod_header = matched_modification_header + experiments

matched_mod.puts matched_mod_header

(c=CSV.open("MatchedModification.txt",'rb',
		{ headers: true, col_sep: "\t" })).each do |line|

	puts "#{Time.now} : Processing record #{c.lineno}/#{record_count} : #{line['Modified sequence']}"

	matched_mod_line = line.fields

	sequence_matches = evidences[line['Modified sequence']]

	experiments.each do |exp|
		puts exp
		intensities = sequence_matches[exp] || []

		#	with "select" use {/} and NOT do/end

		#matches = evidences.select { |e| 
		#	e.sequence == line['Modified sequence'] && e.experiment == exp
		#}
		#matches = evidences[line['Modified sequence']][exp] || []

		puts "Found #{intensities.length} match intensities. Summing ..."

		#matched_mod_line << matches.collect {|e| e.intensity.to_i }.inject(:+) || 0 )
		matched_mod_line << ( intensities.inject(:+) || 0 )

	end	#	experiments.each do |exp|

	matched_mod.puts matched_mod_line

end	#	(c=CSV.open("MatchedModification.txt",'rb',

matched_mod.close

##################################################

puts "Done"
puts Time.now




#	I need to attach intensity sum to the file you previously generated (MatchedModification.txt). I shared those files in the Nevada box with you.
#	
#	For MatchedModification.txt, we have
#	
#	Modified sequence, Leading razor protein, Start position, End position, STY position, STY modification, M position, M modification position
#	
#	 
#	
#	For evidence.txt, there are lots of columns, but the following columns are our interest:
#	
#	Modified sequence, Experiment, Intensity
#	
#	 
#	
#	I want to add sum of intensities of a given modified sequence for each experiment.
#	
#	For example, we have experiment C1, C2, C3, C4, C5, C6, ….
#	
#	For the modified _(ac)AAAAAAAAAAGAAGGRGS(ph)GPGRR_, find intensities of this modified sequence in evidence file for C1 and add all the corresponding intensities and attached it to MatchedModification.txt. This sum will be an entry for C1.
#	
#	 
#	
#	So at the end, the output file contains the following columns
#	
#	Modified sequence, Leading razor protein, Start position, End position, STY position, STY modification, M position, M modification position, C1, C2, C3, C4, C5, C6, E15_1, E15_2, E15_3, E15_4, N1, N2, N3, …
#	
#	 
#	
#	If there is no intensity found for a particular peptide for a given experiment, then you can enter zero.
