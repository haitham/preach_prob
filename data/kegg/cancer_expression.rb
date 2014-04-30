#This script aggregates expression values for every probe id from samples available for the given cancer subtype
#And outputs them in one file, mapping the probe ids to gene names, and subtype number to subtype name in the filename
typenames = {"1" => "Hyperdiploid",
"2" => "TCF3-PBX1",
"3" =>	"ETV6_RUNX1",
"4" => "MLL",
"5" => "Ph",
"6" => "Hypo",
"8" => "T-ALL"}

#First, we load a probe id => gene name conversion map
name_map = {}
open "../CancerData/probe_name_map.csv" do |f|
	f.gets #ignore first line
	until (line = f.gets).nil?
		parts = line.strip.split(",")
		next if parts[3].strip.empty?
		probe = parts[0].strip
		gene = parts[3].strip.split(":").first.delete("\"").strip.split("///").first.strip
		name_map[probe] = gene
	end
end

#Now we do the aggregation for every type
typenames.keys.each do |type|
	probes = {}
	samples = Dir.glob "../CancerData/#{type}/*.tab"
	#Here just collecting probe names
	open samples.first do |f|
		until (line = f.gets).nil?
			probes[line.strip.split.first.strip.delete("\"")] = []
		end
	end
	raise if probes.keys.uniq.size != probes.keys.size #sanity check
	samples.each do |sample|
		open sample do |f|
			until (line = f.gets).nil?
				probe, value = line.strip.split
				probes[probe.delete("\"")] << value.to_f
			end
		end
	end
	raise if probes.values.map{|v| v.size != samples.size}.any? #sanity check
	#Now we aggregate probe values under gene names
	genes = {}
	probes.each do |p, v|
		gene = name_map[p]
		genes[gene] = [] if genes[gene].nil?
		genes[gene] << v
	end
	#Now we collapse multiple values for the same gene name by taking the higher value
	genes_uniq = {}
	genes.each{|g, v| genes_uniq[g] = v.transpose.map{|a| a.max}}
	#Now we output to a file
	open("leukemia/#{typenames[type]}.ge", "w"){ |f| genes_uniq.each{|g, v| f.puts "#{g}  #{v.join("  ")}"} }
end










