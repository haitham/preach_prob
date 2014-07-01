#This script reads all prob output files in subdir
#and produces some stats about every edge across all leukemia subtypes
require 'net/http'
OUTLIER_THRESHOLD = 2.7

def stdev(list)
	mean = list.reduce(:+)/list.size
	Math.sqrt(list.map{|x| (x - mean)*(x - mean)}.reduce(:+)/list.size)
end

def outlier(record, dev)
	sum = record.values.reduce(:+)
	max_key, max_val = "", 0.0
	record.each do |key, value|
		mean = (sum - value)/(record.size-1)
		distance = (value - mean)/dev
		max_key, max_val = key, distance if distance.abs > max_val.abs
	end
	return max_key, max_val
end

def kegg_to_name(kegg)
	Net::HTTP.get(URI.parse "http://rest.kegg.jp/get/#{kegg}").split("\n").select{|l| l =~ /^NAME/}.first.strip.gsub("NAME", "").gsub(/\s+/, "").split(",").first
end

subdirs = ARGV
subdirs.each do |subdir|
	probs = {}
	init = true
	Dir.glob("#{subdir}/prob_*.out").each do |filename|
		subtype = filename.split("prob_").last.split(".out").first
		open filename do |f|
			until (line = f.gets).nil?
				next if line.strip.empty? or line =~ /(^Result)|(UNKNOWN)/
				s, t, p = line.strip.split
				probs["#{s}==>#{t}"] = {} if init
				probs["#{s}==>#{t}"][subtype] = p.to_f
			end
		end
		init = false
	end

	stats = {}
	probs.each do |edge, record|
		stats[edge] = {:stdev => stdev(record.values)}
		stats[edge][:outlier_key], stats[edge][:outlier_value] = outlier record, stats[edge][:stdev]
	end

	open "#{subdir}/edge_stats.csv", "w" do |f|
		f.puts "#edge,stdev,outlier_subtype,outlier_distance"
		stats.each{|e, s| f.puts "#{e},#{s[:stdev]},#{s[:outlier_key]},#{s[:outlier_value]}"}
	end
	
	puts subdir
	tops = []
	stats.each do |edge, record|
		if record[:outlier_value].abs >= OUTLIER_THRESHOLD
			tops << {:sign => (record[:outlier_value] > 0 ? "H." : "L."),
					 :from => kegg_to_name(edge.split("==>").first.split("|").first),
					 :to => kegg_to_name(edge.split("==>").last.split("|").first),
					 :subtype => record[:outlier_key],
					 :value => record[:outlier_value].abs}
		end
	end
	
	open "outliers.txt", "a" do |f|
		f.puts subdir
		f.puts subdir.size.times.map{|i| "-"}.join
		tops.sort{|a,b| b[:value] <=> a[:value]}.each{|t| f.puts "#{t[:sign]}     #{sprintf "%-8s", t[:from]} --> #{sprintf "%8s", t[:to]}     #{t[:subtype]}"}
		f.puts
	end
end





