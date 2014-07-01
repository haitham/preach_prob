#This script reads centrality matrix with column and row labels
#and produces some stats about every node across all leukemia subtypes
require 'net/http'
OUTLIER_THRESHOLD = 2.5

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
	nodes = open("#{subdir}/col_labels.out"){|f| f.read.strip.split.map{|v| v.gsub("'", "")}}
	subtypes = open("#{subdir}/row_labels.out"){|f| f.read.strip.split.map{|s| s.gsub("'", "")}}
	matrix = open("#{subdir}/matrix.out"){|f| f.read.strip.split("\n").map{|l| l.strip.split.map{|v| v.to_f}}}
	centrality = {}
	nodes.each_with_index do |v, i|
		centrality[v] = {}
		subtypes.each_with_index do |s, j|
			centrality[v][s] = matrix[j][i]
		end
	end

	stats = {}
	centrality.each do |node, record|
		stats[node] = {:stdev => stdev(record.values)}
		stats[node][:outlier_key], stats[node][:outlier_value] = outlier record, stats[node][:stdev]
	end

	open "#{subdir}/node_stats.csv", "w" do |f|
		f.puts "#node,stdev,outlier_subtype,outlier_distance"
		stats.each{|v, s| f.puts "#{v},#{s[:stdev]},#{s[:outlier_key]},#{s[:outlier_value]}"}
	end
	
	nodeids = {}
	open "#{subdir}/#{subdir}_hsa_nodes.txt" do |f|
		counter = 1
		until (line = f.gets).nil?
			nodeids[counter.to_s] = line.strip.split("|").first
			counter = counter + 1
		end
	end
	
	puts subdir
	tops = []
	stats.each do |node, record|
		if record[:outlier_value].abs >= OUTLIER_THRESHOLD
			tops << {:sign => (record[:outlier_value] > 0 ? "H." : "L."),
					 :node => kegg_to_name(nodeids[node]),
					 :subtype => record[:outlier_key],
					 :value => record[:outlier_value].abs}
		end
	end
	
	open "outliers.txt", "a" do |f|
		f.puts subdir
		f.puts subdir.size.times.map{|i| "-"}.join
		tops.sort{|a,b| b[:value] <=> a[:value]}.each{|t| f.puts "#{t[:sign]}     #{sprintf "%-8s", t[:node]}     #{t[:subtype]}"}
		f.puts
	end
end