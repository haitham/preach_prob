#This script prepares the GSEA experiment files
#It generates one full expression file, and multiple phenotype label files
#It also prepares gene lists of every phenotype/network pair
require 'net/http'
def kegg_to_names(kegg)
	begin
		Net::HTTP.get(URI.parse "http://rest.kegg.jp/get/#{kegg}").split("\n").select{|l| l =~ /^NAME/}.first.strip.gsub("NAME", "").gsub(/\s+/, "").split(",")
	rescue
		puts ">>>>KEGGMISS"
		[]
	end
end

Z_THRESHOLD = 2.0
networks = ["apoptosis", "cellcycle", "p53", "ccc", "chemokine", "nfkappab", "erbb", "wnt"]
subtypes = []
labels = []
samples = {}
init_genes = true
Dir.glob("leukemia/*.ge").each do |subfile|
	subtype = subfile.split("/").last.split(".ge").first
	subtypes << subtype
	puts "reading: #{subtype}"
	open subfile do |f|
		sample_size = 0
		until (line = f.gets).nil?
			parts = line.strip.split
			samples[parts.first] = [] if init_genes
			sample_size = parts.size-1
			samples[parts.first] = samples[parts.first] + parts[1..parts.size-1]
		end
		labels << {:subtype => subtype, :size => sample_size}
	end
	init_genes = false
end

=begin
#Output the expression file
open "gsea/expr_all.gct", "w" do |f|
	f.puts "#1.2"
	f.puts "#{samples.size}\t#{labels.map{|l| l[:size]}.reduce(:+)}"
	f.puts "NAME\tDescription\t#{labels.map{|l| l[:size].times.map{|i| "#{l[:subtype]}-#{i}"}.join("\t")}.join("\t")}"
	samples.each{|gene, values| f.puts "#{gene}\tNA\t#{values.join("\t")}"}
end

#Output the label files
subtypes.each_with_index do |subtype, order|
	open "gsea/labels_#{subtype}.cls", "w" do |f|
		f.puts "#{labels.map{|l| l[:size]}.reduce(:+)} 2 1"
		if order.zero?
			f.puts "# #{subtype} OTHERS"
		else
			f.puts "# OTHERS #{subtype}"
		end
		f.puts "#{labels.map{|l| l[:size].times.map{|i| l[:subtype] == subtype ? subtype : "OTHERS"}.join(" ")}.join(" ")}"
		f.puts
	end
end

#Output the edge_outliers gene lists
puts "Working on edge outlier lists:"
genelists = networks.map{ |n| {n => subtypes.map{ |s| {s => []} }.reduce(:merge)} }.reduce(:merge)
networks.each do |network|
	puts "compliling gene lists: #{network}"
	open "#{network}/edge_stats.csv" do |f|
		f.gets #header
		until (line = f.gets).nil?
			edge, dev, subtype, zscore = line.strip.split ","
			if zscore.to_f.abs >= Z_THRESHOLD
				genelists[network][subtype] = genelists[network][subtype] + edge.split(/\||\=\=\>/).map{|kegg| kegg_to_names kegg}.reduce(:+).uniq
			end
		end
	end
end

open "gsea/edge_outliers/stats.out", "w" do |fstats|
	fstats.puts "#{sprintf "%-15s", "Subtype"}#{networks.map{|n| sprintf "%-15s", n}.join}"
	subtypes.each do |subtype|
		fstats.print sprintf("%-15s", subtype)
		open "gsea/edge_outliers/genelists_#{subtype}.gmt", "w" do |flists|
			networks.each do |network|
				puts "outputting list: #{subtype} - #{network}"
				genes = genelists[network][subtype].uniq & samples.keys
				fstats.print sprintf("%-15s", genes.size)
				flists.puts "#{network}\tNA\t#{genes.join("\t")}"
			end
			fstats.puts
		end
	end
end
=end

#Output the centrality_outliers gene lists
puts "Working on centrality outlier lists:"
genelists = networks.map{ |n| {n => subtypes.map{ |s| {s => []} }.reduce(:merge)} }.reduce(:merge)
networks.each do |network|
	puts "compliling gene lists: #{network}"
	nodeids = {}
	open "centrality/#{network}/#{network}_hsa_nodes.txt" do |f|
		counter = 1
		until (line = f.gets).nil?
			nodeids[counter.to_s] = line.strip.split("|").first
			counter = counter + 1
		end
	end
	open "centrality/#{network}/node_stats.csv" do |f|
		f.gets #header
		until (line = f.gets).nil?
			node, dev, subtype, zscore = line.strip.split ","
			if zscore.to_f.abs >= Z_THRESHOLD
				genelists[network][subtype] = genelists[network][subtype] + nodeids[node].split("|").map{|kegg| kegg_to_names kegg}.reduce(:+).uniq
			end
		end
	end
end

open "gsea/centrality_outliers/stats.out", "w" do |fstats|
	fstats.puts "#{sprintf "%-15s", "Subtype"}#{networks.map{|n| sprintf "%-15s", n}.join}"
	subtypes.each do |subtype|
		fstats.print sprintf("%-15s", subtype)
		open "gsea/centrality_outliers/genelists_#{subtype}.gmt", "w" do |flists|
			networks.each do |network|
				puts "outputting list: #{subtype} - #{network}"
				genes = genelists[network][subtype].uniq & samples.keys
				fstats.print sprintf("%-15s", genes.size)
				flists.puts "#{network}\tNA\t#{genes.join("\t")}"
			end
			fstats.puts
		end
	end
end














