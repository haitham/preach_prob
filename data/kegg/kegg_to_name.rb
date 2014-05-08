#This script produces a name mapping file from the KEGG ID of the source and target nodes to their gene names
require 'net/http'
network = ARGV[0]
@done = {}

def translate fin, fout
	until (line = fin.gets).nil?
		node = line.strip
		next if node.empty?
		node.split("|").each do |kegg|
			next if @done[kegg]
			@done[kegg] = true
			fout.puts "#{kegg} #{Net::HTTP.get(URI.parse "http://rest.kegg.jp/get/#{kegg}").split("\n").select{|l| l =~ /^NAME/}.first.strip.gsub("NAME", "").gsub(/\s+/, "").split(",").uniq.join("|")}"
			print "."
		end
		print "--"
	end
end

open "#{network}/#{network}_id_map.txt", "w" do |fout|
	open("#{network}/#{network}_hsa_sources.txt"){|fin| translate fin, fout}
	print ">>>"
	open("#{network}/#{network}_hsa_targets.txt"){|fin| translate fin, fout}
end