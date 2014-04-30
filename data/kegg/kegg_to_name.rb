#This script produces a name mapping file from the KEGG ID of the source and target nodes to their gene names
require 'net/http'
network = ARGV[0]
@done = {}

def translate fin, fout
	until (line = fin.gets).nil?
		node = line.strip
		next if node.empty? or @done[node]
		@done[node] = true
		fout.print "#{node} "
		print "--"
		fout.puts node.split("|").map{ |id|
			print "."
			Net::HTTP.get(URI.parse "http://rest.kegg.jp/get/#{id}").split("\n").select{|l| l =~ /^NAME/}.first.strip.delete("NAME").delete(" ").split(",")
		}.flatten.uniq.join("|")
	end
end

open "#{network}/#{network}_id_map.txt", "w" do |fout|
	open("#{network}/#{network}_hsa_sources.txt"){|fin| translate fin, fout}
	print ">>>"
	open("#{network}/#{network}_hsa_targets.txt"){|fin| translate fin, fout}
end