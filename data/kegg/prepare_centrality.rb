#This script prepares the centrallity directory with copies of all networks
#and copies of each network for every subtype, and copies of every subtype
#for every missing node
require 'fileutils'
subdirs = ARGV
subdirs.each do |subdir|
	#Copy files
	FileUtils.mkdir "centrality/#{subdir}"
	open "#{subdir}/#{subdir}_hsa_sources.txt" do |fin|
		counter = 0
		until (line = fin.gets).nil?
			next if line.strip.empty?
			open("centrality/#{subdir}/#{subdir}_hsa_source#{counter}.txt", "w"){|fout| fout.puts line}
			counter = counter + 1
		end
	end
	open "#{subdir}/#{subdir}_hsa_targets.txt" do |fin|
		counter = 0
		until (line = fin.gets).nil?
			next if line.strip.empty?
			open("centrality/#{subdir}/#{subdir}_hsa_target#{counter}.txt", "w"){|fout| fout.puts line}
			counter = counter + 1
		end
	end
	
	#generate nodes file
	nodes = {}
	open "#{subdir}/#{subdir}_hsa.txt" do |fin|
		open "centrality/#{subdir}/#{subdir}_hsa_nodes.txt", "w" do |fout|
			node_counter = 1
			until (line = fin.gets).nil?
				next if line.strip.empty?
				s, t, p = line.strip.split
				if nodes[s].nil?
					nodes[s] = node_counter
					node_counter = node_counter + 1
					fout.puts s
				end
				if nodes[t].nil?
					nodes[t] = node_counter
					node_counter = node_counter + 1
					fout.puts t
				end
			end
		end
	end
	
	#create the subtype network file + files for missing nodes
	Dir.glob("#{subdir}/prob_*.out").each do |subfile|
		subtype = subfile.split("prob_").last.split(".out").first
		FileUtils.mkdir "centrality/#{subdir}/#{subtype}"
		edges = []
		open subfile do |f|
			until (line = f.gets).nil?
				next if line.strip.empty?
				edges << line.strip unless line =~ /(UNKNOWN)|(Result)/
			end
		end
		open("centrality/#{subdir}/#{subtype}/#{subdir}_#{subtype}.txt", "w"){|f| edges.each{|e| f.puts e}}
		nodes.each do |node, index|
			open "centrality/#{subdir}/#{subtype}/#{subdir}_#{subtype}_minus#{index}.txt", "w" do |f|
				edges.each do |e|
					s, t, p = e.strip.split
					f.puts e unless s == node or t == node
				end
			end
		end
	end
end










