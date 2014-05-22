#This script reads the network, and current sources and targets,
#and augment the sources with the the zero-indegree nodes,
#and targets with zero-outdegree nodes

def load_list filename
	list = []
	open filename do |f|
		until (line = f.gets).nil?
			next if line.strip.empty?
			list << line.strip
		end
	end
	list
end

dataset = ARGV[0]
sources = load_list "#{dataset}/#{dataset}_hsa_sources.old.txt"
targets = load_list "#{dataset}/#{dataset}_hsa_targets.old.txt"

indeg, outdeg = {}, {}
open "#{dataset}/#{dataset}_hsa.txt" do |f|
	until (line = f.gets).nil?
		next if line.strip.empty?
		s, t = line.strip.split
		outdeg[s] = 0 if outdeg[s].nil?
		outdeg[t] = 0 if outdeg[t].nil?
		indeg[s] = 0 if indeg[s].nil?
		indeg[t] = 0 if indeg[t].nil?
		outdeg[s] = outdeg[s] + 1
		indeg[t] = indeg[t] + 1
	end
end

sources = (sources + indeg.select{|k,v| v == 0}.map{|k, v| k}).uniq
targets = (targets + outdeg.select{|k,v| v == 0}.map{|k, v| k}).uniq

open("#{dataset}/#{dataset}_hsa_sources.txt", "w"){|f| sources.each{|s| f.puts s}}
open("#{dataset}/#{dataset}_hsa_targets.txt", "w"){|f| targets.each{|t| f.puts t}}