#This script uses the expression values of leukemia subtypes
#to output correlation values between sources and targets of the given dataset
@dataset = ARGV[0]
@id_map = {}
@missing = {}

def fill_id_map
	open "#{@dataset}/#{@dataset}_id_map.txt" do |f|
		until (line = f.gets).nil?
			k, v = line.strip.split
			@id_map[k] = v.split "|"
		end
	end
end

def load_expression filename
	map = {}
	open filename do |f|
		until (line = f.gets).nil?
			parts = line.strip.split
			map[parts.first] = parts[1..parts.size-1].map{|v| v.to_f}
		end
	end
	map
end

def nodes_expression nodesfile, exp_map
	map = {}
	open nodesfile do |f|
		until (line = f.gets).nil?
			next if line.strip.empty?
			keggs = line.strip.split "|"
			records = []
			miss = []
			keggs.each do |kegg|
				kegg_records = @id_map[kegg].map{|g| exp_map[g]}
				if kegg_records.any?
					records = records + kegg_records.select{|r| not r.nil?}
				else
					miss << kegg
				end
			end
			@missing["#{nodesfile}:#{line.strip}"] = miss unless miss.empty?
			map[line.strip] = records.transpose.map{|a| a.reduce(:+)/a.size}
		end
	end
	map
end

def pearson x, y
	xbar, ybar = x.reduce(:+)/x.size, y.reduce(:+)/y.size
	numerator = [x.map{|i| i-xbar}, y.map{|i| i-ybar}].transpose.map{|pair| pair.reduce(:*)}.reduce(:+)
	denomerator = Math.sqrt(x.map{|i| (i-xbar)*(i-xbar)}.reduce(:+)) * Math.sqrt(y.map{|i| (i-ybar)*(i-ybar)}.reduce(:+))
	corr = (numerator/denomerator).abs
	raise "pearson value out of bounds: #{corr}" if corr > 1.0000001
	corr
end



fill_id_map
Dir.glob("leukemia/*.ge").each do |exp_file|
	exp = load_expression exp_file
	source_exp = nodes_expression "#{@dataset}/#{@dataset}_hsa_sources.txt", exp
	target_exp = nodes_expression "#{@dataset}/#{@dataset}_hsa_targets.txt", exp
	open "#{@dataset}/coexp_#{exp_file.split("/").last.split(".").first}.txt", "w" do |f|
		source_exp.keys.each do |s|
			next if source_exp[s].empty?
			target_exp.keys.each do |t|
				next if target_exp[t].empty?
				f.puts "#{s} #{t} #{pearson source_exp[s], target_exp[t]}"
			end
		end
	end
end

open("#{@dataset}/missing.txt", "w"){|f| @missing.each{|k,v| f.puts "#{k}: #{v.join ","}"}}


