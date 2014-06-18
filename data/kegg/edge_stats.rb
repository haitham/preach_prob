#This script reads all prob output files in subdir
#and produces some stats about every edge across all leukemia subtypes

def stdev(list)
	mean = list.reduce(:+)/list.size
	Math.sqrt(list.map{|x| (x - mean)*(x - mean)}.reduce(:+)/list.size)
end

def outlier(record)
	sum = record.values.reduce(:+)
	max_key, max_val = "", 0.0
	d = 0.0000001
	record.each do |key, value|
		mean = (sum - value)/(record.size-1)
		distance = Math.log (value+d)/(mean+d)
		max_key, max_val = key, distance if distance.abs > max_val.abs
	end
	return max_key, max_val
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
		stats[edge][:outlier_key], stats[edge][:outlier_value] = outlier record
	end

	open "#{subdir}/edge_stats.csv", "w" do |f|
		f.puts "#edge,stdev,outlier_subtype,outlier_distance"
		stats.each{|e, s| f.puts "#{e},#{s[:stdev]},#{s[:outlier_key]},#{s[:outlier_value]}"}
	end
end