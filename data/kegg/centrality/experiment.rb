dataset = ARGV[0]
sources = Dir.glob("#{dataset}/#{dataset}_hsa_source*.txt").map{|file| file.split("source").last.split(".").first.to_i}.sort
targets = Dir.glob("#{dataset}/#{dataset}_hsa_target*.txt").map{|file| file.split("target").last.split(".").first.to_i}.sort

nodes = open("#{dataset}/#{dataset}_hsa_nodes.txt"){|f| f.read.split("\n").size}.times.map{|i| i+1}
col_labels = nodes.map{|i| "'#{i}'"}.join " "
row_labels = ""
matrix = ""

open "#{dataset}/matrix.out", "w" do |f|
end

Dir.glob("#{dataset}/*").select{|f| File.directory? f}.map{|f| f.split("/").last.split}.each do |subtype|
	row_labels = "#{row_labels} '#{subtype}'"
	reference = {}
	sources.each do |source|
		reference[source] = {}
		targets.each do |target|
			puts "Running #{subtype}: reference from #{source} to #{target}"
			output = `./Graph #{dataset}/#{subtype}/#{dataset}_#{subtype}.txt #{dataset}/#{dataset}_hsa_source#{source}.txt #{dataset}/#{dataset}_hsa_target#{target}.txt pmc pre`
			reference[source][target] = output.split(">>").last.strip.to_f
		end
	end

	values = {}
	nodes.each do |node|
		values[node] = 0.0
		sources.each do |source|
			targets.each do |target|
				puts "Running #{subtype}: minus#{node} from #{source} to #{target}"
				output = `./Graph #{dataset}/#{subtype}/#{dataset}_#{subtype}_minus#{node}.txt #{dataset}/#{dataset}_hsa_source#{source}.txt #{dataset}/#{dataset}_hsa_target#{target}.txt pmc pre`
				new_value = output.split(">>").last.strip.to_f
				values[node] = values[node] + reference[source][target] - new_value
			end
		end
	end
	
	line = "#{nodes.map{|n| values[n]}.join(" ")}"
	open("#{dataset}/matrix.out", "a"){|f| f.puts line}
	matrix = "#{matrix} #{line} \n"

end

open "#{dataset}/centralgram.m", "w" do |fout|
	fout.puts "RowLabels = {#{row_labels}}"
	fout.puts "ColumnLabels = {#{col_labels}}"
	fout.puts "Matrix = [#{matrix}]"
	fout.puts "CG = clustergram(Matrix, 'RowLabels', RowLabels, 'ColumnLabels', ColumnLabels, 'Symmetric', false, 'ColorMap', colormap('JET'))"
end