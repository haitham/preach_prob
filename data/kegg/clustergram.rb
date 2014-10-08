#This script reads all prob output files in subdir
#And produces a clustergram matlab file
subdir = ARGV[0]
matrix = ""
row_labels = ""
Dir.glob("#{subdir}/prob_*.out").each do |filename|
	row_labels = "#{row_labels} '#{filename.split("prob_").last.split(".out").first}'"
	open filename do |f|
		until (line = f.gets).nil?
			next if line.strip.empty? or line =~ /(^Result)|(UNKNOWN)/
			matrix = "#{matrix} #{line.strip.split.last}"
		end
	end
	matrix = "#{matrix}\n"
end

open "#{subdir}/clusterprob.m", "w" do |fout|
	fout.puts "RowLabels = {#{row_labels}}"
	fout.puts "Matrix = [#{matrix}]"
	fout.puts "CG = clustergram(Matrix, 'RowLabels', RowLabels, 'Symmetric', false, 'ColorMap', colormap('JET'))"
	fout.puts "plot(CG)"
end