sourcesfile, targetsfile, exprfile = ARGV
sources, targets = [], []
open sourcesfile do |f|
	until (line = f.gets).nil?
		sources << line.strip
	end
end
open targetsfile do |f|
	until (line = f.gets).nil?
		targets << line.strip
	end
end
open exprfile, "w" do |f|
	sources.each{|s| targets.each{|t| f.puts "#{s} #{t} #{rand}"}}
end