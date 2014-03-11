netfile, outfile = ARGV
lines = []
open netfile do |f|
	until (line = f.gets).nil?
		lines << line.strip
	end
end
open outfile, "w" do |f|
	lines.each{|l| f.puts l.split[0..1].join(" ")}
end