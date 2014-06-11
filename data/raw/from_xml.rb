require 'active_support/core_ext'
ARGV.each do |filename|
	f = open filename
	xml = Hash.from_xml f.read
	f.close
	ids = {}
	xml["pathway"]["entry"].each do |e|
		if e["name"] == "undefined"
			ids[e["id"]] = e["component"].map{|c| ids[c["id"]]}.join "|"
		elsif e["name"] =~ /hsa/
			ids[e["id"]] = e["name"].split.join "|"
		else
			raise "Undefined entry name: #{e["name"]}"
		end
	end
	open "#{filename.split(".xml").first}_hsa.txt", "w" do |f|
		xml["pathway"]["relation"].each do |r|
			f.puts "#{ids[r["entry1"]]} #{ids[r["entry2"]]}"
		end
	end
end