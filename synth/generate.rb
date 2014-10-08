edges_per_node, steps, step_size, versions = ARGV.map{|a| a.to_i}

@connected = {}
@edges = []
@no_in = []
@no_out = []
@chances = []

def connect node1, node2
	@connected["#{node1}-#{node2}"] = true
	@connected["#{node2}-#{node1}"] = true
	@chances << node1 << node2
	if rand > 0.5
		@edges << [node1, node2]
		@no_in.delete node2
		@no_out.delete node1
	else
		@edges << [node2, node1]
		@no_in.delete node1
		@no_out.delete node2
	end
end

def output filename
	open "#{filename}.txt", "w" do |f|
		@edges.each{|e| f.puts e.join(" ")}
	end
	open "#{filename}_coexp.txt", "w" do |f|
		@no_in.each do |s|
			@no_out.each do |t|
				f.puts "#{s} #{t} #{rand.round 5}"
			end
		end
		if @no_in.empty? or @no_out.empty?
			s = @edges[(rand*@edges.size).floor].first
			t = @edges[(rand*@edges.size).floor].last
			while (s == t)
				t = @edges[(rand*@edges.size).floor].last
			end
			f.puts "#{s} #{t} #{rand.round 5}"
		end
	end
end

1.upto versions do |version|
	@connected = {}
	@edges = []
	@no_in = []
	@no_out = []
	@chances = []
	steps.times do |step|
		1.upto step_size do |i|
			node1 = step * step_size + i
			@no_out << node1
			@no_in << node1
			if node1 <= edges_per_node + 1
				1.upto node1-1 do |node2|
					connect node1, node2
				end
			else
				edges_per_node.times do
					node2 = @chances[(rand*@chances.size).floor]
					while node2 == node1 or @connected["#{node1}#{node2}"]
						node2 = @chances[(rand*@chances.size).floor]
					end
					connect node1, node2
				end
			end
		end
		output "BA_#{edges_per_node}_#{(step+1)*step_size}_#{version}"
	end
end