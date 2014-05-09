#This script runs preach_prob on all correlation files in a dataset directory
dataset = ARGV[0]
Dir.glob("#{dataset}/coexp_*.txt").each do |coexp_file|
	puts "../../bin/Release/preach_prob.exe #{dataset}/#{dataset}_hsa.txt #{coexp_file} #{coexp_file.gsub("coexp", "prob").gsub(".txt", ".out")}"
	`../../bin/Release/preach_prob.exe #{dataset}/#{dataset}_hsa.txt #{coexp_file} #{coexp_file.gsub("coexp", "prob").gsub(".txt", ".out")}`
end