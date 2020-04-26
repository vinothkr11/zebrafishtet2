import re

# The program is for quantifying the different types of CRISPR-based ssODN knock-in 
# mutation outcomes

# 0. Define some important sequence variables for matching patterns in the SAM alignment files
files_sgRNAs = {}
wt_file  = "wt_full.sam"
files_sgRNAs["tet2_sgR1.sam"] = {}
files_sgRNAs["tet2_sgR1.sam"]["seq"] = "TTCCTGATACCATCACCTCC"
files_sgRNAs["tet2_sgR1.sam"]["start"] = 50
files_sgRNAs["tet2_sgR1.sam"]["end"] = 150
files_sgRNAs["tet2_sgR2.sam"] = {}
files_sgRNAs["tet2_sgR2.sam"]["seq"] = "CCCATTTGCCAGACAGAACC"
files_sgRNAs["tet2_sgR2.sam"]["start"] = 65
files_sgRNAs["tet2_sgR2.sam"]["end"] = 165
files_sgRNAs["tet2_sgR3.sam"] = {}
files_sgRNAs["tet2_sgR3.sam"]["seq"] = "GGATAGAACCAACCATGTTG"
files_sgRNAs["tet2_sgR3.sam"]["start"] = 25
files_sgRNAs["tet2_sgR3.sam"]["end"] = 125

sgRNAs = {}
sgRNAs['sgR1'] = "TTCCTGATACCATCACCTCC"
sgRNAs['sgR2'] = "CCCATTTGCCAGACAGAACC"
sgRNAs['sgR3'] = "GGATAGAACCAACCATGTTG"



# categories to which to assign the alignments
CRISPR_cats = {}
CRISPR_cats["tet2_sgR1.sam"] = {}
CRISPR_cats["tet2_sgR1.sam"]["wt"] = 0
CRISPR_cats["tet2_sgR1.sam"]["del"] = 0
CRISPR_cats["tet2_sgR1.sam"]["ins"] = 0
CRISPR_cats["tet2_sgR2.sam"] = {}
CRISPR_cats["tet2_sgR2.sam"]["wt"] = 0
CRISPR_cats["tet2_sgR2.sam"]["del"] = 0
CRISPR_cats["tet2_sgR2.sam"]["ins"] = 0
CRISPR_cats["tet2_sgR3.sam"] = {}
CRISPR_cats["tet2_sgR3.sam"]["wt"] = 0
CRISPR_cats["tet2_sgR3.sam"]["del"] = 0
CRISPR_cats["tet2_sgR3.sam"]["ins"] = 0


# open file for writing the coordinates
coord_file = open("first_match_coordinates.csv", "a")
summary_file = open("summary.csv", "a")

coord_file.write('Result' + ',' + 'sgRNA' + ',' + 'FirstCoordinate' + '\n')
summary_file.write('Sample' + ',' + 'Result' + ',' + 'Count' + '\n')

# iterate over all the results after sgRNA CRISPR targeting
for file in files_sgRNAs.keys():
	with open(file) as f:
	
		Counter = 0
		
		# file for writing individual insertion events
		ins_file_name = file + "_insertions.txt"
		ins_file = open(ins_file_name, "a")

		# file for writing individual deletion events
		del_file_name = file + "_deletions.txt"
		del_file = open(del_file_name, "a")
		
		for line in f:
		
			Counter += 1
			
			line_data = line.split("\t")
			
			# assign relevant parts of the data to variables
			CIGAR = line_data[5]
			cigar_matches = re.findall(r'(\d+)([A-Z]{1})', CIGAR)
			seq = line_data[9]			

			# case 1: the 
			if re.search(files_sgRNAs[file]["seq"], seq):
				CRISPR_cats[file]["wt"] += 1;
				
				# write a line to the coordinate file
				if(len(cigar_matches) >= 1):
					coord_file.write('wt' + ',' + file[5:9] + ',' + cigar_matches[0][0] + '\n')
				
				
			else:
					
		# quantify the total net insertions or deletions in the target region 
				pos = 0                 # current location
				netIndel = 0            # net indel length
				
				# go through all CIGAR matches
				for m in cigar_matches:
					pos += int(m[0])
					
					# check that you are still within the 100-bp window around the sgRNA cut site
					if pos >= files_sgRNAs[file]["end"]:
						break
					else:
						if pos > files_sgRNAs[file]["start"] and m[1] == 'I':
							netIndel += int(m[0])
						elif pos > files_sgRNAs[file]["start"] and m[1] == 'D':
							netIndel -= int(m[0])
						
				# classify the corresponding indel events:
				if netIndel > 0:
					CRISPR_cats[file]["ins"] += 1
					
					# writing this read out to a file
					ins_file.write(">" + "ins_" + CIGAR + "-" + str(Counter) + "\n" + seq + "\n")

					# write the coordinate
					coord_file.write('ins' + ',' + file[5:9] + ',' + cigar_matches[0][0] + '\n')
					
				elif netIndel < 0:
					CRISPR_cats[file]["del"] += 1
					# writing this read out to a file
					del_file.write(">" + "del_" + CIGAR + "-" + str(Counter)+ "\n" + seq + "\n")	

					# write the coordinate
					coord_file.write('del' + ',' + file[5:9] + ',' + cigar_matches[0][0] + '\n')
					
		# close all files
		del_file.close()
		ins_file.close()

# iterate over the wild-type sequencing with different sgRNAs

wt_cats = {}
wt_cats["wt_sgR1"] = {}
wt_cats["wt_sgR1"]["wt"] = 0
wt_cats["wt_sgR1"]["del"] = 0
wt_cats["wt_sgR1"]["ins"] = 0
wt_cats["wt_sgR2"] = {}
wt_cats["wt_sgR2"]["wt"] = 0
wt_cats["wt_sgR2"]["del"] = 0
wt_cats["wt_sgR2"]["ins"] = 0
wt_cats["wt_sgR3"] = {}
wt_cats["wt_sgR3"]["wt"] = 0
wt_cats["wt_sgR3"]["del"] = 0
wt_cats["wt_sgR3"]["ins"] = 0

for sgrna in sgRNAs.keys():
	with open(wt_file ) as f:
	
		label = 'wt_' + sgrna
		file = 'tet2_' + sgrna + '.sam'
		
		for line in f:
					
			line_data = line.split("\t")
			
			# assign relevant parts of the data to variables
			CIGAR = line_data[5]
			cigar_matches = re.findall(r'(\d+)([A-Z]{1})', CIGAR)
			seq = line_data[9]			

			# case 1: the 
			if re.search(sgRNAs[sgrna], seq):
				wt_cats[label]["wt"] += 1;
				
				# write a line to the coordinate file
				if(len(cigar_matches) >= 1):
					coord_file.write('wt' + ',' + label + ',' + cigar_matches[0][0] + '\n')
			
			else:
					
		# quantify the total net insertions or deletions in the target region 
				pos = 0                 # current location
				netIndel = 0            # net indel length
				
				# go through all CIGAR matches
				for m in cigar_matches:
					pos += int(m[0])
					
					# check that you are still within the 100-bp window around the sgRNA cut site
					if pos >= files_sgRNAs[file]["end"]:
						break
					else:
						if pos > files_sgRNAs[file]["start"] and m[1] == 'I':
							netIndel += int(m[0])
						elif pos > files_sgRNAs[file]["start"] and m[1] == 'D':
							netIndel -= int(m[0])
						
				# classify the corresponding indel events:
				if netIndel > 0:
					wt_cats[label]["ins"] += 1
					
					# write the coordinate
					coord_file.write('ins' + ',' + label + ',' + cigar_matches[0][0] + '\n')
					
				elif netIndel < 0:
					wt_cats[label]["del"] += 1
	

					# write the coordinate
					coord_file.write('del' + ',' + label + ',' + cigar_matches[0][0] + '\n')


# write the summary results to the file
for file in CRISPR_cats.keys():
	summary_file.write(file[5:9] + ',' + 'wt' + ',' + str(CRISPR_cats[file]["wt"]) + '\n')
	summary_file.write(file[5:9] + ',' + 'ins' + ',' + str(CRISPR_cats[file]["ins"]) + '\n')
	summary_file.write(file[5:9] + ',' + 'del' + ',' + str(CRISPR_cats[file]["del"]) + '\n')
	
for label in wt_cats.keys():
	summary_file.write(label + ',' + 'wt' + ',' + str(wt_cats[label]["wt"]) + '\n')
	summary_file.write(label + ',' + 'ins' + ',' + str(wt_cats[label]["ins"]) + '\n')
	summary_file.write(label + ',' + 'del' + ',' + str(wt_cats[label]["del"]) + '\n')
		
# close main files
coord_file.close()
summary_file.close()