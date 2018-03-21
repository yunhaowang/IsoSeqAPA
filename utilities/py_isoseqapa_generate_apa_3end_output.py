#!/usr/bin/env python
import sys,time,re,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	if args.end_bed:
		dic_chr_strand_peak = parse_3end_bed(args.end_bed)
	else:
		dic_chr_strand_peak = {}
	dic_lr_pa = parse_long_read_polya(args.lr_gpd)
	generate_apa_output(args.input,args.output,dic_chr_strand_peak,dic_lr_pa,args.dis_group)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def parse_3end_bed(input_3end_bed):
	dic_chr_strand_peak = {}
	for line in input_3end_bed:
		chr,start,end,peak_name,score,strand = line.strip().split("\t")[:6]
		chr_strand = chr + "&" + strand
		if chr_strand not in dic_chr_strand_peak.keys():
			dic_chr_strand_peak[chr_strand] = {}
			dic_chr_strand_peak[chr_strand]["start"] = []
			dic_chr_strand_peak[chr_strand]["end"] = []
			dic_chr_strand_peak[chr_strand]["start"].append(int(start)+1)
			dic_chr_strand_peak[chr_strand]["end"].append(int(end))
		else:
			dic_chr_strand_peak[chr_strand]["start"].append(int(start)+1)
			dic_chr_strand_peak[chr_strand]["end"].append(int(end))
	input_3end_bed.close()
	return dic_chr_strand_peak



def parse_long_read_polya(input_lr_gpd):
	dic_lr_pa = {}
	for line in input_lr_gpd:
		read_id,chr,strand,tss,tts,mapq,sf_flag = line.strip().split("\t")[1:8]
		if read_id.split("_")[-1].endswith("T1"): # check 3 end completeness
			
			if strand == "+":
				if sf_flag.split("_")[1] == "0": # check soft-clip
					dic_lr_pa[read_id] = tts
			elif strand == "-":
				if sf_flag.split("_")[0] == "0": # check soft-clip
					dic_lr_pa[read_id] = (int(tss)+1)
			else:
				pass
	input_lr_gpd.close()
	return dic_lr_pa

def generate_apa_output(input_con_gpd,output_apa_gpd,dic_chr_strand_peak,dic_lr_pa,gather_dis):
	for line in input_con_gpd:
		chr,strand = line.strip().split("\t")[2:4]
		read_set = line.strip().split("\t")[-1].split(",")
		dic_lr_pa_set = {}
		for read in read_set:
			if read.endswith("T0"): continue # remove the read with truncated 3'end
			if read in dic_lr_pa.keys():
				lr_pa = int(dic_lr_pa[read])
				if lr_pa not in dic_lr_pa_set.keys():
					dic_lr_pa_set[lr_pa] = 1
				else:
					dic_lr_pa_set[lr_pa] += 1
		if dic_lr_pa_set == {}: # No polyA information
			print >>output_apa_gpd, "\t".join(["\t".join(line.strip().split("\t")[:-1]),"NA","NA"])
		else: # Have valid polyA
			pos_list = dic_lr_pa_set.keys()
			pos_list.sort()
			dic_pa_lr_3end = {}
			chr_strand = chr+"&"+strand
			if chr_strand in dic_chr_strand_peak.keys(): # 3end-seq
				for pa in pos_list:
					dic_pa_lr_3end[pa] = str(pa) + "_" + str(dic_lr_pa_set[pa]) + "_0"
					for i in range(0,len(dic_chr_strand_peak[chr_strand]["start"])):
						if pa >= dic_chr_strand_peak[chr_strand]["start"][i] and pa <= dic_chr_strand_peak[chr_strand]["end"][i]:
							dic_pa_lr_3end[pa] = str(pa) + "_" + str(dic_lr_pa_set[pa]) + "_1"
							break
			else:
				for pa in pos_list:
					dic_pa_lr_3end[pa] = str(pa) + "_" + str(dic_lr_pa_set[pa]) + "_0"

			if len(pos_list) == 1: # one pa
				idv_pa_set = dic_pa_lr_3end[pos_list[0]]
				grp_pa_set = idv_pa_set
				print >>output_apa_gpd, "\t".join(["\t".join(line.strip().split("\t")[:-1]),idv_pa_set,grp_pa_set])
			else:

				# individual pa
				idv_pa_list = []
				for idv_pa in pos_list:
					idv_pa_list.append(dic_pa_lr_3end[idv_pa])
				idv_pa_set = ",".join(idv_pa_list)

				# gather pa
				group_idx_set = []
				group_idx = [pos_list[0]]
				for i in range(1,len(pos_list)):
					if (pos_list[i] - gather_dis) > group_idx[-1]:
						group_idx_set.append(group_idx)
						group_idx = [pos_list[i]]
						if pos_list[i] == pos_list[-1]:
							group_idx_set.append(group_idx)
					else:
						group_idx.append(pos_list[i])
						if pos_list[i] == pos_list[-1]:
							group_idx_set.append(group_idx)

				# grouped pa
				grp_pa_list = []
				for gp in group_idx_set:
					sum_pos_frq_lr = 0 # long read data
					sum_pos_frq_sr = 0 # 3end data
					dic_pos_frq_pa = {}
					for pa in gp:
						pos_frq = int(dic_pa_lr_3end[pa].split("_")[1]) + int(dic_pa_lr_3end[pa].split("_")[2])
						sum_pos_frq_lr += int(dic_pa_lr_3end[pa].split("_")[1])
						sum_pos_frq_sr += int(dic_pa_lr_3end[pa].split("_")[2])
						if pos_frq not in dic_pos_frq_pa.keys():
							dic_pos_frq_pa[pos_frq] = []
							dic_pos_frq_pa[pos_frq].append(pa)
						else:
							dic_pos_frq_pa[pos_frq].append(pa)
					max_pos_frq_pa = dic_pos_frq_pa[max(dic_pos_frq_pa.keys())][0]
					grp_pa_list.append(str(max_pos_frq_pa)+"_"+str(sum_pos_frq_lr)+"_"+str(sum_pos_frq_sr))
				grp_pa_set = ",".join(grp_pa_list)
				print >>output_apa_gpd, "\t".join(["\t".join(line.strip().split("\t")[:-1]),idv_pa_set,grp_pa_set])

def do_inputs():
	output_gpd_format = '''
1. gene id
2. isoform id
3. chromosome id
4. strand
5. TSS (+)
6. TTS (+)
7. number of support full-length long reads
8. number of support total long reads
9. exon count
10. exon start set
11. exon end set
12. For novel isoform, derived genic locus
13. For novel isoform, overlap percentage with derived genic locus
14. For novel singleton isoform, if it is located at the last exon of any known isoform. If yes, isoform ID otherwise '-'
15. For novel singleton isoform, the overlap percentage with the the last exon
16. For novel multi-exon isoform, number of splice sites are detected by anno and/or short reads; and the total number of splice sites
17. For novel multi-exon isoform, if the multi-exon isoform is the subset (based on splice junction combination) of known multi-exon isoform, isoform ID if yes otherwise '-'
18. For novel isoform, maximal length of polyA track in defined region
19. For novel isoform, maximal percentage of nucleotide A in defined region
20. Individual polyA site information: polyA site postion; number of supporting long reads; number of supporting 3End-Seq position
21. Grouped polyA site information: polyA site postion; number of supporting long reads; number of supporting 3End-Seq position'''

	parser = argparse.ArgumentParser(description="Function: generate final output file (gpd format)",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: constructed isoform (GPD file) generated by 'py_isoseqapa_generate_output.py'")
	parser.add_argument('-l','--lr_gpd',type=argparse.FileType('r'),required=True,help="Input: polished long read (GPD file) generated by 'py_isoseqcon_polish.py'")
	parser.add_argument('-3','--end_bed',type=argparse.FileType('r'),help="Optional input: 3End-Seq data (BED format)")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: constructed isoform with polyA information (GPD file)")
	parser.add_argument('-d','--dis_group',type=int,default=5,help="If the distance (bp) between two polyA sites is <= the set, group them together")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
