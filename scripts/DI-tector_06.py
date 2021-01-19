# -*- coding: utf-8 -*-
# ====================================================================================================
# DI-tector tool written by Beauclair Guillaume (2017-06-14)
# Feel free to contact me at guillaume.beauclair@pasteur.fr for questions or to improve the script
# ====================================================================================================

#=================================
import os
import os.path
import sys
import time
import argparse
import gzip
import re
import math
from subprocess import Popen, PIPE, check_output, call, getoutput
import csv

#=================================
Lib_Host = ''
Lib_virus = ''
Fastq_File = ''
FileTag = ''
DeDup = ''
Output_Dir = ''
RED   = "\033[0;31m"  
GREEN = "\033[0;32m"
RESET = "\033[0;0m"
algo_start = time.time()
nb_reads = 0
input_nb_reads = 0
error_txt = ""
too_short_reads = 0

#=================================
## Import arguments from command line
#=================================
if __name__ =='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-g", "--Host_Ref", help="Host genome reference sequence in FASTA format.")
	parser.add_argument("Virus_Ref", help="Virus genome reference sequence in FASTA format.")
	parser.add_argument("Input_Data", help="File containing single reads in FASTQ format.")
	parser.add_argument("-s", "--Min_Segment", help="Minimum segment length. Default is 15.", type=int, default=15)
	parser.add_argument("-m", "--Min_MAPQ", help="Skip alignments with MAPQ smaller than INT. Default is 25.", type=int, default=25)
	parser.add_argument("-n", "--Nb_Reads", help="Show only DVGs with counts reads > or egal to INT. Default is 2.))", type=int, default=2)
	parser.add_argument("-o", "--Output_Directory", help="Enter a directory name that all compiled output files will be saved in.")
	parser.add_argument("-t", "--Tag", help="Enter a tag name that will be appended before each output file. Default is 'DI-tector'.", default="DI-tector")
	#parser.add_argument("-x", "--Index_Ref", action='store_true', help="Index Host and Viral reference sequence. Default is (OFF).")
	#parser.add_argument("-u", "--DeDup", action='store_true', help="Remove potential PCR duplicates. Default is (OFF).")
	parser.add_argument("-d", "--DVG_sequences", action='store_true', help="Generate multi-fasta file with DVG sequences. Default is (OFF).")
	parser.add_argument("-l", "--InDel_Length", help="Skip alignments with size of InDels smaller or egal to INT. Default size is 1.", type=int, default=1)
	parser.add_argument("-f", "--Fasta", action='store_true', help="Select '-f' if data is in FASTA format fasta. Default is FASTQ.")
	parser.add_argument("-p", "--Polarity", help="[0] Positive strand genome / [1] Negative strand genome. Default is 0.", type=int, default=0)
	parser.add_argument("-q", "--No_Quantification", action='store_true', help="Inactive percentage quantification. Quantification need bedtools Default is (ON).")
	parser.add_argument("-k", "--Keep_files", action="store_true", help="Keep intermediaite files (i.e. alignment etc...). Default is (OFF).")
	parser.add_argument("-x", "--Nb_threads", help="Number of threads. Default is 1.", type=int, default=1)
	args = parser.parse_args()
	
	file_tag = args.Tag
	if args.Output_Directory:	  
		out_dir = args.Output_Directory
		if not out_dir.endswith("/"):
		    out_dir = out_dir+'/'
		os.makedirs(out_dir, exist_ok=True)
	else:
		out_dir = ""	
#=================================

#=================================
## Function
#=================================
def reverse_complement(dna):
	complement = {'A': 't', 'C': 'g', 'G': 'c', 'T': 'a', 'N': 'n'}
	return ''.join([complement[base] for base in dna[::-1]])

def draw_progress_bar(text_bar, percent, start, barLen=20):
	sys.stdout.write("\r")
	progress = ""
	for i in range(barLen):
		if i < int(barLen * percent):
			progress += "="
		else:
			progress += " "
	elapsedTime = time.time() - start;
	estimatedRemaining = int(elapsedTime * (1.0/percent) - elapsedTime)

	if (percent == 1.0):
		sys.stdout.write("%s[%s] %.1f%% Elapsed: %im %02is ETA: Done!\n" % 
			(text_bar, progress, percent * 100, int(elapsedTime)/60, int(elapsedTime)%60))
		sys.stdout.flush()
		return
	else:
		sys.stdout.write("%s[%s] %.1f%% Elapsed: %im %02is ETA: %im%02is " % 
			(text_bar, progress, percent * 100, int(elapsedTime)/60, int(elapsedTime)%60,
			estimatedRemaining/60, estimatedRemaining%60))
		sys.stdout.flush()
		return
#=================================

#=================================
## Print info
#=================================
sys.stdout.write(GREEN)
print ('=================================')
print ('Program: DI-tector')
print ('Version: 0.6 - Last modified 25/05/2018')
print ('=================================')
sys.stdout.write(RED)
print ('Requirement: (must be in your $PATH)\n-bwa\n-samtools\n')
print ('Optional: (must be in your $PATH)\n-bedtools\n')
sys.stdout.write(RESET)

print ("Input file:               {}".format(args.Input_Data))
print ("Host reference:	          {}".format(args.Host_Ref))
print ("Virus reference:          {}".format(args.Virus_Ref))
print ("Remove segment:	        < {} nt".format(args.Min_Segment))
print ("Remove reads with MAPQ: < {}".format(args.Min_MAPQ))
print ("Allow InDel. length:    > {} nt".format(args.InDel_Length))
print ("\n", end="")

no_step = 0
total_step = 5
if args.Host_Ref:
	total_step += 1

if args.Fasta:
	factor_nb_reads = 2
else:
	factor_nb_reads = 4

if os.path.splitext(args.Input_Data)[1]==".gz":
	print ("=================================\nReading gziped input file...\n=================================")
	open_file = gzip.open(args.Input_Data, 'rb')			
else:
	print ("=================================\nReading input file...\n=================================")
	open_file = open(args.Input_Data)

no_line = 0
for line in open_file:
	no_line+=1
	if no_line%(factor_nb_reads*100000) == 0:
		sys.stdout.write("\r")
		sys.stdout.write("Number of reads in input file: {:,.0f}".format(no_line/factor_nb_reads))
		sys.stdout.flush()
sys.stdout.write("\r")
sys.stdout.write("Number of reads in input file: {:,.0f}\n".format(no_line/factor_nb_reads))
sys.stdout.flush()
open_file.close()
print ("=================================\n")
input_nb_reads = no_line/factor_nb_reads

if args.Host_Ref:
	no_step += 1
	print("=================================\nStep %i/%i : Alignment against Host reference...Filtering...\n================================="%(no_step,total_step))
	timer_start = time.time()
	with open(out_dir+file_tag+'_temp_file_woHost'+'.fq', 'w') as out_file_woHost:
		p1 = Popen(["bwa", "mem", "-t", str(args.Nb_threads), args.Host_Ref, args.Input_Data], stdout=PIPE)
		p2 = Popen(["samtools", "view", "-S", "-f", "4", "-"], stdin=p1.stdout, stdout=PIPE)
		p1.stdout.close()
		p3 = Popen(["samtools", "fastq", "-"], stdin=p2.stdout, stdout=out_file_woHost) 
		p2.stdout.close()
		p3.communicate()[0]
	print("=================================\nStep {}/{} : Done! ({:.0f}m {:.0f}s)\n=================================\n".format(no_step,total_step,int(time.time()-timer_start)/60,int(time.time()-timer_start)%60))  
	temp_file_woHost = out_dir+file_tag+'_temp_file_woHost'+'.fq'
else:
	temp_file_woHost = args.Input_Data

no_step += 1
print("=================================\nStep %i/%i : Alignment against Viral reference...Filtering...\n================================="%(no_step,total_step))
timer_start = time.time()
if args.No_Quantification:
	with open(out_dir+file_tag+'_temp_file_woVirus'+'.sam', 'w') as out_file_woVirus:
		p1 = Popen(["bwa", "mem", "-t", str(args.Nb_threads), args.Virus_Ref, temp_file_woHost], stdout=PIPE)
		p2 = Popen(["awk", r'$6 ~ /S|H/ || $2 == 4'], stdin=p1.stdout, stdout=out_file_woVirus)
		p1.stdout.close()
		p2.communicate()[0]
else:
	with open(out_dir+file_tag+'_temp_file_onVirus'+'.sam', 'w') as out_file_onVirus:
		Popen(["bwa", "mem", "-t", str(args.Nb_threads), args.Virus_Ref, temp_file_woHost], stdout=out_file_onVirus).communicate()[0]
	with open(out_dir+file_tag+'_temp_file_onVirus'+'.sam', 'r') as out_file_onVirus, open(out_dir+file_tag+'_temp_file_woVirus'+'.sam', 'w') as out_file_woVirus:
		Popen(["awk", r'$6 ~ /S|H/ || $2 == 4'], stdin=out_file_onVirus, stdout=out_file_woVirus).communicate()[0]
	with open(out_dir+file_tag+'_temp_file_onVirus'+'.sam', 'r') as out_file_onVirus, open(out_dir+file_tag+'_temp_file_Virus_unsorted'+'.bam', 'w') as out_file_Virus:
		Popen(["samtools", "view", "-Sb", "-F", "4"], stdin=out_file_onVirus, stdout=out_file_Virus).communicate()[0]
	with open(out_dir+file_tag+'_temp_file_Virus_unsorted'+'.bam', 'r') as out_file_Virus_unsorted, open(out_dir+file_tag+'_temp_file_Virus'+'.bam', 'w') as out_file_Virus:
		Popen(["samtools", "sort"], stdin=out_file_Virus_unsorted, stdout=out_file_Virus).communicate()[0]
	try:
		os.remove(out_dir+file_tag+'_temp_file_Virus_unsorted'+'.bam')
	except OSError:
			pass
	file_bam = out_dir+file_tag+'_temp_file_Virus'+'.bam'
	with open(out_dir+file_tag+'_temp_file_Virus_Coverage'+'.txt', 'w') as out_file_Virus_Cov:
		Popen(["samtools", "depth", "-a", file_bam], stdout=out_file_Virus_Cov).communicate()[0]
	with open(out_dir+file_tag+'_temp_file_Virus_Coverage'+'.txt', 'r') as out_file_Virus_Cov:
		rows = ( line.split('\t') for line in out_file_Virus_Cov )
		d = { row[1]:row[2].replace("\n", "") for row in rows }

no_line = 0
for line in open(out_dir+file_tag+'_temp_file_woVirus'+'.sam', 'r'):
	no_line+=1
nb_reads_after_filtering = no_line
print("=================================\nStep {}/{} : Done! {:,.0f} reads left after filtering ({:.2%}) ({:.0f}m {:.0f}s)\n=================================\n".format(no_step,total_step,nb_reads_after_filtering,nb_reads_after_filtering/input_nb_reads,int(time.time()-timer_start)/60,int(time.time()-timer_start)%60))
if nb_reads_after_filtering == 0:
	print("=================================\nAnalyse is stoped as no reads left after filtering\n=================================\n")
	sys.exit()

no_step += 1
print("=================================\nStep %i/%i : Segmentation...\n================================="%(no_step,total_step))
no_line = 0
timer_start = time.time()
nb_seg = 0

max_SEQ = 0
with open(out_dir+file_tag+'_temp_seqment'+'.fq', 'w') as file_segment:
		for line in open(out_dir+file_tag+'_temp_file_woVirus'+'.sam', 'r'):
			no_line+=1
			splited_line = line.split('\t')
			QNAME = splited_line[0]
			SEQ = splited_line[9]
			len_SEQ = len(SEQ)
			if max_SEQ < len_SEQ:
				max_SEQ = len_SEQ
			if len_SEQ>=2*args.Min_Segment:
				if args.Fasta:
					for no_seg in range(args.Min_Segment, len_SEQ-args.Min_Segment+1):
						file_segment.write("@SEGF"+str(no_seg).zfill(3)+":"+str(len_SEQ).zfill(3)+"::"+QNAME+"\n"+SEQ[:no_seg]+"\n")
						file_segment.write("@SEGL"+str(no_seg).zfill(3)+":"+str(len_SEQ).zfill(3)+"::"+QNAME+"\n"+SEQ[no_seg:]+"\n")
						nb_seg += 2
				else:
					for no_seg in range(args.Min_Segment, len_SEQ-args.Min_Segment+1):
						file_segment.write("@SEGF"+str(no_seg).zfill(3)+":"+str(len_SEQ).zfill(3)+"::"+QNAME+"\n"+SEQ[:no_seg]+"\n+\n"+splited_line[10][:no_seg]+"\n")
						file_segment.write("@SEGL"+str(no_seg).zfill(3)+":"+str(len_SEQ).zfill(3)+"::"+QNAME+"\n"+SEQ[no_seg:]+"\n+\n"+splited_line[10][no_seg:]+"\n")
						nb_seg += 2
			else: 
				error_txt += "Read named {} is too short ({} nt)\n".format(QNAME,len(SEQ))
				too_short_reads += 1
			draw_progress_bar("Segmentation: ", no_line/nb_reads_after_filtering, timer_start)
print("=================================\nStep {}/{} : Done! {:,} reads ({:.2%}) were excluded (length < {} nt) ({:.0f}m {:.0f}s)\n=================================\n".format(no_step,total_step,too_short_reads,too_short_reads/nb_reads_after_filtering,args.Min_Segment*2,int(time.time()-timer_start)/60,int(time.time()-timer_start)%60))

no_step += 1
print("=================================\nStep {}/{} : Second alignment against Viral reference...Filtering...\n=================================".format(no_step,total_step))
timer_start = time.time()

with open(out_dir+file_tag+'_temp_aln'+'.sam', 'w') as file_aln:
	p1 = Popen(["bwa", "aln", "-t", str(args.Nb_threads), args.Virus_Ref, out_dir+file_tag+'_temp_seqment'+'.fq'], stdout=PIPE)
	p2 = Popen(["bwa", "samse", args.Virus_Ref, "-", out_dir+file_tag+'_temp_seqment'+'.fq'], stdin=p1.stdout, stdout=file_aln)
	p1.stdout.close()
	p2.communicate()[0]
print("=================================\nStep {}/{} : Done! ({:.0f}m {:.0f}s)\n=================================\n".format(no_step,total_step,int(time.time()-timer_start)/60,int(time.time()-timer_start)%60)) 

no_step += 1
print("=================================\nStep {}/{} : Filtering...Analyses...\n=================================".format(no_step,total_step))
timer_start = time.time()
sam_headers = ['@HD', '@SQ', '@RG', '@PG', '@CO']
no_line = 0
error_msg = ""
error_ali = ""

prev_DI_type = "DVG's type"
prev_QNAME_F = "*************QNAME_F"
prev_RNAME_F = "RNAME_F"
prev_RNAME_L = "RNAME_L"
prev_BP_POS = "BP_Pos"
prev_RI_POS = "RI_Pos"
prev_MAPQ_F = "MAPQ_F"
prev_MAPQ_L = "MAPQ_L"
prev_CIGAR_F = "CIGAR_F"
prev_CIGAR_L = "CIGAR_L"
prev_MD_CIGAR_F = "MD_CIGAR_F"
prev_MD_CIGAR_L = "MD_CIGAR_L"
prev_POS_F = "POS_F"
prev_POS_L = "POS_L"
prev_SEQ_FL_ori = "SEQ_FL_ori"
prev_SEQ_FL_ref = "SEQ_FL_ref"
prev_size_DI = "Length"
prev_delta_pos = "Delta_Positions"
prev_pos_seg = "Segmentation"

fasta_viral_seq = ""
with open(args.Virus_Ref, 'r') as fasta_viral_file:
	for line in fasta_viral_file:
		if ">" in line:
			id_viral_seq = line[1:]
			continue
		else:
			fasta_viral_seq += line.replace("\n", "")
			continue

nb_seq_tot=0
multi_fasta_seq={}
with open(args.Virus_Ref, 'r') as f:
	for l in f:
		if l.startswith(">"):
			if nb_seq_tot>0:
				multi_fasta_seq.update({name:multi_fasta_seq[name]})
			name=l[1:].strip().split(' ')[0]
			multi_fasta_seq.update({name:""})
			nb_seq_tot+=1
		else:
			multi_fasta_seq.update({name:multi_fasta_seq[name]+l.strip().upper()})
	if nb_seq_tot>0:
		multi_fasta_seq.update({name:multi_fasta_seq[name]})

nb_ins_reads=0
nb_del_reads=0
nb_3cb_reads=0
nb_5cb_reads=0
nb_err=0
nb_seg_analysed=0

with open(out_dir+file_tag+'_temp_aln'+'.sam', 'r') as file_aln, open(out_dir+file_tag+'_output'+'.txt', 'w') as file_output:
	DI_type = "Start"
	size_cut_off = "Too small"
	for line in file_aln:
		if line[:3] in sam_headers:
			continue
		if (str(line[:4]) == "SEGF"):
			nb_seg_analysed += 2
			splited_line_F = line.split('\t')
			if int(splited_line_F[1]) != 4 and int(splited_line_F[4])>=args.Min_MAPQ:
				splited_line_L = next(file_aln).split('\t')
				if int(splited_line_L[1]) != 4 and int(splited_line_L[4])>=args.Min_MAPQ: 
					DI_type = False
					size_cut_off = "To_define"
					SEQ_FL_ref = ""
					SEQ_FL_ori = ""
					insert_jt_nt_F = 0
					insert_jt_str_F = ""
					insert_jt_nt_L = 0
					insert_jt_str_L = ""
					nb_I_F = 0
					nb_I_L = 0

					QNAME_F = splited_line_F[0]
					FLAG_F = int(splited_line_F[1])
					RNAME_F = splited_line_F[2]
					POS_F = int(splited_line_F[3])
					MAPQ_F = splited_line_F[4]
					CIGAR_F = splited_line_F[5]
					SEQ_F = splited_line_F[9]
					len_SEQ_F = len(SEQ_F)
					MD_CIGAR_F = splited_line_F[18].replace("\n", "")[5:]

					QNAME_L = splited_line_L[0]
					FLAG_L = int(splited_line_L[1])
					RNAME_L = splited_line_L[2]
					POS_L = int(splited_line_L[3])
					MAPQ_L = splited_line_L[4]
					CIGAR_L = splited_line_L[5]
					SEQ_L = splited_line_L[9]
					len_SEQ_L = len(SEQ_L)
					MD_CIGAR_L = splited_line_L[18].replace("\n", "")[5:]

					#==========
					# !!! Insertion/Deletion (Forward strand)
					#==========
					if FLAG_F == 0 and FLAG_L == 0:
						#Pour F
						if not CIGAR_F[:-1].isnumeric():
							if 'I' in CIGAR_F:
								nb_I_F = sum(list(map(int, re.findall(r'(\d+)I',CIGAR_F))))
						if not MD_CIGAR_F.isnumeric():
							MD_temp_F = MD_CIGAR_F
							while MD_temp_F[-1] == str(0) and MD_temp_F[-2:-1].isalpha():
								insert_jt_nt_F += 1
								insert_jt_str_F += MD_temp_F[-2:-1]
								MD_temp_F = MD_temp_F[:-2]
						#Pour L
						if not MD_CIGAR_L.isnumeric():
							MD_temp_L = MD_CIGAR_L
							while (MD_temp_L[0] == str(0) and MD_temp_L[1:2].isalpha()):
								insert_jt_nt_L += 1
								insert_jt_str_L += MD_temp_L[1:2]
								MD_temp_L = MD_temp_L[2:]
						BP_POS = POS_F+len_SEQ_F-1-insert_jt_nt_F+insert_jt_nt_L-nb_I_F
						RI_POS = POS_L-insert_jt_nt_F+insert_jt_nt_L
						size_DI = BP_POS+len(fasta_viral_seq)+1-RI_POS
						delta_pos = int(math.fabs(RI_POS-BP_POS-1))
						SEQ_FL_ori = SEQ_F.upper()+SEQ_L.lower()
						SEQ_FL_ref = "aaa"
						
						if delta_pos <= args.InDel_Length:
							size_cut_off = "Too small"
						else:
							size_cut_off = "Ok"

						if BP_POS+1 < RI_POS:
							DI_type = "Deletion DVG (Fwd. strand)"
						elif BP_POS+1 > RI_POS:
							DI_type = "Insertion DVG (Fwd. strand)"
						elif BP_POS+1 == RI_POS:
							error_ali += "Correctly align segments:\t{}\n".format(QNAME_F)
							DI_type = "Correctly align"
						else:
							error_msg += "ERROR POS for read:\t{}\n".format(QNAME_F)
							DI_type = "Error FLAG or POS"
							nb_err += 1
			
					#==========
					# !!! Insertion/Deletion (Reverse strand)
					#==========
					elif FLAG_F == 16 and FLAG_L == 16:
						#Pour F
						if not MD_CIGAR_F.isnumeric(): 
							MD_temp_F = MD_CIGAR_F
							while MD_temp_F[0] == str(0) and MD_temp_F[1:2].isalpha():
								insert_jt_nt_F += 1
								insert_jt_str_F += MD_temp_F[1:2]
								MD_temp_F = MD_temp_F[2:]
						#Pour L
						if not CIGAR_L[:-1].isnumeric():
							if 'I' in CIGAR_L:
								nb_I_L = sum(list(map(int, re.findall(r'(\d+)I',CIGAR_L))))
						if not MD_CIGAR_L.isnumeric():
							MD_temp_L = MD_CIGAR_L
							while MD_temp_L[-1] == str(0) and MD_temp_L[-2:-1].isalpha():
								insert_jt_nt_L += 1
								insert_jt_str_L += MD_temp_L[-2:-1]
								MD_temp_L = MD_temp_L[:-2]

						BP_POS = POS_F+insert_jt_nt_F-insert_jt_nt_L
						RI_POS = POS_L+len_SEQ_L-1+insert_jt_nt_F-insert_jt_nt_L-nb_I_L
						size_DI = len(fasta_viral_seq)+1-BP_POS+RI_POS
						delta_pos = int(math.fabs(RI_POS-BP_POS+1))
						SEQ_FL_ori = reverse_complement(SEQ_F).upper()+reverse_complement(SEQ_L).lower()
						SEQ_FL_ref = "aaa"

						if delta_pos <= args.InDel_Length:
							size_cut_off = "Too small"
						else:
							size_cut_off = "Ok"

						if BP_POS-1 > RI_POS:
							DI_type = "Deletion DVG (Rev. strand)"
						elif BP_POS-1 < RI_POS:
							DI_type = "Insertion DVG (Rev. strand)"
						elif BP_POS-1 == RI_POS:
							error_ali += "Correctly align read:\t{}\n".format(QNAME_F)
							DI_type = "Correctly align"
						else:
							error_msg += "ERROR POS for read:\t{}\n".format(QNAME_F)
							DI_type = "Error FLAG or POS"
							nb_err += 1 

					#==========
					# 5' cb/sb DVG
					#==========
					elif (FLAG_F == 0 and FLAG_L == 16): 
						#Pour F
						if not CIGAR_F[:-1].isnumeric(): 
							if 'I' in CIGAR_F:
								nb_I_F = sum(list(map(int, re.findall(r'(\d+)I',CIGAR_F))))
						if not MD_CIGAR_F.isnumeric():
							MD_temp_F = MD_CIGAR_F
							while MD_temp_F[-1] == str(0) and MD_temp_F[-2:-1].isalpha(): 
								insert_jt_nt_F += 1
								insert_jt_str_F += MD_temp_F[-2:-1]
								MD_temp_F = MD_temp_F[:-2]
						#Pour L
						if not CIGAR_L[:-1].isnumeric(): 
							if 'I' in CIGAR_L:
								nb_I_L = sum(list(map(int, re.findall(r'(\d+)I',CIGAR_L))))
						if not MD_CIGAR_L.isnumeric():
							MD_temp_L = MD_CIGAR_L
							while MD_temp_L[-1] == str(0) and MD_temp_L[-2:-1].isalpha(): 
								insert_jt_nt_L += 1
								insert_jt_str_L += MD_temp_L[-2:-1]
								MD_temp_L = MD_temp_L[:-2]

						BP_POS = POS_F+len_SEQ_F-1-insert_jt_nt_F+insert_jt_nt_L-nb_I_F
						RI_POS = POS_L+len_SEQ_L-1+insert_jt_nt_F-insert_jt_nt_L-nb_I_L
						size_DI = BP_POS+RI_POS
						delta_pos = int(math.fabs(RI_POS-BP_POS))
						size_cut_off = "NA"
						SEQ_FL_ori = SEQ_F.upper()+reverse_complement(SEQ_L).lower()
						SEQ_FL_ref = "aaa"
						if args.Polarity==1: 
							DI_type = "3' cb/sb DVG"
						else:
							DI_type = "5' cb/sb DVG"
					
					#==========
					# 3' cb/sb DVG
					#==========
					elif (FLAG_F == 16 and FLAG_L == 0): #cbDVG and sbDVG
						#Pour F
						if not MD_CIGAR_F.isnumeric(): 
							MD_temp_F = MD_CIGAR_F
							while MD_temp_F[0] == str(0) and MD_temp_F[1:2].isalpha():
								insert_jt_nt_F += 1
								insert_jt_str_F += MD_temp_F[1:2]
								MD_temp_F = MD_temp_F[2:]
						#Pour L
						if not MD_CIGAR_L.isnumeric():
							MD_temp_L = MD_CIGAR_L
							while MD_temp_L[0] == str(0) and MD_temp_L[1:2].isalpha(): 
								insert_jt_nt_L += 1
								insert_jt_str_L += MD_temp_L[1:2]
								MD_temp_L = MD_temp_L[2:]
								
						BP_POS = POS_F+insert_jt_nt_F-insert_jt_nt_L
						RI_POS = POS_L-insert_jt_nt_F+insert_jt_nt_L
						size_DI = 2*(len(fasta_viral_seq)+1)-BP_POS-RI_POS
						delta_pos = int(math.fabs(RI_POS-BP_POS))
						size_cut_off = "NA"
						SEQ_FL_ori = reverse_complement(SEQ_F).upper()+SEQ_L.lower()
						SEQ_FL_ref = "aaa"
						if args.Polarity==1: 
							DI_type = "5' cb/sb DVG"
						else:
							DI_type = "3' cb/sb DVG"
					
					#==========
					# Error
					#==========
					else: 
						error_msg += "ERROR FLAG for read:\t{}\n".format(QNAME_F)
						DI_type = "Error FLAG or POS"
						nb_err +=1


					#==========
					# Writing in file
					#==========
					if QNAME_F[7:]!=prev_QNAME_F[7:] and DI_type!="Correctly align" and DI_type!="Error FLAG or POS" and size_cut_off!="Too small": 
						if DI_type=="3' cb/sb DVG":
							nb_3cb_reads += 1
						elif DI_type=="5' cb/sb DVG":
							nb_5cb_reads += 1
						elif DI_type=="Deletion DVG (Fwd. strand)" or DI_type=="Deletion DVG (Rev. strand)":
							nb_del_reads += 1
						elif DI_type=="Insertion DVG (Fwd. strand)" or DI_type=="Insertion DVG (Rev. strand)":
							nb_ins_reads += 1
						#Writing in the file
						file_output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(prev_DI_type,prev_size_DI,prev_BP_POS,prev_RI_POS,prev_delta_pos,prev_pos_seg,prev_MAPQ_F,prev_MAPQ_L,prev_RNAME_F,prev_RNAME_L,prev_CIGAR_F,prev_CIGAR_L,prev_MD_CIGAR_F,prev_MD_CIGAR_L,prev_POS_F,prev_POS_L,prev_QNAME_F[13:],prev_SEQ_FL_ori))
						prev_DI_type = DI_type
						prev_QNAME_F = QNAME_F
						prev_RNAME_F = RNAME_F
						prev_RNAME_L = RNAME_L
						prev_BP_POS = BP_POS
						prev_RI_POS = RI_POS
						prev_MAPQ_F = MAPQ_F
						prev_MAPQ_L = MAPQ_L
						prev_CIGAR_F = CIGAR_F
						prev_CIGAR_L = CIGAR_L
						prev_MD_CIGAR_F = MD_CIGAR_F
						prev_MD_CIGAR_L = MD_CIGAR_L
						prev_POS_F = POS_F
						prev_POS_L = POS_L
						prev_SEQ_FL_ori = SEQ_FL_ori
						prev_SEQ_FL_ref = SEQ_FL_ref
						prev_size_DI = size_DI
						prev_delta_pos = delta_pos
						prev_pos_seg = str(int(QNAME_F[4:7]))+"|"+str(int(QNAME_F[8:11])-int(QNAME_F[4:7]))
					elif DI_type!="Correctly align" and DI_type!="Error FLAG or POS" and size_cut_off!="Too small" and QNAME_F[7:]==prev_QNAME_F[7:] and ((MD_CIGAR_F.isnumeric() and MD_CIGAR_F > prev_MD_CIGAR_F) or (not MD_CIGAR_F.isnumeric() and not prev_MD_CIGAR_F.isnumeric() and MD_CIGAR_L.isnumeric() and MD_CIGAR_L > prev_MD_CIGAR_L)):
						prev_DI_type = DI_type
						prev_QNAME_F = QNAME_F
						prev_RNAME_F = RNAME_F
						prev_RNAME_L = RNAME_L
						prev_BP_POS = BP_POS
						prev_RI_POS = RI_POS
						prev_MAPQ_F = MAPQ_F
						prev_MAPQ_L = MAPQ_L
						prev_CIGAR_F = CIGAR_F
						prev_CIGAR_L = CIGAR_L
						prev_MD_CIGAR_F = MD_CIGAR_F
						prev_MD_CIGAR_L = MD_CIGAR_L
						prev_POS_F = POS_F
						prev_POS_L = POS_L
						prev_SEQ_FL_ori = SEQ_FL_ori
						prev_SEQ_FL_ref = SEQ_FL_ref
						prev_size_DI = size_DI
						prev_delta_pos = delta_pos
						prev_pos_seg = str(int(QNAME_F[4:7]))+"|"+str(int(QNAME_F[8:11])-int(QNAME_F[4:7]))

			if nb_seg_analysed%50000 == 0:
				draw_progress_bar("Analyse: ", nb_seg_analysed/nb_seg, timer_start)
	draw_progress_bar("Analyse: ", 1, timer_start)
	print("=================================\n")
		
	if DI_type == "Start":
		file_output.write("No Data left after filtering...")
	else:
		file_output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(prev_DI_type,prev_size_DI,prev_BP_POS,prev_RI_POS,prev_delta_pos,prev_pos_seg,prev_MAPQ_F,prev_MAPQ_L,prev_RNAME_F,prev_RNAME_L,prev_CIGAR_F,prev_CIGAR_L,prev_MD_CIGAR_F,prev_MD_CIGAR_L,prev_POS_F,prev_POS_L,prev_QNAME_F[13:],prev_SEQ_FL_ori))
if error_msg != "":
	with open(out_dir+file_tag+'_Errors'+'.txt', 'w') as file_error:
		file_error.write(error_msg)
if error_ali != "":
	with open(out_dir+file_tag+'_Ali'+'.txt', 'w') as file_error:
		file_error.write(error_ali)

print("5' cb/sb DVGs       : {} reads\n3' cb/sb DVGs       : {} reads\nDVGs with deletion  : {} reads\nDVGs with insertion : {} reads\nNb. Errors          : {}".format(nb_5cb_reads,nb_3cb_reads,nb_del_reads,nb_ins_reads,nb_err))
print("=================================\nStep {}/{} : Done! ({:.0f}m {:.0f}s)\n=================================\n".format(no_step,total_step,int(time.time()-timer_start)/60,int(time.time()-timer_start)%60))
with open(out_dir+file_tag+'_summary'+'.txt', 'w') as file_summary:
	file_summary.write("*** Parameters used: ***\n\
-g --Host_Ref         : {}\n\
     Virus_Ref        : {}\n\
     Input_Data       : {}\n\
-s --Min_Segment      : {}\n\
-m --Min_MAPQ         : {}\n\
-n --Nb_Reads         : {}\n\
-o --Output_Directory : {}\n\
-t --Tag              : {}\n\
-d --DVG_sequences    : {}\n\
-l --InDel_Length     : {}\n\
-f --Fasta            : {}\n\
-p --Polarity         : {}\n\
-k --Keep_files       : {}\n\
-x --Nb_threads       : {}\n\
".format(args.Host_Ref,\
args.Virus_Ref,\
args.Input_Data,\
args.Min_Segment,\
args.Min_MAPQ,\
args.Nb_Reads,\
out_dir,\
args.Tag,\
args.DVG_sequences,\
args.InDel_Length,\
args.Fasta,\
args.Polarity,\
args.Keep_files,\
args.Nb_threads))
	file_summary.write("\n*** Results ***\n\
Input                 : {:,.0f} reads\n\
After filtering       : {:,.0f} reads\n\
\n\
5' cb/sb DVGs         : {:,.0f} reads\n\
3' cb/sb DVGs         : {:,.0f} reads\n\
DVGs with deletion    : {:,.0f} reads\n\
DVGs with insertion   : {:,.0f} reads\n\
Errors                : {:,.0f}".format(input_nb_reads,nb_reads_after_filtering,nb_5cb_reads,nb_3cb_reads,nb_del_reads,nb_ins_reads,nb_err))

no_step += 1
print("=================================\nStep {}/{} : Sorting...Export...\n=================================".format(no_step,total_step))
timer_start = time.time()
with open(out_dir+file_tag+'_output_sorted'+'.txt', 'w') as file_sorted:
	p1 = Popen(["head", "-1", out_dir+file_tag+'_output'+'.txt'], stdout=file_sorted)
	p2 = Popen(["tail", "-n", "+2", out_dir+file_tag+'_output'+'.txt'], stdout=PIPE)
	p3 = Popen(["sort", "-t\t", "-k1,1", "-k2,2g", "-k3,3g", "-k4,4g", "-k7,7gr", "-k8,8gr"], stdin=p2.stdout, stdout=file_sorted)
	p2.stdout.close()
	p3.communicate()[0]

prev_0 = "DVG's type"
prev_1 = "Length"
prev_2 = "BP_Pos"
prev_3 = "RI_Pos"
prev_4 = "Delta_Positions"
prev_5 = "Ref"
prev_6 = "Counts"
prev_7 = "%_to_Virus"
prev_none = "None or Reads number below cut-off\n"
	 
with open(out_dir+file_tag+'_output_sorted'+'.txt', 'r') as file_sorted, \
	 open(out_dir+file_tag+'_counts'+'.txt', 'w') as file_counts, \
	 open(out_dir+file_tag+'_fasta'+'.fa', 'w') as file_fasta:
	next(file_sorted) 
	for line in file_sorted:
		splited_line = line.split('\t')
		if splited_line[0]==prev_0 and splited_line[1]==prev_1 and splited_line[2]==prev_2 and splited_line[3]==prev_3 and splited_line[4]==prev_4 and (splited_line[8]+"|"+splited_line[9])==prev_5: 
			prev_6 += 1
		else:
			if prev_6 != "Counts" and prev_6 >= args.Nb_Reads:
				prev_none = ""
				if args.No_Quantification or int(prev_2) > len(multi_fasta_seq[prev_5a]) or int(prev_2) <= 0 or int(prev_3) > len(multi_fasta_seq[prev_5b]) or int(prev_3) <= 0:
					prev_7 = "ND"
				else:
					if int(d[str(prev_2)]) > 0:
						str_per_2 = str(round((prev_6/int(d[str(prev_2)])*100),2))+"%"
					else:
						str_per_2 = "ND"
					if int(d[str(prev_3)]) > 0:
						str_per_3 = str(round((prev_6/int(d[str(prev_3)])*100),2))+"%"
					else:
						str_per_3 = "ND"
					prev_7 = str_per_2+"|"+str_per_3
				file_counts.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(prev_0,prev_1,prev_2,prev_3,prev_4,prev_5,prev_6,prev_7))
				
				if args.DVG_sequences:
					new_seq_DVG=""
					if prev_0=="3' cb/sb DVG":
						new_seq_DVG=(multi_fasta_seq[prev_5a][:int(prev_2)]+reverse_complement(multi_fasta_seq[prev_5b][:int(prev_3)])).upper()
					elif prev_0=="5' cb/sb DVG":
						new_seq_DVG=(reverse_complement(multi_fasta_seq[prev_5a][int(prev_2)-1:])+multi_fasta_seq[prev_5b][int(prev_3)-1:]).upper()
					elif prev_0=="Deletion DVG (Fwd. strand)" or prev_0=="Insertion DVG (Fwd. strand)":
						new_seq_DVG=(multi_fasta_seq[prev_5a][:int(prev_2)]+multi_fasta_seq[prev_5b][int(prev_3)-1:]).upper()
					elif prev_0=="Deletion DVG (Rev. strand)" or prev_0=="Insertion DVG (Rev. strand)":
						new_seq_DVG=(reverse_complement(multi_fasta_seq[prev_5a][int(prev_2)-1:])+reverse_complement(multi_fasta_seq[prev_5b][:int(prev_3)])).upper()
					else:
						new_seq_DVG="Error"
					file_fasta.write(">{}\n{}\n".format(("DVG_"+prev_1+"_BP_"+prev_2+"_"+prev_5a+"_RI_"+prev_3+"_"+prev_5b),new_seq_DVG))
				
			if prev_0!=splited_line[0]:
				if prev_6 != "Counts":
					file_counts.write("{}\n".format(prev_none))
				file_counts.write("=================================\n= {}\n=================================\nDVG's type\tLength\tBP_Pos\tRI_Pos\tDelta_Positions\tRef\tCounts\t%_to_Virus\n".format(splited_line[0]))
				prev_none = "None or Reads number below cut-off\n"
			prev_0 = splited_line[0]
			prev_1 = splited_line[1]
			prev_2 = splited_line[2]
			prev_3 = splited_line[3]
			prev_4 = splited_line[4]
			prev_5 = splited_line[8]+"|"+splited_line[9]
			prev_5a = splited_line[8]
			prev_5b = splited_line[9]
			prev_6 = 1
			if args.No_Quantification or int(prev_2) > len(multi_fasta_seq[prev_5a]) or int(prev_2) <= 0 or int(prev_3) > len(multi_fasta_seq[prev_5b]) or int(prev_3) <= 0:
				prev_7 = "ND"
			else:
				if int(d[str(prev_2)]) > 0:
					str_per_2 = str(round((prev_6/int(d[str(prev_2)])*100),2))+"%"
				else:
					str_per_2 = "ND"
				if int(d[str(prev_3)]) > 0:
					str_per_3 = str(round((prev_6/int(d[str(prev_3)])*100),2))+"%"
				else:
					str_per_3 = "ND"
				prev_7 = str_per_2+"|"+str_per_3
	if prev_6 != "Counts" and prev_6 >= args.Nb_Reads:
		file_counts.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(prev_0,prev_1,prev_2,prev_3,prev_4,prev_5,prev_6,prev_7))
	else:
		file_counts.write("{}".format(prev_none))

	if args.DVG_sequences and prev_6.isnumeric():
		new_seq_DVG=""
		if prev_0=="3' cb/sb DVG":
			new_seq_DVG=(multi_fasta_seq[prev_5a][:int(prev_2)]+reverse_complement(multi_fasta_seq[prev_5b][:int(prev_3)])).upper()
		elif prev_0=="5' cb/sb DVG":
			new_seq_DVG=(reverse_complement(multi_fasta_seq[prev_5a][int(prev_2)-1:])+multi_fasta_seq[prev_5b][int(prev_3)-1:]).upper()
		elif prev_0=="Deletion DVG (Fwd. strand)" or prev_0=="Insertion DVG (Fwd. strand)":
			new_seq_DVG=(multi_fasta_seq[prev_5a][:int(prev_2)]+multi_fasta_seq[prev_5b][int(prev_3)-1:]).upper()
		elif prev_0=="Deletion DVG (Rev. strand)" or prev_0=="Insertion DVG (Rev. strand)":
			new_seq_DVG=(reverse_complement(multi_fasta_seq[prev_5a][int(prev_2)-1:])+reverse_complement(multi_fasta_seq[prev_5b][:int(prev_3)])).upper()
		else:
			new_seq_DVG="Error"
		file_fasta.write(">{}\n{}\n".format(("DVG_"+prev_1+"_BP_"+prev_2+"_"+prev_5a+"_RI_"+prev_3+"_"+prev_5b),new_seq_DVG))

try:
    os.remove(out_dir+file_tag+"_output.txt")
except OSError:
	pass
print("Output files:\n/{}{}_summary.txt\n/{}{}_counts.txt\n/{}{}_output_sorted.txt".format(out_dir,file_tag,out_dir,file_tag,out_dir,file_tag))
if args.Keep_files:
	print("/{}{}_temp_file_woHost.fq\n/{}{}_temp_file_onVirus.sam\n/{}{}_temp_file_woVirus.sam\n/{}{}_temp_file_Virus.bam\n/{}{}_temp_aln.sam\n/{}{}_temp_seqment.fq\n/{}{}_temp_file_Virus_Coverage.txt".format(out_dir,file_tag,out_dir,file_tag,out_dir,file_tag,out_dir,file_tag,out_dir,file_tag,out_dir,file_tag,out_dir,file_tag))
else:
	try:
		os.remove(out_dir+file_tag+"_temp_file_woHost.fq")
	except OSError:
		pass
	try:
		os.remove(out_dir+file_tag+"_temp_file_woVirus.sam")
	except OSError:
		pass
	try:
		os.remove(out_dir+file_tag+"_temp_file_onVirus.sam")
	except OSError:
		pass
	try:
		os.remove(out_dir+file_tag+"_temp_file_Virus.bam")
	except OSError:
		pass
	try:
		os.remove(out_dir+file_tag+"_temp_file_Virus_Coverage.txt")
	except OSError:
		pass
	try:
		os.remove(out_dir+file_tag+"_temp_aln.sam")
	except OSError:
		pass
	try:
		os.remove(out_dir+file_tag+"_temp_seqment.fq")
	except OSError:
		pass
if args.DVG_sequences:
	print("/{}{}_fasta.fa".format(out_dir,file_tag))
else:
	try:
		os.remove(out_dir+file_tag+"_fasta.fa")
	except OSError:
		pass
print("=================================\nStep {}/{} : Done! ({:.0f}m {:.0f}s)\n=================================\n".format(no_step,total_step,int(time.time()-timer_start)/60,int(time.time()-timer_start)%60))
sys.stdout.write(GREEN)
print ('=================================\n{:,.0f} reads analysed in {:.0f}m {:.0f}s\n================================='.format(input_nb_reads, int(time.time()-algo_start)/60, int(time.time()-algo_start)%60))
sys.stdout.write(RESET)