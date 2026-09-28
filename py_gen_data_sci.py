#!/usr/bin/env python3
import sys
import argparse

#this file is a mess of a program I used to complete the Python for Genomic Data Science course from coursera... I hope nobody steals this and cheats -_-'

"""
This program is going to be controlled via command line
Hence, there must be options connected and different functions to be called potentially...

This module/program is for the final exam for the python for genomic data science
course on coursera.

This final exam asks me to do the following:
	1. Input a multi-fasta file
	2. Store all records in the fasta file by their identifiers and nucleotide sequences
	3. Find the lengths of every record via their identifier
	4. Find all longest and shortest records and output their sequences and identifiers
	5. Given a reading frame:
		a) Find out all ORFs for every record in the fasta file (maybe)
		b) Sort and identify which record has the longest ORF and output it'd identifier
	6. Given an identifier:
		a) What is the longest ORF, their reading frame and starting position?
	7. Given an integer 'k':
		a) identify all repeated 'k-mers'
		b) count how often they repeat
		c) sort and identify the most repeated k-mer
"""


def FASTA_dict(fasta_file):
	try:
		fasta=open(fasta_file,'r') #this allows us to read the fasta file
	except IOError:
		print("%s aint here, bud..." %fasta)
	fasta_dictionary={}     #this sets up the fasta dictionary
	for line in fasta:      #this loop will go through all lines in our inp>
		if line[:1]=='>': #this line will make sure that the only lines>
			split_line = line.split(" ",1)
			identifier = split_line[0][1:]
			fasta_dictionary[identifier]=fasta.readline().rstrip()
		else:
			latest_key=list(fasta_dictionary.keys())[-1] #this is the last added key
			fasta_dictionary[latest_key]=fasta_dictionary[latest_key]+line.rstrip() #continues adding to the sequence since every new line not containing '>' is a part of the sequence we did not add yet
	fasta.close() #this closes the input file so we free up memory
	return fasta_dictionary #this returns but does not print our dictionary

"""
the following function takes a fasta dictionary and then measures the length of the sequences
and stores those in a dictionary made using the identifiers in the fasta dictionary:
fasta_dictionary={identifier:sequence} and fasta_length_dictionary={identifier:length}
I used python list comprehension on a dictionary to measure and store all lengths to a dictionary
"""

def FASTA_dict_lengths(fasta_dict):
	fasta_length_dict={}
	fasta_length_dict={identifier: len(fasta_dict[identifier]) for identifier in fasta_dict} #list comprehension!!!
	sorted_fasta_length_dict=dict(sorted(fasta_length_dict.items(), key = lambda item: item[1]))
		#the above line first uses the sorted function (not method), to create an ordered set of tuples based on the values of fasta_length_dict rather than the keys. 
		#This is done using lambda, which essentially tells the sorted function by which
		#key or reference the sorted function will sort.
	return sorted_fasta_length_dict

"""
this function will work on a single sequence. and reading_frame (default is 1) is an optional argument
It will provide all ORFs on the forward strand, given a reading frame (1,2,3) ONLY
It will return a list of all ORFs
"""


def codon_indexer(sequence):
	#start_codon="ATG"
	codons=["ATG","TAA","TAG","TGA"]
	codons_index=[] #this will contain all the codons we search for and their positions (NOT index)
	for codon in codons:	#searches for every codon we are looking for in our list of codons
		moving_index = sequence.find(codon)
		while moving_index > -1:	#will continue looping until i is more than -1 because the .find() method returns a -1 if nothing is found
			codons_index.append([codon,moving_index+1]) #stores position and not index!!
			moving_index = sequence.find(codon,moving_index+1)
			#print(codon + str(i))
	#['ATG 125', 'ATG 176', 'ATG 464', 'TAG 23', 'TAG 77', 'TAG 89', 'TAG 281']
	#print(codons_index)
	ORF_1=[]
	ORF_2=[]
	ORF_3=[]
	for codon_pos_set in codons_index:
		position = codon_pos_set[1]
		if position%3 == 0:
			ORF_3.append(codon_pos_set)
		elif position%3 == 2:
			ORF_2.append(codon_pos_set)
		elif position%3 == 1:
			ORF_1.append(codon_pos_set)
#	print("ORF_1 IS ASDDDDAAAAAAAAAAAAAAAAASSSSSSSSSSSS:",ORF_1)
	ORF_1.sort(key = lambda codon_pos_set: codon_pos_set[1])	#these three lines sort my ORF dictionaries
	ORF_2.sort(key = lambda codon_pos_set: codon_pos_set[1])
	ORF_3.sort(key = lambda codon_pos_set: codon_pos_set[1])

	ORF_dict={"ORF_1":ORF_1,"ORF_2":ORF_2,"ORF_3":ORF_3} #this makes it easier to loop through all of my open reading frames
#everything below, I'm trying to basically figure out how to make my ORFs...
	return ORF_dict #dictionary where the key is the open reading frame, and the value is a set of all codons and their positions (sorted by position)
	#ORF_dict looks like: {ORF_1:[[codon, codon position]..], ORF_2:[[codon, codon position]..], ...} where the codons are ordered in increasing position


#AFTER WE INDEX AND RECORD ALL POSITIONS, WE NEED TO BE ABLE TO USE THOSE INDEXES TO BUILD THE ACTUAL ORF SEQUENCES!!!

def ORF_sequences(fasta_sequence, orf_codon_pos_dict): #ORF indexes is a dictionary that looks like: {ORF_1:[[codon, codon position]..]} (ordered by pos)
	orf_seq_dict = {}
	open_reading_frame_positions = []
	for reading_frame in orf_codon_pos_dict:
		list_index1=0
		orf_codon_list = orf_codon_pos_dict[reading_frame]
		while list_index1 < len(orf_codon_list):
			codon1 = orf_codon_list[list_index1][0]
			codon1_pos = orf_codon_list[list_index1][1]
			list_index2 = list_index1 + 1
			while list_index2 < len(orf_codon_list):
				codon2 = orf_codon_list[list_index2][0]
				codon2_pos = orf_codon_list[list_index2][1]
				if codon1 != "ATG":
					break
				elif codon2 in ["TAG","TGA","TAA"]:
					orf_sequence = fasta_sequence[codon1_pos-1:codon2_pos+2]
					orf_length = codon2_pos+2 - codon1_pos-1
					open_reading_frame_positions.append([orf_sequence, orf_length, reading_frame, codon1_pos, codon2_pos])
					break
				list_index2 += 1
			list_index1 += 1
	open_reading_frame_positions.sort(key = lambda open_reading_frame:open_reading_frame[1]) #sorts by length of the orf
	#orf_seq_dict[reading_frame]=open_reading_frame_positions
	print(open_reading_frame_positions)
	return open_reading_frame_positions #this is sorted for every reading frame

def biggest_orf_per_reading_frame(orf_seq_list): #this just returns our largest open reading frames for all three forward reading frames
	biggest_orfs_per_rf = []
	reading_frames = ["ORF_1", "ORF_2", "ORF_3"]
	#print("my orf_seq_dict is WOAH", orf_seq_dict)
	orf_seq_list_index = len(orf_seq_list)-1
	while len(reading_frames)>0 and orf_seq_list_index >= 0:
		orf_reading_frame = orf_seq_list[orf_seq_list_index][2] #reading frame for an orf in the orf_seq_list

		if orf_reading_frame in reading_frames: #assumes list is sorted by length
			reading_frames.remove(orf_seq_list[orf_seq_list_index][2])
			print("my reading frames are:", reading_frames)
			biggest_orfs_per_rf.append(orf_seq_list[orf_seq_list_index]) #adds to the final big list
		orf_seq_list_index -= 1

	return biggest_orfs_per_rf

def ORF_reader(fasta_sequence): #this will combine all the above functions into one little function to produce:
						#our list of open reading frames
						#our largest open reading frames per reading frame
	index = codon_indexer(fasta_sequence)
	orf_seqs = ORF_sequences(fasta_sequence, index)
	biggest_orfs = biggest_orf_per_reading_frame(orf_seqs)
	#print("This is a dictionary of our orfs per reading frame:", orf_seqs)
	print("\nThis is a dictionary of our biggest orfs per reading frame:", biggest_orfs)
	return biggest_orfs



def kmer_finder(fasta_dict,k=2):
	kmers=[]
	#first find all possible kmers in entire fasta file
	for identifier in fasta_dict:
		sequence = fasta_dict[identifier]
#		sequence_list=list(seq) #converst sequence to a list
		for i in range(len(sequence)):
			kmer=sequence[i:i+k]
			if kmer not in kmers and len(kmer)==k:
				kmers.append(kmer)
	return kmers

def kmer_repeats_searcher2(fasta_dict,kmers):
	#find all repeats kmer by kmer
	repeats={}
	for kmer in kmers:
		repeats[kmer] = 0

	for identifier, sequence in fasta_dict.items(): #simplifies the whole ordeal using tuple unpacking!

		for kmer in kmers:
			kmer_index = sequence.find(kmer) #gives me the index of the very first kmer in the fasta file

			while kmer_index > -1 and kmer_index <= len(sequence) - len(kmer):
				repeats[kmer] = repeats[kmer]+1
				kmer_index = sequence.find(kmer, kmer_index+1)
#	print(repeats)
	sorted_repeats = sorted(repeats.items(), key=lambda item: item[1])
	print(sorted_repeats)

parser = argparse.ArgumentParser()
parser.add_argument("--echo", help="echo the string u put here")
parser.add_argument("--Long", help="gives longest sequences, their lengths, and their identifiers", action="store_true")
parser.add_argument("filename", help="gimme da file name plz")
args = parser.parse_args()
fasta_input_dict=FASTA_dict(args.filename)


fasta_seq_test = fasta_input_dict["gi|142022655|gb|EQ086233.1|16"]

ORF_reader(fasta_seq_test)



"""
codon_index_dict = codon_indexer(fasta_seq_test)
print("This is our codon positions dictionary:" , codon_index_dict)
print("\nThis is our open reading frames positions:\n")
orf_seq = ORF_sequences(fasta_seq_test,codon_index_dict)
biggest_orf_per_reading_frame(orf_seq)
#print(fasta_input_dict)

"""

print("the following is a sorted list of all of our fasta identifiers: \n", FASTA_dict_lengths(fasta_input_dict))


##WHAT I NEED TO DO HERE IS TO REORDER EVERYTHING BY MAKING A NEW LIST OR DICTIONARY TO FLATTEN OUT OUR DICTIONARY AND THEN SORT
#IT WOULD PROBABLY BE EASIER TO FLATTEN IT ALL INTO A LIST OF LISTS FOR EACH LONGEST ORF TO HAVE A SORT OF "OBJECT"
longest_orfs_ordered=[]
for identifier in fasta_input_dict:
	sequence=fasta_input_dict[identifier]
	orf_reader = ORF_reader(sequence)
	for orf in orf_reader:
		longest_orfs_ordered.append([identifier,orf[0],orf[1],orf[2],orf[3],orf[4]])
#need to find a way to sort this list of list pairs with a list in the second index...
longest_orfs_ordered.sort(key=lambda item: item[2])

print("\nPINGAS,",longest_orfs_ordered,"\nPINGAS\n")
#longest_orfs_ordered.sort(key=lambda orf_object: orf_object[1][3])

#print("\nThe following is the longest ORF in the fasta file:\n" , longest_orfs_ordered[-1])

#ORF_reader(fasta_input_dict["gi|142022655|gb|EQ086233.1|16"])


k=7
kmers = kmer_finder(fasta_input_dict,k)
kmer_repeats_searcher2(fasta_input_dict, kmers)

