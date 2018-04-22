#/usr/bin/python
#coding:utf-8

import time
import numpy as np
from sys import argv
from math import log
from pandas import DataFrame
from collections import defaultdict

global res_dict, df
res_list = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",\
            "M", "F", "P", "S", "T", "W", "Y", "V"]
res_dict = {res:0 for res in res_list}
df = DataFrame(np.zeros([20, 20]), index=res_list, columns=res_list) #"df" will save the calulate result

#Procese the input data to blocks
def procese_text(blocks):
	for block in blocks:
		block = block.strip().split("\n")[4:]
		while '' in block:
			block.remove('')
		for i in range(len(block)):
			block[i] = block[i].strip().split(' ')
			while '' in block[i]:
				block[i].remove('')
			block[i] = block[i][3] if len(block[i])==5 else block[i][2]
			block[i] = list(block[i])
		yield(block)

#Calulate the identity of pairs sequences
def compare_seqs(seq1, seq2):
	count, leng = 0, len(seq1)
	for i in range(leng):
		if seq1[i] == seq2[i]:
			count += 1
	return(int((count/leng)*100))

#Cluster the Block
def cluster(block):
	group_num = 0
	pass_list = []
	group_list = []
	group = DataFrame()
	for i, row1 in block.iterrows(): #First loop, traverse every row of block
		if i in pass_list: 
			continue
		count = tmp =1
		group_num += 1
		pass_list.append(i)
		group = group.append(row1)
		block = block.drop(i, axis=0)
		while tmp <= count: #While loop; Each group of clustering completes,jump loop
			for j, row2 in group.iterrows(): #Second loop, traverse every row of group
				if j == "#":
					continue
				tmp += 1
				row_name = list(group.index)
				m = row_name.index(j)
				row_name[m] = "#" #The "#" is a flag
				group.index = row_name
				for k, row3 in block.iterrows(): #Thridly loop
					if compare_seqs(row2, row3) >= 62:
						count += 1
						pass_list.append(k)
						group = group.append(row3)
						block = block.drop(k, axis=0)
		group_list += [count]*count
	return(group, group_list, group_num)
#Calulate one columns of block
def calulate(resduces, factor):
	global res_dict, df
	dict_block = {res:0 for res in res_list}
	for i in range(len(resduces)):
		try:
			res_dict[list(resduces)[i]] += 1/factor[i]
			dict_block[list(resduces)[i]] += 1/factor[i]
		except KeyError as e:
			pass
	for x in range(len(res_list)):
		for y in range(len(res_list)):
			i = res_list[x]
			j = res_list[y]
			if x <= y:
				if i  == j:
					df[i][j]+=(dict_block[i]*(dict_block[i]-1))/2
				else:
					df[i][j]+=dict_block[i]*dict_block[j]
			else:
				df[i][j] = df[j][i]

#Prcess the block
def procese_block(block, block_list, block_num):
	K = block_num
	L = len(block.columns)
	T = K * L
	P = int((K*(K-1)*L)/2)
	for name, row in block.iteritems():
		calulate(row, block_list) 
	return(T, P)

#Calulate score
def blosum(T, P):
	global res_dict, df
	df = df/P
	for key in res_dict.keys():
		res_dict[key] /= T
	for x in range(len(res_list)):
		for y in range(len(res_list)):
			i = res_list[x]
			j = res_list[y]
			N = df[i][j]
			M = res_dict[i]*res_dict[j]
			m = round(N/M, 3)
			if x <= y:
				try:
					df[i][j] = round(2*log(m, 2))
				except ValueError as e:
					df[i][j] = 0
			else:
				df[i][j] = df[j][i]
	return(df)
	
def tran_time(times):
	s = times % 60
	tmp = int(times/60)
	m = tmp if tmp<60 else tmp%60
	tmp = int(tmp/60)
	h = tmp if tmp<24 else tmp%24
	d = int(tmp/24)
	return(d, h, m, s) 
	
if __name__ == "__main__":
	Count, Total, Pairs =0, 0, 0
	infile = argv[1]
	data = open(infile, 'r')
	blocks = data.read()
	blocks = blocks.strip().split("//")
	Leng = len(blocks)
	block_gener = procese_text(blocks)
	for block in block_gener:
		Count += 1
		rate = Count*100/Leng
		block = DataFrame(block)
		block, block_list, block_num = cluster(block)
		nums = procese_block(block, block_list, block_num)
		Total, Pairs = Total+nums[0], Pairs+nums[1]
		times = tran_time(time.clock())
		print("<<Done/Total:{0}/{1}>> || <<Running schedule:{2}%>> || <<Running time:{3}d {4}h {5}m {6}s>>".\
		format(Count, Leng, round(rate, 2),times[0], times[1], times[2], round(times[3], 2)), end="\r")
	print("\n")
	out_data = blosum(Total, Pairs)
	out_data.to_csv("out_data_62_mul.csv")
	data.close()
