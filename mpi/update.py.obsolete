#!/usr/bin/python
import os, sys, string

S = open('../write.i', 'r').readlines()
state = 0
var = {'complex':[], 'real':[]}
arr = []

def fill(s):
	T = string.split(s,',')
	for i in range(len(T)):
		T[i] = string.strip(T[i])
	return T

for s in S:
	if s[0] in 'cC!' and s[:4]!='CMPI':
		continue
	pos = string.find(s,'double complex')
	if pos!=-1:
		var['complex'] = var['complex']+fill(s[pos+14:-1])
		continue
	pos = string.find(s,'double precision')
	if pos!=-1:
		var['real'] = var['real']+fill(s[pos+16:-1])
		continue
	pos = string.find(s,'common/write/')
	if pos!=-1:
		state = 1
		arr = arr+fill(s[pos+13:-1])
		continue
	if s[:15]=='CMPI MAINOUTPUT':
		break
	if state==1:
		arr = arr+fill(s[pos+7:-1])
i = 0
while i<(len(arr)):
    if arr[i]=='':
        del arr[i]
        i = i-1
    else:
        arr[i] = arr[i][:string.find(arr[i],'(')]
    i = i+1

src = open('mpi.ins','r')
dst = open('mpi.new','w')
state = 0
s = src.readline()[:-1]

optmpl = {}
optmpl['send'] = '      call MPI_SEND($,nrayelt,#,0,\n     +     tag,MPI_COMM_WORLD,ierr)'
optmpl['recv'] = '      call MPI_RECV($,nrayelt,#,src,\n     +     tag,MPI_COMM_WORLD,status,ierr)'
def prn(op):
    for v in arr:
        s = string.replace(optmpl[op],'$',v)
	type = ''
        if v in var['complex']:
            type = 'MPI_DOUBLE_COMPLEX'
        else:
            type = 'MPI_DOUBLE_PRECISION'
        s = string.replace(s,'#',type)
#        print s
        dst.write(s+'\n')

while s!='':
    if state!=1 and state!=3:
#        print s
        dst.write(s+'\n')
    if s=='CVARCOMM':
        state = state+1
        if state==1:
            prn('send')
#            print 'CVARCOMM'
            dst.write('CVARCOMM\n')
        if state==3:
            prn('recv')
#            print 'CVARCOMM'
            dst.write('CVARCOMM\n')
    s = src.readline()[:-1]

os.system('cp mpi.ins mpi.ins.old')
os.system('mv mpi.new mpi.ins')
