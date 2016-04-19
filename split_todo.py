#!/usr/bin/python

import sys
import cPickle as pickle

# first argument is the saved todo file
todo_file=sys.argv[1]
# second argument is the number of parts to generate
parts=int(sys.argv[2])

assert todo_file.endswith('.pickle')

f=open(todo_file,"rb")
m,todo=pickle.load(f)
f.close()
todo_file=todo_file[:-7]
perpart=len(todo) // parts
if len(todo) % parts:
    perpart+=1

index=0

def save():
    global index
    f=open(todo_file+str(index)+'.pickle',"wb")
    pickle.dump((m,todo1),f,-1)
    f.close()
    index+=1

todo1=[]
for elem in todo:
    todo1.append(elem)
    if len(todo1)==perpart:
        save()
        todo1=[]
if len(todo1):
    save()
