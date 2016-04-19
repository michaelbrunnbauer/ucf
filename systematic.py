#!/usr/bin/python

import sys,copy,random,os
import cPickle as pickle
from family import familymember,emptymember,family

# Search all separating families for a given universe for families where
# the sum of (frequency - n/2 + 1) for all abundant elements is lower than m.
# Only families with the empty set are considered.
# Supply size of universe or todo file as command line argument.
#
# If track_progress=True, you can stop calculation and save the current todos by
# creating a file named "stop". Calculation can be resumed by removing this file
# and calling systematic.py with the todo file as command line argument. The
# todo file can be split for paralell processing with split_todo.py.

# tracking progress is slower and needs more mem but enables stopping, 
# resuming, paralell processing and estimates on runtime. useful for m>5.
track_progress=True

# number of todos (search space branches) to collect for progress meter
min_todo=1000

# parent_tracking is much faster but needs more mem
# mandatory if track_progress=True
parent_tracking=True

if track_progress:
    assert parent_tracking
    assert not os.path.exists('stop')

# do not search branches isomorphic to ones already searched?
avoid_isomorphisms=True

# recursion limit may need to be increased for higher m
#rlimit=sys.getrecursionlimit()
#sys.setrecursionlimit(2*rlimit)

# member compare for sorting (by cardinality first)
def member_cmp(A,B):
    # ascending cardinality causes small members and isomorphisms to be
    # added to dontremove early, which reduces the search space
    if len(A) > len(B):
        return 1
    elif len(A) < len(B):
        return -1
    sA=sum(A)
    sB=sum(B)
    if sA > sB:
        return -1
    elif sA < sB:
        return 1
    if str(A) > str(B):
        return -1
    else:
        return 1

class RuntimeException(Exception):
    pass

# search for separating families better than bestscore by recursively removing
# basis sets from A. dontremove is a list of members that should not be removed
def search(A,dontremove,max_cnt=0):
    global bestscore,cnt

    # dontremove will be extended by the paths we have followed so that
    # paths do not converge later via different removal order of the same
    # members
    dontremove=list(dontremove)

    # list of families that have been searched
    if avoid_isomorphisms:
        done=[]

    # list of members to remove so that the family stays separating
    # and has the same universe
    todo=list(A.basis_sets)
    todo.sort(cmp=member_cmp)
    for member in todo:
        # removing the empty set is pointless
        if not member:
            continue
        if member in dontremove:
            continue
        if not parent_tracking:
            B=copy.deepcopy(A)
        else:
            B=A
        B.remove(member)
        if not B.is_separating():
            if parent_tracking:
                B.unremove()
            continue
        if len(B.elem_count) < m:
            if parent_tracking:
                B.unremove()
            continue

        cnt+=1
        if max_cnt and cnt > max_cnt:
            if parent_tracking:  
                B.unremove()
            raise RuntimeException

        if track_progress and not max_cnt and not cnt % 1000:
            if os.path.exists('stop'):
                if parent_tracking:
                    B.unremove()   
                raise RuntimeException

        if not parent_tracking:
            B.unionclose()

        score=B.abundant_elements_total()
        # we have a new optimum, print it
        if score < bestscore:
            B.stats("best")
            bestscore=score

        if avoid_isomorphisms:
            isomorphic=False
            for C in done:
                if B.is_isomorphic_to(C):
                    isomorphic=True
                    break
            if isomorphic:
                # every possibility without member is isomorphic to a
                # possibility we have already explored so removing member is
                # not necessary from now on
                dontremove.append(member)
                if parent_tracking:
                    B.unremove()
                continue

            if parent_tracking:
                # make a copy without the parent tracking information to check
                # for isomorphisms later
                C=family()
                C.as_list=list(B.as_list)
                C.as_dict=dict(B.as_dict)
                C.basis_sets=set(B.basis_sets)
                C.elem_count=dict(B.elem_count)
            else:
                # no need to make a copy, B is already a deep copy
                C=B
            done.append(C)

        # recurse
        try:
            # no need to recurse if empty set is only candidate to remove
            if len(dontremove) < len(B.basis_sets) - 1:
                search(B,dontremove,max_cnt=max_cnt)
        finally:
            if parent_tracking:
                B.unremove()  

        # we have now searched every possibility without member so removing
        # member is not necessary from now on
        dontremove.append(member)

# this is like search, but it will construct a todo list of branches to
# search instead of recursing
def get_todos(B,dontremove):
    global bestscore,cnt

    assert parent_tracking

    # list of families that have been searched
    if avoid_isomorphisms:
        done=[]

    # list of members to remove so that the family stays separating
    # and has the same universe
    todo=list(B.basis_sets)
    todo.sort(cmp=member_cmp)
    todo_recursion=[]
    for member in todo:
        # removing the empty set is pointless
        if not member:
            continue
        if member in dontremove:
            continue
        B.remove(member)
        if not B.is_separating():
            B.unremove()
            continue
        if len(B.elem_count) < m:
            B.unremove()
            continue

        cnt+=1

        score=B.abundant_elements_total()
        # we have a new optimum, print it
        if score < bestscore:
            B.stats("best")
            bestscore=score

        if avoid_isomorphisms:
            isomorphic=False
            for C in done:
                if B.is_isomorphic_to(C):
                    isomorphic=True
                    break
            if isomorphic:
                # every possibility without member is isomorphic to a
                # possibility we have already explored so removing member is
                # not necessary from now on
                dontremove.append(member)
                B.unremove()
                continue

            # make a copy without the parent tracking information to check
            # for isomorphisms later
            C=family()
            C.as_list=list(B.as_list)
            C.as_dict=dict(B.as_dict)
            C.basis_sets=set(B.basis_sets)
            C.elem_count=dict(B.elem_count)
            done.append(C)

        # no need to recurse if empty set is only candidate to remove
        if len(dontremove) < len(B.basis_sets) - 1:
            # save a reference to B and the member that should be removed
            todo_recursion.append((B,member,list(dontremove),len(B.basis_sets)))

        # restore original state of B
        B.unremove()

        # we have now searched every possibility without member so removing
        # member is not necessary from now on
        dontremove.append(member)

    return todo_recursion

# expand todo list from get_todos by calling get_todos for every todo
def expand_todo(todo):
    assert parent_tracking
    new_todo=[]
    for A,member,dontremove,basis_sets_size in todo:
        A.remove(member)
        # the same A may be referenced several times - so we need a copy
        B=copy.deepcopy(A)
        new_todo+=get_todos(B,dontremove)
        # restore original state of A
        A.unremove()
    return new_todo

# spend some time on removing easy candidates
def remove_easy_todos(todo):
    global cnt
    new_todo=[]
    for B,member,dontremove,basis_sets_size in todo:
        B.remove(member)
        old_cnt=cnt
        try:
            search(B,dontremove,max_cnt=cnt+100)
        except RuntimeException:
            new_todo.append((B,member,dontremove,basis_sets_size))
            cnt=old_cnt
        B.unremove()
    return new_todo

# same as search(), but with progress output
# get min_todo todos (branches of search space) and work through them
def search_progress(todo_recursion):
    assert parent_tracking

    while len(todo_recursion) < min_todo:
        todo_recursion=expand_todo(todo_recursion)
        if not len(todo_recursion):
            break
        if len(todo_recursion) >= min_todo:
            todo_recursion=remove_easy_todos(todo_recursion)

    print len(todo_recursion),"todos"
    sys.stdout.flush()

    # should distribute inaccuracy
    random.shuffle(todo_recursion)

    # estimate search space size from number of basis sets and size of
    # dontremove. a joke, but better than nothing
    def search_space(basis_sets_size,dontremove_size):
        n=basis_sets_size-dontremove_size-1
        assert n>=0
        return n

    search_space_size=0
    for B,member,dontremove,basis_sets_size in todo_recursion:
        search_space_size+=search_space(basis_sets_size,len(dontremove))
    search_space_size=float(search_space_size)

    searched=0
    for B,member,dontremove,basis_sets_size in list(todo_recursion):
        searched+=search_space(basis_sets_size,len(dontremove))
        B.remove(member)
        try:
            search(B,dontremove)
            todo_recursion.remove((B,member,dontremove,basis_sets_size))
        except RuntimeException:
            B.unremove()
            filename='systematic.pickle'
            if sys.argv[1].endswith('.pickle'):
                filename=sys.argv[1]
            print "saving",filename
            f=open(filename,"wb")
            pickle.dump((m,todo_recursion),f,-1)
            return
        B.unremove()
        pct=searched/search_space_size
        print pct*100,"%"
        sys.stdout.flush()

cnt=0
if sys.argv[1].endswith('.pickle'):
    # called with todo file as argument
    assert track_progress
    f=open(sys.argv[1],"rb")
    m,todo_recursion=pickle.load(f)
    bestscore=m
    f.close()
    search_progress(todo_recursion)
else:
    # start new calculation
    m=int(sys.argv[1]) 
    # build the power set of the universe and start the search with it
    A=family(parent_tracking=parent_tracking)
    A.add(emptymember)
    for elem in range(1,m+1):
        A.add(familymember([elem]))
    A.unionclose()
    bestscore=A.abundant_elements_total()
    assert bestscore==m
    cnt+=1
    if track_progress:
        todo_recursion=get_todos(A,[])
        search_progress(todo_recursion)
    else:
        search(A,[])

print cnt,"separating families searched"
