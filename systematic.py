#!/usr/bin/python

import sys,copy,random
from family import familymember,emptymember,family

# search all separating families for a given universe for families where
# the sum of (frequency - n/2 + 1) for all abundant elements is lower than m
# only families with the empty set are considered

# supply size of universe as command line argument
m=int(sys.argv[1])

# tracking progress is slower and needs more mem but estimates on runtime are
# crucial
track_progress=False

# number of todos (search space branches) to collect for progress meter
min_todo=1000

# parent_tracking is much faster but needs more mem
# mandatory if track_progress=True
parent_tracking=True

if track_progress:
    assert parent_tracking

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

# search for separating families better than bestscore by recursively removing
# basis sets from A. dontremove is a list of members that should not be removed
def search(A,dontremove):
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
    if avoid_isomorphisms:
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
        search(B,dontremove)

        # we have now searched every possibility without member so removing
        # member is not necessary from now on
        dontremove.append(member)

        if parent_tracking:
            B.unremove()

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
    if avoid_isomorphisms:
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

# same as search(), but with progress output
# get min_todo todos (branches of search space) and work through them
def search_progress(A,dontremove):
    assert parent_tracking

    todo_recursion=get_todos(A,dontremove)
    while len(todo_recursion) and len(todo_recursion) < min_todo:
        todo_recursion=expand_todo(todo_recursion)

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
    for B,member,dontremove,basis_sets_size in todo_recursion:
        searched+=search_space(basis_sets_size,len(dontremove))
        B.remove(member)
        search(B,dontremove)
        B.unremove()
        pct=searched/search_space_size
        print pct*100,"%"
        sys.stdout.flush()

# build the power set of the universe and start the search with it
A=family(parent_tracking=parent_tracking)
A.add(emptymember)
for elem in range(1,m+1):
    A.add(familymember([elem]))
A.unionclose()
bestscore=A.abundant_elements_total()
cnt=1
if track_progress:
    search_progress(A,[])
else:
    search(A,[])

print cnt,"separating families searched"
