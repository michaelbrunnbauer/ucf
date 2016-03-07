#!/usr/bin/python

import sys,copy,random
from family import familymember,family

# search all separating families for a given universe for families where
# the sum of (frequency - n/2 + 1) for all abundant elements is lower than m

# supply size of universe as command line argument
m=int(sys.argv[1])

# parent_tracking is much faster but needs more mem
parent_tracking=True

# do not search branches isomorphic to ones already searched?
# if you want a reproducible count of separating families, set to False.
# otherwise branches will be searched or not depending on the order of the
# family members. This order determines dontremove, which then overlaps more 
# or less with isomorphic candidates
avoid_isomorphisms=True

# recursion limit may need to be increased for higher m
#rlimit=sys.getrecursionlimit()
#sys.setrecursionlimit(2*rlimit)

# sorts members (try bigger members first, then randomly)
def member_cmp(A,B):
    if len(A) > len(B):
        return -1
    elif len(A) < len(B):
        return 1
    return random.choice([-1,1])

# search for separating families better than bestscore by recursively removing
# basis sets from A. dontremove is a list of members that should not be removed
def search(A,dontremove):
    global bestscore,cnt

    #if len(dontremove):
    #    sys.exit(0)

    # this list will now be extended by the paths we have followed so that
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

        if not parent_tracking:
            B.check_union_closed_conjecture()

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

# build the power set of the universe and start the search with it
A=family(parent_tracking=parent_tracking)
for elem in range(1,m+1):
    A.add(familymember([elem]))
A.unionclose()
bestscore=A.abundant_elements_total()
cnt=0
search(A,[])
print cnt,"separating families searched"
