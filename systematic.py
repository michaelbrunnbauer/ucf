#!/usr/bin/python

import sys,copy,random
from family import familymember,family

# search all separating families for a given universe for families where
# the sum of (frequency - n/2 + 1) for all abundant elements is lower than m

# supply size of universe as command line argument
m=int(sys.argv[1])

# search for separating families better than bestscore by recursively removing
# basis sets from A. dontremove is a list of members that should not be removed
def search(A,dontremove):
    global bestscore
    # make a list of members to remove so that the family stays separating
    # and has the same universe
    candidates=[]
    for member in A.basis_sets:
        # removing the empty set is pointless
        if not member:
            continue
        if member in dontremove:
            continue
        B=copy.deepcopy(A)
        B.remove(member)
        if not B.is_separating():
            continue
        if len(B.elem_count) < m:
            continue
        # B is already union closed but we have to recalculate the basis sets
        l=len(B)
        B.check_union_closed_conjecture()
        assert len(B)==l
        score=B.abundant_elements_total()
        # we have a new optimum, print it
        if score < bestscore:
            B.stats("best")
            bestscore=score
        candidates.append((score,random.random(),member,B))
    # try promising candidates first, then randomly
    candidates.sort()
    # this list will now be extended by the paths we have followed so that
    # paths do not converge later via different removal order of the same
    # members
    dontremove=list(dontremove)
    # list of families that have been searched
    done=[]
    for score,dummy,member,A in candidates:
        # if two possibilities are isomorphic, search only one
        isomorphic=False
        for B in done:
            if A.is_isomorphic_to(B):
                isomorphic=True
                break
        if isomorphic:
            # every possibility without member is isomorphic to a
            # possibility we have already explored so removing member is not
            # necessary from now on
            dontremove.append(member)
            continue
        search(A,dontremove)
        # we have now searched every possibility without member so removing
        # member is not necessary from now on
        dontremove.append(member)
        done.append(A)

# build the power set of the universe and start the search with it
A=family()
for elem in range(1,m+1):
    A.add(familymember([elem]))
A.unionclose()
bestscore=A.abundant_elements_total()
search(A,[])
