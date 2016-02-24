import sys,copy

# in case we want to replace frozenset with something else later
familymember=frozenset

# returns the indexes of the 1-bits of n. used for the worst case union
# close algorithm, which does the union of all subsets of the family.
bitindex_cache={}
def bitindexes(n):
    indexes=bitindex_cache.get(n)
    if indexes is not None:
        return indexes
    indexes=[]
    n1=n
    i=0
    while n1:
        if n1 & 1:
            indexes.append(i)
        n1>>=1
        i+=1
    bitindex_cache[n]=tuple(indexes)
    return indexes

# returns the powerset of a set of cardinality n as list of lists of member
# indexes - excluding the empty set. used for the worst case union close 
# algorithm as described above.
pset_cache={}
def pset(n):
    indexes=pset_cache.get(n)
    if indexes is not None:
        return indexes
    indexes=[]
    for subsetindex in xrange(1,2**n):
        indexes.append(bitindexes(subsetindex))
    pset_cache[n]=tuple(indexes)
    return indexes

class family(object):

    def __init__(self):
        # the members as list. accessing the members by index can be useful
        self.as_list=[]
        # the members as dictionary
        self.as_dict={}
        # optional element frequencies
        self.elem_count=None
        # basis sets will be computed when doing the union close
        self.basis_sets=None

    # make sure that element frequencies are in self.elem_count
    def require_elem_count(self):
        if self.elem_count is not None:
            return
        self.elem_count={}
        for member in self.as_list:
            for elem in member:
                self.elem_count[elem]=self.elem_count.get(elem,0)+1

    # is the average frequency >= 1/2 ?
    def avg_size_too_big(self):
        self.require_elem_count()
        target=len(self)*len(self.elem_count)
        target=target/2 + target%2
        return sum(self.elem_count.values()) >= target

    # will output a message, statistics and the optional basis sets
    def stats(self,message):
        self.require_elem_count()
        flags = ''
        if len(self.elem_count) == len(self) -1:
            flags+=' normalized'
        flags+=' '+str(self.basis_sets)
        print message+': n='+str(len(self))+' m='+str(len(self.elem_count))+', '+str(len(self.abundant_elements()))+'/'+str(self.abundant_elements_total())+' of '+str(self.total_size())+' abundant'+flags
        sys.stdout.flush()

    # size of the family
    def __len__(self):
        return len(self.as_list)

    # add a member
    def add(self,member):
        i=self.as_dict.get(member)
        if i is not None:
            return i
        assert type(member) is familymember
        if self.elem_count is not None:
            self.elem_count=None
        i=len(self.as_list)
        self.as_list.append(member)
        self.as_dict[member]=i
        return i

    # remove a member (not fast)
    def remove(self,member):
        i=self.as_dict[member]
        if self.elem_count is not None:
            self.elem_count=None
        del self.as_list[i]
        del self.as_dict[member]
        while i < len(self.as_list):
            member=self.as_list[i]
            self.as_dict[member]=i
            i+=1

    # is the family separating?
    def is_separating(self):
        self.require_elem_count()
        l=list(self.elem_count.keys())
        for i,elem1 in enumerate(l):
            for elem2 in l[i+1:]:
                ok=False
                for member in self.as_list:
                    if elem1 in member and elem2 not in member:
                        ok=True
                        break
                    if elem2 in member and elem1 not in member:
                        ok=True
                        break
                if not ok:
                    return False
        return True

    # total size (sum of frequency of all elements)
    def total_size(self):
        self.require_elem_count()
        return sum(self.elem_count.values())

    # total abundance (sum of frequency - n/2 + 1 of all abundant elements)
    def abundant_elements_total(self):
        self.require_elem_count()
        min_len=len(self)
        min_len = min_len / 2.0
        sum=0
        for count in self.elem_count.values():
            if count >= min_len:
                sum+=count-min_len+1
        return sum

    # returns abundant elements
    def abundant_elements(self):
        self.require_elem_count()
        min_len=len(self)
        min_len = min_len / 2 + min_len % 2
        result=[]
        for elem,count in self.elem_count.items():
            if count >= min_len:
                result.append(elem)
        return result

    # do the union close and check the union closed conjecture
    def check_union_closed_conjecture(self):
        self.unionclose()
        self.require_elem_count()
        min_len=len(self)
        min_len = min_len / 2 + min_len % 2
        for count in self.elem_count.values():
            if count >= min_len:
                return
        print
        print self.as_list
        sys.stdout.flush()
        print
        self.stats("counterexample found")
        sys.exit(0)

    # worst case union close (used only if necessary)
    # do the union of the members for all subsets and calculate the basis sets
    def unionclose1(self):
        self.basis_sets=set(self.as_list)
        n=len(self)
        for indexes in pset(n):
            # combine members defined by indexes with union and add the result
            union=None
            for setindex in indexes:
                if union is None:
                    union=self.as_list[setindex]
                else:
                    union=union.union(self.as_list[setindex])
            assert union is not None
            i=self.add(union)
            # i is the index of union in the family
            # i < n means union was in the initial family
            # i in indexes means union is one of the members that were combined
            if i < n and i not in indexes:
                self.basis_sets.discard(union)

    # optimistic union close
    # tries recursive approach on a copy first. switches to unionclose1() as 
    # soon as it will be faster
    # calculates the basis sets
    # will add the empty set if missing
    def unionclose(self,A=None):
        # first call. make a copy of ourself in A and work with the copy
        if A is None:
            self.basis_sets=set(self.as_list)
            A=copy.deepcopy(self)

        # if we have more pairs than in the powerset of our inital version,
        # switch to worst case approach
        n=len(A)
        if (n*(n-1)/2) >= 2**len(self):
            self.unionclose1()
            self.add(familymember([]))
            return

        # do the union for all pairs and recurse if we got new members
        for i1 in range(n):
            for i2 in range(i1+1,n):
                first_set=A.as_list[i1]
                second_set=A.as_list[i2]
                union=first_set.union(second_set)
                if union==first_set or union==second_set:
                    continue
                A.add(union)
                self.basis_sets.discard(union)

        if len(A)==n:
            # no new members, we are finished and copy the results from A
            self.as_list=A.as_list
            self.as_dict=A.as_dict
            self.elem_count=A.elem_count
            self.add(familymember([]))
        else:
            # new members -> recursion
            self.unionclose(A)
