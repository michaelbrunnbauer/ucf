import sys,copy,itertools

# in case we want to replace frozenset with something else later
familymember=frozenset

# without empty set!
def powerset_without_0(A):
    s=list(A)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(1,len(s)+1))

class family(object):

    # parent_tracking enables fast recalculation of basis sets when removing
    # basis sets but needs more memory and slows down union close.
    def __init__(self,members=None,parent_tracking=False):
        self.parent_tracking=parent_tracking
        # the members as list. accessing the members by index can be useful
        self.as_list=[]
        # the members as dictionary
        self.as_dict={}
        # optional element frequencies
        self.elem_count=None
        # basis sets and parents will be computed when doing the union close
        self.basis_sets=None
        # self.remove() can be revoked with self.unremove()
        # this is the stack used for it
        self.remove_stack=[]
        # optional parent tracking for members (mother union father = child)
        # which is mother and which is father is determined by self.as_list
        self.parents=None # list of (mother,father) by child
        self.childs_by_mother=None # list of (child,father) by mother
        self.childs_by_father=None # list of (child,mother) by father
        
        if members is not None:
            for member in members:
                self.add(member)

    def __contains__(self, member):
        assert type(member) is familymember
        return member in self.as_dict

    def __iter__(self):
        return self.as_list.__iter__()

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
        target=target//2 + target%2
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

    # remove a member
    # after union close, only basis sets can be removed
    # will recalculate basis sets if self.parent_tracking==True
    def remove(self,member):
        i=self.as_dict[member]
        if self.elem_count is not None:
            self.elem_count=None
        if self.basis_sets is not None:
            # removing non basis set not allowed after union close
            assert member in self.basis_sets
            if self.parent_tracking:
                self.basis_sets.remove(member)
                parent_tracking_actions=self.remove_parent(member)
            else:
                parent_tracking_actions=None
                self.basis_sets=None
        self.remove_stack.append((member,parent_tracking_actions))
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
        min_len = min_len // 2 + min_len % 2
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
        min_len = min_len // 2 + min_len % 2
        for count in self.elem_count.values():
            if count >= min_len:
                return
        print
        print self.as_list
        sys.stdout.flush()
        print
        self.stats("counterexample found")
        sys.exit(0)

    # add parents as parents of child
    def add_parents(self,child,mother,father):
        if child not in self.parents:
            self.parents[child]=set()
        self.parents[child].add((mother,father))
        if mother not in self.childs_by_mother:
            self.childs_by_mother[mother]=set()
        self.childs_by_mother[mother].add((child,father))
        if father not in self.childs_by_father:
            self.childs_by_father[father]=set()
        self.childs_by_father[father].add((child,mother))

    # remove parent, recalculate basis sets, return tuple of actions taken
    def remove_parent(self,parent):
        orphaned=[]
        childs_by_mother_deleted=[]
        childs_by_father_deleted=[]
        parentslist_deleted=[]
        if parent in self.childs_by_mother:
            for child,father in self.childs_by_mother[parent]:
                childs_by_mother_deleted.append((parent,child,father))
                parentslist=self.parents[child]
                assert parentslist
                parentslist.remove((parent,father))
                parentslist_deleted.append((child,parent,father))
                self.childs_by_father[father].remove((child,parent))
                childs_by_father_deleted.append((father,child,parent))
                if not parentslist:
                    orphaned.append(child)

        if parent in self.childs_by_father:
            for child,mother in self.childs_by_father[parent]:
                childs_by_father_deleted.append((parent,child,mother))
                parentslist=self.parents[child]
                assert parentslist
                parentslist.remove((mother,parent))
                parentslist_deleted.append((child,mother,parent))
                self.childs_by_mother[mother].remove((child,parent))
                childs_by_mother_deleted.append((mother,child,parent))
                if not parentslist:
                    orphaned.append(child)

        for orphan in orphaned:
            self.basis_sets.add(orphan)

        return (tuple(orphaned),tuple(childs_by_mother_deleted),tuple(childs_by_father_deleted),tuple(parentslist_deleted))

    # reverse last call of self.remove(). nice for search trees
    def unremove(self):
        member,parent_tracking_actions=self.remove_stack.pop()
        self.add(member)
        if not self.parent_tracking:
            if self.basis_sets is not None:
                self.basis_sets=None
            return
        self.basis_sets.add(member)
        orphaned,childs_by_mother_deleted,childs_by_father_deleted,parentslist_deleted=parent_tracking_actions
        for member in orphaned:
            self.basis_sets.remove(member)
        for child,mother,father in parentslist_deleted:
            self.parents[child].add((mother,father))
        for mother,child,father in childs_by_mother_deleted:
            self.childs_by_mother[mother].add((child,father))
        for father,child,mother in childs_by_father_deleted:
            self.childs_by_father[father].add((child,mother))
        
    # worst case union close (used only if necessary)
    # do the union of the members for all subsets and calculate the basis sets
    # will add the empty set if missing
    def unionclose1(self):
        self.basis_sets=set(self.as_list)
        assert not self.parent_tracking
        n=len(self)
        empty=familymember([])
        for members in powerset_without_0(self.as_list):
            union=empty.union(*members)
            i=self.add(union)
            # i is the index of union in the family
            # i < n means union was in the initial family
            if i < n and union not in members:
                self.basis_sets.discard(union)
        self.add(empty)
        self.basis_sets.add(empty)

    # optimistic union close
    # tries all pairs on a copy first. switches to unionclose1() as 
    # soon as it will be faster
    # calculates the basis sets
    # will add the empty set if missing
    def unionclose(self):

        def pairs(n):
            return (n*(n-1))/2

        self.basis_sets=set(self.as_list)
        if self.parent_tracking:
            self.parents={}
            self.childs_by_mother={}
            self.childs_by_father={}
        A=copy.deepcopy(self)
        worst_case=2**len(self)
        i1=1
        while i1 < len(A):
            # switch to unionclose1 if we have more pairs left than subsets
            # of the original family
            if not self.parent_tracking:
                if pairs(len(A)) - pairs(i1) >= worst_case:
                    self.unionclose1()
                    return
            first_set=A.as_list[i1] 
            for i2 in xrange(i1):
                second_set=A.as_list[i2]
                union=first_set.union(second_set)
                if union==first_set or union==second_set:
                    continue
                A.add(union)
                self.basis_sets.discard(union)
                if self.parent_tracking:
                    self.add_parents(union,first_set,second_set)
            i1+=1
        self.as_list=A.as_list
        self.as_dict=A.as_dict
        self.elem_count=A.elem_count 
        self.add(familymember([]))  
        self.basis_sets.add(familymember([]))

    # returns a dictionary mapping frequencies to elements
    def count_to_elem(self):
        self.require_elem_count()
        result={}
        for elem,count in self.elem_count.items():
            if count not in result:
                result[count]=[]
            result[count].append(elem)
        return result

    # if self isomorphic to other via mapping?
    def is_isomorphic_with_mapping(self,other,mapping):
        def map(member):
            result=[]
            for elem in member:
                result.append(mapping[elem])
            return familymember(result)
        if len(self)!=len(other):
            return False
        if self.basis_sets and other.basis_sets:
            if len(self.basis_sets) != len(other.basis_sets):
                return False
            for member in self.basis_sets:
                if map(member) not in other.basis_sets:
                    return False
        else:
            for member in self.as_list:
                if map(member) not in other:
                    return False
        return True

    # is self isomorphic to other?
    def is_isomorphic_to(self,other):
        assert type(other) is family
        # check size first
        if len(self)!=len(other):
             return False
        if self.basis_sets and other.basis_sets:
            if len(self.basis_sets) != len(other.basis_sets):
                return False
        # compare element frequencies
        self.require_elem_count()
        other.require_elem_count()
        elem_count1=self.elem_count.values()
        elem_count2=other.elem_count.values()
        elem_count1.sort()
        elem_count2.sort()
        if elem_count1 != elem_count2:
            return False

        # map frequencies to elements for other
        count_to_elem=other.count_to_elem()

        # generator for possible mappings
        def maps(elems1,elems2,todo):
            todo=list(todo)
            # map another element and recurse
            elem,count=todo.pop()
            elems1.append(elem)
            # candidates must have the same frequency
            for candidate in count_to_elem[count]:
                if candidate in elems2:
                    continue
                elems2.append(candidate)
                if not todo:
                    assert len(elems1)==len(elems2)
                    yield dict(zip(elems1,elems2))
                else:
                    for mapping in maps(elems1,elems2,todo):
                        yield mapping
                elems2.pop()
            elems1.pop()
        # check possible mappings
        for mapping in maps([],[],self.elem_count.items()):
            if self.is_isomorphic_with_mapping(other,mapping):
                return True
        return False
