import sys,random
from family import familymember,emptymember,family

# elementwise addition of tuples A and B
def tuple_add(A,B):
    A=list(A)
    for i,v in enumerate(B):
        A[i]+=v
    return tuple(A)

# elementwise division of tuples A by b
def tuple_div(A,b):
    A=list(A)
    for i,v in enumerate(A):
        A[i]/=b
    return tuple(A)

# add a random element to a random member of family A while respecting
# maximum member cardinality. will return the new family as list
def add_elem(A,max_member_card):
    A=list(A)
    member=random.choice(A)
    elem=random.choice(list(member))
    random.shuffle(A)
    for member in A:
        if elem not in member and len(member) < max_member_card:
            A.remove(member)
            member=list(member)
            member.append(elem)
            A.append(familymember(member))
            return A
    return None

# remove a random member from A
def remove_member(A):
    A=list(A)
    member=random.choice(A)
    A.remove(member)
    return A

# remove a random element from a random member of A while respecting minimum
# member cardinality
def remove_elem(A,min_member_card):
    A=list(A)
    random.shuffle(A)
    for member in A:
        if len(member) > min_member_card:
            break
    if len(member) <= min_member_card:
        return None
    A.remove(member)
    member=list(member)
    member.remove(random.choice(member))
    A.append(familymember(member))
    return A

# swap two random distinct elements from two random members where each
# element is not already present in the other member
def swap_elem(A):
    A=list(A)
    srcmember=random.choice(A)
    while True:
        dstmember=random.choice(A)
        if dstmember!=srcmember:
            break
    A.remove(srcmember)
    A.remove(dstmember)
    srcmember=list(srcmember)
    dstmember=list(dstmember)

    random.shuffle(srcmember)
    for srcelem in srcmember:
        if srcelem not in dstmember:
            break
    if srcelem in dstmember:
        return None

    random.shuffle(dstmember)
    for dstelem in dstmember:
        if dstelem not in srcmember:
            break
    if dstelem in srcmember:
        return None

    srcmember.remove(srcelem)
    dstmember.remove(dstelem)
    srcmember.append(dstelem)
    dstmember.append(srcelem)
    A.append(familymember(srcmember))
    A.append(familymember(dstmember))
    return A

# move a random element from a random member to another random member where
# this element is not yet present while respecting minimum and maximum
# member cardinality
def move_elem(A,min_member_card,max_member_card):
    A=list(A)
    random.shuffle(A)
    for dstmember in A:
        if len(dstmember) < max_member_card:
            break
    if len(dstmember) >= max_member_card:
        return None
    for srcmember in A:
        if dstmember!=srcmember and len(srcmember) > min_member_card:
            break
    if dstmember==srcmember or len(srcmember) <= min_member_card:
        return None
    A.remove(srcmember)
    A.remove(dstmember)
    srcmember=list(srcmember)
    dstmember=list(dstmember)
    random.shuffle(srcmember)
    for srcelem in srcmember:
        if srcelem not in dstmember:
            break
    if srcelem in dstmember:
        return None
    srcmember.remove(srcelem)
    dstmember.append(srcelem)
    A.append(familymember(srcmember))
    A.append(familymember(dstmember))
    return A

# a union closed family that can reproduce using the functions defined above
# BEWARE: the empty set will not be counted as basis set but will be implicitly
# added to the union closed variant
class fertilefamily(object):

    # initialize from a list of members
    def __init__(self,B):
        self.fitness=None

        A=family()
        for member in B:
            A.add(member)

        # we want only separating families
        if not A.is_separating():
            return

        m=len(A.elem_count)

        # do the union close, add the empty set and get the basis sets 
        # (without empty set)
        A.check_union_closed_conjecture()

        B=frozenset([member for member in A.basis_sets if member!=emptymember])

        # see if this is a case for which the union closed conjecture is proven
        n=len(A)
        self.conjecture_proven=False
        if n <= 50:
            self.conjecture_proven=True
        elif n >= ((2**m)*2)/3:
            self.conjecture_proven=True
        elif n <= 2*m:
            self.conjecture_proven=True

        # reproduction attempts
        self.age=0
        # number of viable childs
        self.childs=0

        self.n=n
        self.m=m
        self.members=B
        self.abundant_elements=len(A.abundant_elements())
        self.abundant_elements_total=A.abundant_elements_total()
        self.totalsize=A.total_size()
        # fitness is a tuple so we can sort by several parameters
        # currently only the total abundance is used
        self.fitness=(self.abundant_elements_total,0)
        # accumulated fitness of childs
        self.childfitness=(0,0)

    # so we can have the structure in a dictionary or set
    def __hash__(self):
        return self.members.__hash__()

    # size of the family
    def __len__(self):
        return len(self.members)

    # member in self ?
    def __contains__(self, key):
        return key in self.members

    # for member in self
    def __iter__(self):
        return self.members.__iter__()

    # compare by fitness as first parameter and then randomly
    def __cmp__(self,other):
        # lower values are better!

        if self.fitness < other.fitness:
            return -1
        if self.fitness > other.fitness:
            return 1

        # does not help. avgchildfitness also does not help
        #fertility=(self.age+1)/float(self.childs+1)
        #fertility1=(other.age+1)/float(other.childs+1)
        #if fertility < fertility1:
        #    return -1
        #if fertility > fertility1:
        #    return 1

        return random.choice([-1,1])

    # self == other
    def __eq__(self,other):
        return self.members==other.members

    # self != other
    def __ne__(self,other):
        return self.members!=other.members

    # string representation: some statistics and members
    def __repr__(self):
        if self.childs:
            avgchildfitness=tuple_div(self.childfitness,float(self.childs))
        else:
            avgchildfitness=None
        fertility=(self.age+1)/float(self.childs+1)
        return str(self.abundant_elements_total)+' of '+str(self.totalsize)+' total elements abundant, '+str(self.abundant_elements)+' elements abundant, card '+str(len(self))+', n '+str(self.n)+', age '+str(self.age)+', fertility '+str(fertility)+', avgchildfitness '+str(avgchildfitness)+', members '+repr(self.members)

    # minimum member cardinality (without empty set)
    def min_cardinality(self):
        min=None
        for member in self.members:
            l=len(member)
            # do not count empty set
            if not l:
                continue
            if min is None or min > l:
                min=l
        return min

    # maximum member cardinality
    def max_cardinality(self):
        max=None
        for member in self.members:
            l=len(member)
            if max is None or max < l:
                max=l
        return max

    # average member cardinality (without empty set)
    def avg_cardinality(self):
        sum=0.0
        for member in self.members:
            sum+=len(member)
        return sum/len(self.members)

    # attempt reproduction with same universe and constraints on
    # minimum and maximum member cardinality. if ignore_if_conjecture_proven
    # is True and we get a case where the conjecture is proven in general,
    # reproduction fails
    def reproduce(self,min_member_card,max_member_card,
                  ignore_if_conjecture_proven):
        self.age+=1
        A=self

        # occasionally generate offspring with lower cardinality on purpose
        # does not help
        #if random.choice([0,1])==1:
        #    B=remove_member(A)
        #    B=fertilefamily(B)
        #    if B.fitness is not None and B.m==self.m:
        #        return B

        operations=[]
        operations.append('Removeelem')
        operations.append('Swap')
        operations.append('Move')
        operations.append('Addelem') # seems to be important for good results

        # try up to 3 mutations
        mutations=random.randrange(1,4)
        while mutations:
            operation=random.choice(operations)
            if operation=='Removeelem':
                B=remove_elem(A,min_member_card)
            elif operation=='Swap':
                B=swap_elem(A)
            elif operation=='Move':
                B=move_elem(A,min_member_card,max_member_card)
            elif operation=='Addelem':
                B=add_elem(A,max_member_card)
            else:
                assert False,operation
            if B is not None:
                A=B
            mutations-=1

        if A is None or type(A) is fertilefamily:
            return None

        A=fertilefamily(A)
        if A.fitness is None or A.m!=self.m:
            return None

        if ignore_if_conjecture_proven and A.conjecture_proven:
            return None

        self.childs+=1
        self.childfitness=tuple_add(self.childfitness,A.fitness)
        return A

# population as list and set. sort by fitness and trim to popsize
def trim(pop,popset,popsize):
    pop.sort()
    for A in pop[popsize:]:
        popset.remove(A)
    del pop[popsize:]

# population as several tuples of (list,set) - indexed by the number of
# basis sets in a dictionary. this function will remove all families with
# fitness worse than the best fitness of the highest basis set number index
def trim1(pops):
    n=max(pops.keys())
    pop,popset=pops[n]
    if not pop:
        return
    pop.sort()
    best=pop[0].abundant_elements_total
    for pop,popset in pops.values():
        for A in list(pop):
            if A.abundant_elements_total > best:
                pop.remove(A)
                popset.remove(A)

# population as several tuples of (list,set) - indexed by the number of
# basis sets in a dictionary. train the highest basis set number index until
# convergence. every iteration will attempt a reproduction of the whole
# population and then trim the population to popsize. constraints on minimum
# and maximum member cardinality are respected. if ignore_if_conjecture_proven
# is True, childs with properties for which the union closed conjecture is
# proven in general will be rejected. if stats=True, population stats will
# be printed on every iteration. convergence is the maximum number of 
# iterations for which a change in minimum or maximum fitness is expected.
# convergence=None will loop forever.
def train(pops,popsize,min_member_card,max_member_card,
          ignore_if_conjecture_proven,stats=False,convergence=None):

    # train only the highest basis set number
    n=max(pops.keys())
    pop,popset=pops[n]
    if not pop:
        return

    # trim to popsize
    trim(pop,popset,popsize)

    best_fitness=pop[0].fitness
    worst_fitness=None
    if len(pop)==popsize:
        worst_fitness=pop[-1].fitness
    convergence_count=0
    if stats:
        reference_set=None
    #seen=set()

    # loop until convergence reached
    while convergence is None or convergence_count < convergence:
        if not stats:
            sys.stdout.write('.')

        assert len(pop) == len(popset)

        # update best_fitness
        if pop[0].fitness < best_fitness:
            best_fitness=pop[0].fitness
            convergence_count=0
            if not stats:
                sys.stdout.write('<')

        # update worst_fitness
        if len(pop)==popsize:
            if worst_fitness is None or pop[-1].fitness < worst_fitness:
                worst_fitness=pop[-1].fitness
                convergence_count=0
                if not stats:
                    sys.stdout.write('>')

        convergence_count+=1

        #seencount=0
        #nseencount=0

        # attempt reproduction for whole population
        for A in list(pop):
            B=A.reproduce(min_member_card,max_member_card,
                          ignore_if_conjecture_proven)
            # reproduction failed?
            if B is None:
                continue

            # check number of basis sets of child
            n=len(B)
            if n < len(A):
                # number of basis sets decreased, save in different population
                if n not in pops:
                    pops[n]=([],set([]))
                pop1,popset1=pops[n]
                if B not in popset1:
                    pop1.append(B)
                    popset1.add(B)
                    if len(pop1) > 2*popsize:
                        trim(pop1,popset1,popsize)                    
            # number of basis sets constant. save if not already in population
            elif B not in popset:
                #if B in seen:
                #    seencount+=1
                #else:
                #    nseencount+=1
                #    seen.add(B)
                pop.append(B)
                popset.add(B)

        #print seencount,"seen",nseencount,"not seen"

        # trim population to popsize
        trim(pop,popset,popsize)

        if stats:
            reference_set=popstats(pop,reference_set)
        sys.stdout.flush()

    if not stats:
        print

# number of different elements in two sets of same length. used to calculate
# a measure of diversity
def diff(A,B):
    sum=0
    for member in A:
        if member not in B:
            sum+=1
    return sum

# print population statistics and return new reference set used to calculate
# diversity
def popstats(pop,reference_set=None):
    if not pop:
        return reference_set

    pop.sort()
    if reference_set is None or pop[0].fitness < reference_set.fitness:
        reference_set=pop[0]
    diffsum=0
    for A in pop:
        if A is reference_set:
            continue
        diffsum+=diff(A,reference_set)
    if diffsum:
        diffsum/=float(len(pop)-1)

    print "best",pop[0].fitness,"worst",pop[-1].fitness,"pop size",len(pop),"diversity",diffsum
    print pop[0]
    print
    sys.stdout.flush()

    return reference_set
