#!/usr/bin/python

# will randomly generate union closed and separating families and save them in
# a MySQL database table with the structure defined in create_table.txt
#
# -the field family contains a serialization of the basis sets generated with 
#  repr()
# -abundant_elements_total is the sum of frequency - n/2 + 1 for all abundant
#  elements
# -target_min_member_card and target_max_member_card are the targets for
#  random generation. min_member_card/max_member_card will be in this range.
# -if basis_sets_optimized is null, the family was generated randomly. if not,
#  the family was generated from random families having the specified
#  universe and number of basis sets with a genetic algorithm and
#  basis_sets_optimized is the number of basis sets of the optimized family

import sys,random,MySQLdb
from math import factorial
from family import familymember
from reproduction import fertilefamily,train,popstats,trim1

# generation range for m
min_universe_card=13
max_universe_card=13

# generation range for number of basis sets
min_basis_sets=3
max_basis_sets=14 # union close currently needs up to 2**n operations

# generation range for minimum member cardinality
min_member_card_min=1
min_member_card_max=13

# generation range for maximum member cardinality (basis sets)
max_member_card_min=1
max_member_card_max=13

# MySQL database credentials
dbname=''
dbhost=''
dbuser=''
dbpasswd=''

# database query function
db=None
def query(query,params=None):
    global db
    if (db==None):
        db=MySQLdb.connect(db=dbname,host=dbhost,user=dbuser,passwd=dbpasswd)
    c=db.cursor()
    if params==None:
        c.execute(query)
    else:
        c.execute(query,params)
    db.commit()
    return c

# random seed can optionally be supplied as command line argument
if len(sys.argv)==2:
    seed=long(sys.argv[1])
else:
    seed=random.randrange(2**32)
random.seed(seed)
print "random seed",seed

# calculate number of sets w. regard to size of universe and size of members
search_space_stats={}
for m in range(min_universe_card,max_universe_card+1):
    stats=[]
    for card in range(1,m+1):
        stats.append(factorial(m)/(factorial(m-card)*factorial(card)))
    search_space_stats[m]=stats

# get weighted random choice of member size for universe size m,
# min_member_card, and max_member_card
def random_member_size(m,min_member_card,max_member_card):
    assert max_member_card <= m
    assert min_member_card <= max_member_card
    stats=search_space_stats[m]
    stats=stats[min_member_card-1:max_member_card]
    total=sum(stats)
    setindex=random.randrange(1,total+1)
    for i,count in enumerate(stats):
        setindex-=count
        if setindex <= 0:
            return i+min_member_card
    assert False

def main():
    for m in range(min_universe_card,max_universe_card+1):
        for basis_sets in range(min_basis_sets,max_basis_sets+1):
            for target_min_member_card in range(min_member_card_min,min_member_card_max+1):
                for target_max_member_card in range(max_member_card_min,max_member_card_max+1):
                    genfamilies(m,basis_sets,target_min_member_card,target_max_member_card)

# try to generate up to 1000 random families with the given parameters and
# save them in the database. abort if 1000 consecutive generation attempts 
# fail. run genetic algorithm on the results and save the best optimized
# families for every resulting number of basis sets
def genfamilies(m,basis_sets,target_min_member_card,target_max_member_card):
    # abort if parameters obviously do not make sense
    if target_min_member_card > target_max_member_card:
        return
    if target_max_member_card*basis_sets < m:
        return

    # if the program is run several times, it will retry parameters tuples
    # for which no random family has been generated
    for cnt, in query("select count(*) from families where basis_sets_optimized is null and m=%s and basis_sets=%s and target_min_member_card=%s and target_max_member_card=%s",(m,basis_sets,target_min_member_card,target_max_member_card)):
        if cnt:
            return

    print "m",m,"basis sets",basis_sets,"min_member_card",target_min_member_card,"max_member_card",target_max_member_card
    sys.stdout.flush()

    universe=range(1,m+1)

    # our population as list and set
    pop=[]
    popset=set()

    # keep track of what we have tried
    initial_families=set()
    consecutivefailures=0

    while consecutivefailures < 1000:
        consecutivefailures+=1
        A=[]
        # pick members (unbiased)
        while len(A) < basis_sets:
            # pick weighted random member size
            card=random_member_size(m,target_min_member_card,target_max_member_card)
            # choose elements randomly from universe
            # try it up to basis_sets*10 times in case search space is small
            for dummy in xrange(basis_sets*10):
                new_member=familymember(random.sample(universe,card))
                if new_member not in A:
                    break
            # search space too small
            if new_member in A:
                break
            A.append(new_member)

        A=frozenset(A)
        if A in initial_families:
            continue
        initial_families.add(A)

        A=fertilefamily(A)
        if A.fitness is None:
            continue
        if A.m < m:
            continue
        if len(A)!=basis_sets:
            continue
        if A in popset:
            continue

        consecutivefailures=0

        pop.append(A)
        popset.add(A)

        save(A,basis_sets,target_min_member_card,target_max_member_card,False)

        if not len(pop) % 10:
            sys.stdout.write('-')
            sys.stdout.flush()

        if len(pop)==1000:
            break

    # initial population generated and saved. run genetic algorithm for every
    # number of basis sets
    pops={basis_sets:(pop,popset)}
    while pops:
        n=max(pops.keys())
        pop,popset=pops[n]
        pop.sort()
        train(pops,100,target_min_member_card,target_max_member_card,False,
              convergence=100)
        # save result
        if pop:
            save(pop[0],basis_sets,target_min_member_card,target_max_member_card,True)
        trim1(pops)
        del pops[n]

# save family in database
def save(A,basis_sets,target_min_member_card,target_max_member_card,optimized):
    min_cardinality=A.min_cardinality()
    max_cardinality=A.max_cardinality()
    avg_cardinality=A.avg_cardinality()
    members=repr(A.members)
    if optimized:
        basis_sets_optimized=len(A)
    else:
        assert len(A) == basis_sets
        basis_sets_optimized=None
    query("insert ignore into families (family,m,n,abundant_elements_total,abundant_elements,basis_sets,totalsize,target_min_member_card,target_max_member_card,min_member_card,max_member_card,avg_member_card,basis_sets_optimized) values (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",(members,A.m,A.n,A.abundant_elements_total,A.abundant_elements,basis_sets,A.totalsize,target_min_member_card,target_max_member_card,min_cardinality,max_cardinality,avg_cardinality,basis_sets_optimized))

main()
