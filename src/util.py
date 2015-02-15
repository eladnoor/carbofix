#!/usr/bin/python2.5

import os
import types

def read_simple_mapfile(filename, default_value=""):
    map = {}
    file = open(filename, 'r')
    for line in file.readlines():
        if (line.find('=') == -1):
            map[line.strip()] = default_value
        else:
            (key, value) = line.split('=')
            map[key.strip()] = value.strip()
    file.close()
    return map

def _mkdir(newdir):
    """works the way a good mkdir should :)
        - already exists, silently complete
        - regular file in the way, raise an exception
        - parent directory(ies) does not exist, make them as well
    """
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError("a file with the same name as the desired " \
                      "dir, '%s', already exists." % newdir)
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        if tail:
            os.mkdir(newdir)

def gcd(a,b=None):
    """ Return greatest common divisor using Euclid's Algorithm.
        a - can be either an integer or a list of integers
        b - if 'a' is an integer, 'b' should also be one. if 'a' is a list, 'b' should be None
    """
    if (b == None):
        if (a == []):
            return 1
        g = a[0]
        for i in range(1, len(a)):
            g = gcd(g, a[i])
        return g
    else:
        while b:      
          a, b = b, a % b
        return a

# choose function (number of way to choose k elements from a list of n)
def choose(n, k):
    """ 
        return the binomial coefficient of n over k
    """
    def rangeprod(k, n):
        """
            returns the product of all the integers in {k,k+1,...,n}
        """ 
        res = 1
        for t in range(k, n+1):
            res *= t
        return res
    
    if (n < k):
        return 0
    else:
        return (rangeprod(n-k+1, n) / rangeprod(1, k))

def subsets(items, minsize=0, maxsize=-1):
    def subsets_recursive(items, size, begin_index=0):
        """
            returns a list of all subsets of "items" of size == "size"
            e.g. subsets([1,2,3], 2) = [[1,2],[1,3],[2,3]]
        """
        if (items == []):
            return []
        elif (size == 0):
            return [[]]
        elif (size == 1):
            return [[x] for x in items[begin_index:]]
        else:
            s = []
            for i in range(begin_index, len(items)-1):
                x = items[i]
                for y in subsets_recursive(items, size-1, begin_index=i+1):
                    s.append([x] + y)
            return s
    
    if (maxsize == -1): # by default, return subsets of all sizes
        maxsize = len(items)
    s = []
    for size in range(minsize, maxsize+1):
        s += subsets_recursive(items, size)
    return s
    
def list2pairs(l):
    """
        Turns any list with N items into a list of (N-1) pairs of consecutive items.
    """
    res = []
    for i in range(len(l)-1):
        res.append((l[i], l[i+1]))
    return res

def sum(l):
    """
        Returns the sum of the items in the container class.
        This is more general than the build-in 'sum' function, because it is not specific for numbers.
        This function uses the '+' operator repeatedly on the items in the contrainer class.
        For example, if each item is a list, this will return the concatenation of all of them
    """
    return reduce(lambda x,y: x+y, l)

def flatten(l):
    res = []
    for x in l:
        if (isinstance(x, types.ListType)):
            res += flatten(x)
        else:
            res.append(x)
    return res

def median(v):
    if (len(v) == 0):
        return None
    sv = sorted(v)
    if (len(v) % 2 == 1):
        return (sv[(len(v)-1)/2])
    else:
        return (sv[len(v)/2 - 1] + sv[len(v)/2]) * 0.5

###############################################################
def test():
    print sum([[1,2,3,4,5],[1,2,3]])
    for s in subsets([(1,'a'),(2,'b'),(3,'c'),(4,'d'),(5,'e')], minsize=1, maxsize=2):
        print s
    
    print flatten([[1,2,3],[4,[5,6],[[7]]]])
      
if (__name__ == '__main__'):
    test()
