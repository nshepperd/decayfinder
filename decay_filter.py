from ROOT import Belle2
from ROOT import TDatabasePDG
from basf2 import inspect, B2INFO, Module
from descriptor import DecayDescriptor, AtomicDescriptor
from descriptor import LogicalDescriptor, ListDescriptor, decode

PDB = TDatabasePDG.Instance()

def println(*args):
    B2INFO(' '.join(map(str, args)))

def print_tree(mcparticles, i=0, level=0):
    p = mcparticles[i]
    pid = p.getPDG()
    mass = p.getMass()
    charge = p.getCharge()
    energy = p.getEnergy()

    entry = PDB.GetParticle(pid)
    if entry:
        name = entry.GetName()
    else:
        name = 'unknown'

    if 0 < level < 7:
        prefix = '\033[3' + str(level) + 'm'
        suffix = '\033[m'
    else:
        prefix = ''
        suffix = ''
    print '    ' * level, prefix + '[%i] %s mass=%f energy=%f charge=%f' % (pid, name, mass, energy, charge) + suffix
    if p.getFirstDaughter() > 0:
        m = p.getFirstDaughter()
        n = p.getLastDaughter()
        for j in range(m-1, n):
            print_tree(mcparticles, j, level + 1)

def recursive_matching(left, right, edges):
    """
    The problem of matching up mcparticle children to the decay descriptor tree
    is actually the 'Bipartite Graph Matching' problem. `left` and `right`
    are disjoint subsets of a graph containing only edges that go between an element
    of `left` and an element of `right`.

    I could use a real algorithm with decent asymptotic properties, but since my
    input is usually of size < 5, something like the Hopcroft-Karp algorithm would
    probably be slower than naive brute force, due to complexity. Not that I've tested
    this.

    So essentially I do a depth first search, for each edge in some set of mutually
    incompatible edges, selecting that edge, and recursing on the same problem with
    that edge, its endpoints, and any other edges connected to those endpoints removed.

    The procedure yields each subset of 'edges' that forms a valid bijection between
    'left' (or some subset of 'left', in the case of inclusive decays) and 'right'.
    That is, all valid ways of matching the monte carlo children to the decay descriptor.
    """
    # Normally, require a bijection.
    if (not left) and (not right):
        # Everything that needs to be is matched.
        yield []
        return
    if not edges:
        # Something is unmatched, and we're out of free edges.
        return

    incompatible = list(edges)
    while incompatible:
        (i, j) = incompatible.pop()
        new_left = left.copy()
        new_left.discard(i)
        new_right = right.copy()
        new_right.discard(j)
        new_edges = [(x, y) for (x, y) in edges if (x != i) and (y != j)]
        for matching in recursive_matching(new_left, new_right, new_edges):
            yield matching + [(i, j)]
        incompatible = [(x, y) for (x, y) in incompatible if (x == i) or (y == j)]

def long_matching(left, right, edges, relatives):
    # Normally, require a bijection.
    if (not left) and (not right):
        # Everything that needs to be is matched.
        yield []
        return
    if not edges:
        # Something is unmatched, and we're out of free edges.
        return

    incompatible = list(edges)
    while incompatible:
        (i, j) = incompatible.pop()
        new_left = left.copy()
        new_left.difference_update(relatives[i])
        new_right = right.copy()
        new_right.remove(j)
        new_edges = [(x, y) for (x, y) in edges if (x not in relatives[i]) and (y != j)]
        for matching in long_matching(new_left, new_right, new_edges, relatives):
            yield matching + [(i, j)]
        incompatible = [(x, y) for (x, y) in incompatible if (x in relatives[i]) or (y == j)]


def nonempty(iterable):
    # Return True if the iterable (generator, list, etc) contains any values.
    # If it's a generator, we obviously want to shortcut as soon as any value
    # is produced, so we don't waste too much time computing the rest of the
    # values. (Necessary if it's an infinite generator!)
    for item in iterable:
        return True
    return False

def match_atomic(particle, descriptor):
    """Is it the same particle type as given by the descriptor?"""
    if descriptor.name == 'X':
        return True
    elif descriptor.name == 'X0':
        return particle.getCharge() == 0
    elif descriptor.name == 'X+':
        return particle.getCharge() > 0
    elif descriptor.name == 'X-':
        return particle.getCharge() < 0
    return (PDB.GetParticle(descriptor.name)
            and PDB.GetParticle(descriptor.name).PdgCode() == particle.getPDG())


def get_mcp_children(particle):
    first = particle.getFirstDaughter()
    last = particle.getLastDaughter()
    if first > 0:
        # Child indices in normal array (starting at 0) indexing.
        return list(range(first - 1, last))
    else:
        return []    
    

def match_decay(mcparticles, index, descriptor):
    """
    Does the decay of mcparticles[index] match the descriptor?
    """
    if not match(mcparticles, index, descriptor.origin):
        # Match the left hand side.
        return False

    if descriptor.arrow in ('->', '=>'):
        # Direct decay.
        mcp_children = get_mcp_children(mcparticles[index])
        mcp_num = len(mcp_children)
        mcp_items = range(mcp_num)

        if descriptor.inclusive:
            # Can have extra monte carlo decay products.
            # That is, none of the monte carlo are 'required' to be matched.
            left = set()
        elif descriptor.arrow == '=>':
            # Can have extra monte carlo gammas.
            # 22 is the PDG code for a photon
            left = set(i for i in mcp_items if mcparticles[mcp_children[i]].getPDG() != 22)
        else:
            left = set(mcp_items)

        desc_num = len(descriptor.decays)
        desc_items = range(desc_num)
        right = set(desc_items)

        # basic requirement
        if not (len(left) <= desc_num <= mcp_num):
            return False

        # Find a way of matching up the mcparticle children to the decay descriptor tree.
        edges = [(i, j) for i in mcp_items for j in desc_items
                 if match(mcparticles, mcp_children[i], descriptor.decays[j])]
        return nonempty(recursive_matching(left, right, edges))
    else:
        # Long decay (--> or ==>).
        mcp_children = get_mcp_children(mcparticles[index])

        # Add descendants of our mcp children, and calculate relatives.
        ancestors = [[] for _ in range(len(mcp_children))]
        descendants = [[] for _ in range(len(mcp_children))]        
        i = 0
        while i < len(mcp_children):
            index = mcp_children[i]
            this_children = get_mcp_children(mcparticles[index])
            if this_children:
                start = len(mcp_children)
                mcp_children.extend(this_children)
                # ...and update everyone's relatives.
                for j in range(start, len(mcp_children)):
                    ancestors.append(ancestors[i] + [i])
                    for k in ancestors[j]:
                        descendants[k].append(j)
                    descendants.append([])
            i += 1
        relatives = [ancestors[i] + [i] + descendants[i] for i in range(len(mcp_children))]
                
        mcp_num = len(mcp_children)
        mcp_items = range(mcp_num)

        if descriptor.inclusive:
            left = set()
        elif descriptor.arrow == '=>':
            left = set(i for i in mcp_items if mcparticles[mcp_children[i]].getPDG() != 22)
        else:
            left = set(mcp_items)

        desc_num = len(descriptor.decays)
        desc_items = range(desc_num)
        right = set(desc_items)

        # Find a way of matching up the mcparticle children to the decay descriptor tree.
        edges = [(i, j) for i in mcp_items for j in desc_items
                 if match(mcparticles, mcp_children[i], descriptor.decays[j])]
        return nonempty(long_matching(left, right, edges, relatives))


def match_logical(mcparticles, index, descriptor):
    if descriptor.op == '||':
        return match(mcparticles, index, descriptor.left) or match(mcparticles, index, descriptor.right)
    elif descriptor.op == '&&':
        return match(mcparticles, index, descriptor.left) and match(mcparticles, index, descriptor.right)
    else:
        raise ValueError("Unknown logical operator '{}'.".format(descriptor.op))

def any(iterable):
    for item in iterable:
        if item:
            return True
    return False

def match_list(mcparticles, index, descriptor):
    return any(match(mcparticles, index, d) for d in descriptor.items)

def match(mcparticles, index, descriptor):
    if isinstance(descriptor, AtomicDescriptor):
        return match_atomic(mcparticles[index], descriptor)
    elif isinstance(descriptor, DecayDescriptor):
        return match_decay(mcparticles, index, descriptor)
    elif isinstance(descriptor, LogicalDescriptor):
        return match_logical(mcparticles, index, descriptor)
    elif isinstance(descriptor, ListDescriptor):
        return match_list(mcparticles, index, descriptor)
    else:
        raise ValueError("Unrecognised particle descriptor.")

class DecayFilterModule(Module):
    def __init__(self):
        Module.__init__(self)
        self.setName('DecayFilterModule')
        self.descriptor = None

    def set_decay(self, descriptor):
        self.descriptor = decode(descriptor)

    def event(self):
        mcparticles = Belle2.PyStoreArray('MCParticles')
        mcp = list(mcparticles)

        for (i, p) in enumerate(mcp):
            if match(mcp, i, self.descriptor):
                print_tree(mcp, i)
