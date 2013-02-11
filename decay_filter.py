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
        prefix = '\033[3{}m'.format(level)
        suffix = '\033[m'
    else:
        prefix = ''
        suffix = ''
    print '    ' * level, prefix + '[{}] {} mass={} energy={} charge={}'.format(pid, name, mass, energy, charge) + suffix
    if p.getFirstDaughter() > 0:
        m = p.getFirstDaughter()
        n = p.getLastDaughter()
        for j in range(m-1, n):
            print_tree(mcparticles, j, level + 1)

def recursive_matching(left, right, edges, inclusive=set()):
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
    if (not left.difference(inclusive)) and (not right):
        # Everything that needs to be is matched.
        yield set()
        return
    if not edges:
        # Something is unmatched, and we're out of free edges.
        return

    incompatible = list(edges)
    while incompatible:
        (i, j) = incompatible.pop()
        new_left = left.copy()
        new_left.remove(i)
        new_right = right.copy()
        new_right.remove(j)
        new_edges = {(x, y) for (x, y) in edges if (x != i) and (y != j)}
        for matching in recursive_matching(new_left, new_right, new_edges, inclusive):
            yield matching.union({(i, j)})
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
        new_edges = {(x, y) for (x, y) in edges if (x not in relatives[i]) and (y != j)}
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

def match_decay(mcparticles, index, descriptor):
    """
    Does the decay of mcparticles[index] match the descriptor?
    """
    if not match(mcparticles, index, descriptor.origin):
        # Match the left hand side.
        return False

    first = mcparticles[index].getFirstDaughter()
    last = mcparticles[index].getLastDaughter()
    if first > 0:
        # Child indices in normal array (starting at 0) indexing.
        mcp_decays = list(range(first - 1, last))
        num_mcp = len(mcp_decays)
    else:
        mcp_decays = []
        num_mcp = 0

    if descriptor.inclusive:
        inclusion = 'all'
    elif descriptor.arrow in ('=>', '==>'):
        inclusion = 'gamma'
    else:
        inclusion = 'none'

    if descriptor.arrow in ('-->', '==>'):
        return long_decay(mcparticles, mcp_decays, descriptor.decays, inclusion)

    if inclusion == 'all':
        # Can have extra monte carlo decay products.
        # That is, any decay products are allowed to be unmatched.
        inclusive = set(range(num_mcp))
    elif inclusion == 'gamma':
        # Can have extra monte carlo gammas.
        # 22 is the PDG code for a photon
        inclusive = set(i for i in range(num_mcp)
                        if mcparticles[mcp_decays[i]].getPDG() == 22)
    else:
        inclusive = set()

    num_desc = len(descriptor.decays)
    if not (num_mcp - len(inclusive) <= num_desc <= num_mcp):
        return False

    # Find a way of matching up the mcparticle children to the decay descriptor tree.
    edges = {(i, j) for i in range(num_mcp)
                    for j in range(num_desc)
                    if match(mcparticles, mcp_decays[i], descriptor.decays[j])}
    return nonempty(recursive_matching(set(range(num_mcp)), set(range(num_desc)), edges, inclusive))

def long_decay(mcparticles, mcp_decays, desc_decays, inclusion):
    descendants = {i : [] for i in range(len(mcp_decays))}
    ancestors = {i : [] for i in range(len(mcp_decays))}

    # we're going to do the match all at once on the decay tree
    i = 0
    while i < len(mcp_decays):
        index = mcp_decays[i]
        first = mcparticles[index].getFirstDaughter()
        if first > 0:
            start = len(mcp_decays)
            # Has children, so let's expand!
            last = mcparticles[index].getLastDaughter()
            thischildren = list(range(first - 1, last))

            # Add the children to the end...
            mcp_decays.extend(thischildren)

            # ...and update everyone's relatives.
            for j in range(start, start + len(thischildren)):
                ancestors[j] = ancestors[i] + [i]
                for k in ancestors[j]:
                    descendants[k].append(j)
                descendants[j] = []
        i += 1

    relatives = {i : [i] + descendants[i] + ancestors[i] for i in range(len(mcp_decays))}

    if inclusion == 'all':
        # Can have extra monte carlo decay products.
        # That is, none of the monte carlo are 'required' to be matched.
        left = set()
    elif inclusion == 'gamma':
        # Can have extra monte carlo gammas.
        # 22 is the PDG code for a photon
        left = set(i for i in range(len(mcp_decays)) if mcparticles[mcp_decays[i]].getPDG() != 22)
    else:
        left = set(range(len(mcp_decays)))

    right = set(range(len(desc_decays)))

    edges = {(i, j) for i in range(len(mcp_decays))
                    for j in range(len(desc_decays))
                    if match(mcparticles, mcp_decays[i], desc_decays[j])}

    return nonempty(long_matching(left, right, edges, relatives))

def match_logical(mcparticles, index, descriptor):
    if descriptor.op == '||':
        return match(mcparticles, index, descriptor.left) or match(mcparticles, index, descriptor.right)
    elif descriptor.op == '&&':
        return match(mcparticles, index, descriptor.left) and match(mcparticles, index, descriptor.right)
    else:
        raise ValueError("Unknown logical operator '{}'.".format(descriptor.op))

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
