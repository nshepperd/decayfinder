#!/usr/bin/env python2
# encoding: utf-8
from collections import namedtuple
from pyparsing import Word, Regex, Or, Literal
from pyparsing import Optional, OneOrMore, ZeroOrMore, Forward
from pyparsing import Group, Suppress

# from ROOT import TDatabasePDG
# PDB = TDatabasePDG.Instance()
# # PDB.ReadPDGTable()
# print PDB.GetParticle("anti-D0")
# names = [item.GetName() for item in PDB.ParticleList()]
# exit()

pname_letters = r'[a-zA-Z0-9/*+_\-]'
pname_resonance = r'(?:\({letters}+\))'.format(letters=pname_letters)
pname = Regex(r'[a-zA-Z]{letters}*{resonance}?'.format(letters=pname_letters, resonance=pname_resonance))

# pname = Or([Literal(name) for name in names])

atomic = (pname
          ^ Literal('X')
          ^ Literal('X0')
          ^ Literal('X+')
          ^ Literal('X-'))

arrow = Literal('->') ^ Literal('=>') ^ Literal('-->') ^ Literal('==>')
operator = Literal('||') ^ Literal('&&')

expression = Forward()
expression << Group(atomic
                    ^ (Suppress('(') + expression + arrow + OneOrMore(expression) + Optional('...') + Suppress(')'))
                    ^ (Suppress('(') + expression + operator + expression + Suppress(')'))
                    ^ (Literal('[') + expression + ZeroOrMore(Suppress(",") + expression) + Literal(']')))

top = (atomic
       ^ (expression + arrow + OneOrMore(expression) + Optional('...'))
       ^ (expression + operator + expression)
       ^ (Literal('[') + expression + ZeroOrMore(Suppress(",") + expression) + Literal(']')))

AtomicDescriptor = namedtuple('AtomicDescriptor', ['name'])
DecayDescriptor = namedtuple('DecayDescriptor', ['origin', 'decays', 'arrow', 'inclusive'])
LogicalDescriptor = namedtuple('LogicalDescriptor', ['op', 'left', 'right'])
ListDescriptor = namedtuple('ListDescriptor', ['items'])

# pyparsing outputs a nested list of symbols,
# with level of nesting determined by Group elements.
# We need to turn this into Descriptor objects.
def decode2(parse):
    if len(parse) == 1:
        # Atomic predicate.
        assert isinstance(parse[0], str) 
        return AtomicDescriptor(name=parse[0])
    else:
        assert len(parse) >= 3
        if parse[1] in ('->', '=>', '-->', '==>'):
            # Decay predicate.
            inclusive = False
            if parse[-1] == '...':
                del parse[-1]
                inclusive = True
            return DecayDescriptor(origin=decode2(parse[0]),
                                   decays=[decode2(x) for x in parse[2:]],
                                   arrow=parse[1],
                                   inclusive=inclusive)
        elif parse[1] in ('||', '&&'):
            # Logical operator.
            return LogicalDescriptor(op=parse[1],
                                     left=decode2(parse[0]),
                                     right=decode2(parse[2]))
        elif parse[0] == '[':
            # List, equivalent to multiple ||.
            return ListDescriptor([decode2(x) for x in parse[1:-1]])

def decode(string):
    return decode2(top.parseString(string, parseAll=True).asList())

def _to_json(desc):
    if isinstance(desc, AtomicDescriptor):
        return {'type' : 'Atomic', 'name' : desc.name}
    elif isinstance(desc, DecayDescriptor): # namedtuple('DecayDescriptor', ['origin', 'decays', 'arrow', 'inclusive'])
        return {'type' : 'Decay',
                'origin' : _to_json(desc.origin),
                'decays' : [_to_json(x) for x in desc.decays],
                'arrow' : desc.arrow,
                'inclusive' : desc.inclusive}
    elif isinstance(desc, LogicalDescriptor): # namedtuple('LogicalDescriptor', ['op', 'left', 'right'])
        return {'type' : 'Logical',
                'op' : desc.op,
                'left' : _to_json(desc.left),
                'right' : _to_json(desc.right)}
    elif isinstance(desc, ListDescriptor): # = namedtuple('ListDescriptor', ['items'])
        return {'type' : 'List',
                'items' : [_to_json(x) for x in desc.items]}
    else:
        raise ValueError("Not a descriptor.")

def to_json(desc):
    import json
    return json.dumps(_to_json(desc))

# print decode("(B0 -> X+ X- ...) || (B0 -> D*0_bar ...)")
