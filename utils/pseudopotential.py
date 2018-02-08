from aiida_vasp.data.paw import PawData, get_pseudos_from_structure

def validate_and_prepare_pseudos_inputs(structure, pseudos=None, pseudo_family=None):
    """
    Use the explicitly passed pseudos dictionary or use the pseudo_family in combination with
    the structure to obtain that dictionary.
    The pseudos dictionary should now be a dictionary of PAW nodes with the kind as linkname
    As such, if there are multiple kinds with the same element, there will be duplicate PAW nodes
    but multiple links for the same input node are not allowed. Moreover, to couple the PAW nodes
    to the Calculation instance, we have to go through the use_pseudo method, which takes the kind
    name as an additional parameter. When creating a Calculation through a Process instance, one
    cannot call the use methods directly but rather should pass them as keyword arguments. However, 
    we can pass the additional parameters by using them as the keys of a dictionary
    :param structure: StructureData node
    :param pseudos: a dictionary where keys are tuples of kind name and value are PawData nodes
    :param pseudo_family: string name of the pseudopotential family to use
    :raises: ValueError if neither pseudos or pseudo_family is specified or if no PawData is found for
        every element in the structure
    :returns: a dictionary of PawData nodes where the key is a tuple with the kind name
    """
    result_pseudos = {}
    unique_pseudos = {}

    if pseudos is None and pseudo_family is None:
        raise ValueError('neither an explicit pseudos dictionary nor a pseudo_family was specified')
    elif pseudo_family:
        pseudos = get_pseudos_from_structure(structure, pseudo_family.value)

    for kind in structure.get_kind_names():
        if (kind,) not in pseudos:
            raise ValueError('no pseudo available for element {}'.format(kind))
        elif not isinstance(pseudos[ (kind,) ], PawData):
            raise ValueError('pseudo for element {} is not of type PawData'.format(kind))

#    for kind, pseudo in pseudos.iteritems():
#        unique_pseudos.setdefault(pseudo, []).append(kind)

#    for pseudo, kinds in unique_pseudos.iteritems():
#        result_pseudos[(kinds,)] = pseudo

    return pseudos
