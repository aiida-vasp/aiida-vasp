.. _vasp:

Customizations for VASP developers
==================================

Adding support for additional supported INCAR tags
--------------------------------------------------
Currently we do checks of the tags supplied as parameters. They need to comply with the keys
in the `parameters.yml` file that are place in the `assistant` folder.

We have based this file on the VASP wiki and try to keep this updated.

For VASP developers that have need of additional tags not mentioned in the wiki page of `VASP INCAR tags`_, or if there
are cases where the updates between the wiki and the tags file is out of sync, additional tags can be supplied in the `settings.unsupported_parameters` input of the `VaspWorkChain` (or any workchain calling it). The format of this should be similar to the dictionary representation of `parameters.yml` file that is installed and contains the standard VASP tags as represented on the `VASP INCAR tags`_. You should pass a dictionary of the following form::

  unsupported_parameters = {'mynewtag': {'default': False,
                                         'description': 'A boolean to set enable my new awesome method',
                                         'type': float
                                         'values': [True, False]}}

Currently, only the tag key is used in the code to check that the user supplies a valid tag. If you use custom tags and do not supply this tag in ``parameter_tags``, the calculation will not start.

.. _VASP INCAR tags: https://www.vasp.at/wiki/index.php/Category:INCAR
