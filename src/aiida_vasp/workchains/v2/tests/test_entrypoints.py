from aiida.plugins import WorkflowFactory


def test_entrypoints(aiida_profile):
    """Test that we can instantiate the workchains from the entry points."""
    _ = aiida_profile
    entrypoints = [
        'vasp.vasp', 'vasp.relax', 'vasp.converge', 'vasp.bands', 'vasp.v2.vasp', 'vasp.v2.relax', 'vasp.v2.converge',
        'vasp.v2.bands', 'vasp.v2.hybrid_bands'
    ]
    for point in entrypoints:
        WorkflowFactory(point)
