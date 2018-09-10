# Changelog
All notable changes to this project after version 0.1.1 will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [v0.2.4]

### Added
- 'vasp.base' workchain which transparently calls through to the 'vasp.vasp' calculation and restarts if appropriate.
 - restarting conditions are minimal yet (only submission failed will trigger a retry)
- 'vasp.relax' workchain which specializes in structure relaxation and optionally iterates until the structure converges.
- 'vasp.calc.immigrant' added in order to support VASP import Aiida-external jobs.

### Changed
- PotcarData.get_potcars_from_structure() now returns an entry for each `kind.name` in the structure, not one per `kind.symbol`
- PotcarData.get_potcars_dict() no longer fails if there is more than one PotcarData with the same full name in the family
- Minor changes in parser.

## [v0.2.3]

### Changed
- fixed POSCAR writing (was writing 'direct' followed by cartesian coordinates)

### Removed
- pruned outdated parts of documentation

## [v0.2.2] - 2018-03-15

### Changed
- missing requirement added
- PyPI description and keywords updated

## [v0.2.1] - 2018-03-12

### Added
- Allow Structures which do not have sites of same element grouped
- Allow using different potentials for atoms of same element in a structure
- `VaspCalculation`'s `settings` input now accepts `poscar_precision` key to set maximum precision for coordinates in POSCAR
 - This can be indicated in the structure by adding sites with same `symbol` but different `name`.
- `io.poscar.PoscarIo`, POSCAR writer replacing pymatgen one, to be replaced by parsevasp.

### Changed
- Use always the same order for species in POSCAR and POTCAR
- POSCAR float precision default changed to 10 decimal places

## [v0.2.0] - 2018-03-07

### Added
- `data.potcar.PotcarData` (`vasp.potcar`) & `data.potcar.PotcarFileData`, replacement for PawData
 - `PotcarData` is shareable, holds no licenced data
 - `PotcarFileData` holds the licenced POTCAR file
- `io.potcar.PotcarIo`, handle conversion between File and DB representation of POTCAR information
- `io.potcar.MultiPotcarIo`, handle concatenation and splitting of POTCAR files containing multiple potentials
- `verdi data vasp-potcar`, CLI for importing and exporting POTCAR files
 - `uploadfamily`, loads a set of POTCAR files into the db
 - `exportfamily`, creates a compressed archive of POTCAR files in a family
 - `listfamilies`, list existing POTCAR families
- `io.outcar.OutcarParser`, parses some information from OUTCAR file
- `calcs.VaspCalculation` and descendants:
 - inputs: `settings:ParameterData` additional settings that are not passed to VASP but to parsers etc

### Changed
- `calcs.VaspCalculation` and descendants:
 - inputs: `paw:PawData` -> `potentials:PotcarData`
- `parsers.vasp.VaspParser`, redesigned:
 - parses OUTCAR too (or tries to)
 - takes `parser_options` that can be given in the `settings` input in a `VaspCalculation`
- renamed `data.paw.PawData` -> `data.paw.LegacyPawData`
 - `LegacyPawData` can not be stored or changed (read-only)

### Removed
- `io.potcar.PawParser`

### Deprecated
- `PawData` is deprecated and has been renamed `LegacyPawData`, it can be used to read `PawData` database objects.

## [v0.1.1]

Baseline
