# Changelog
All notable changes to this project after version 0.1.1 will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
