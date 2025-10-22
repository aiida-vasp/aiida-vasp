"""
VASP Calculation Monitoring Functions.

This module provides monitoring functions for detecting and handling problems that may occur
during VASP calculations running on remote machines. These monitors are designed to identify
issues early and take preventive actions to avoid system crashes or resource waste.

The monitoring functions can detect:

* Stdout file overflow that could crash the AiiDA daemon
* Electronic loop timing issues that indicate inefficient calculations
* Stalled calculations that are no longer making progress

These monitors are typically used by AiiDA's calculation monitoring system to provide
real-time feedback about running VASP calculations and automatically handle problematic
situations.

.. note::
   These monitors operate on the remote machine where the VASP calculation is running
   and require a transport connection to access files and execute commands.
"""

import time
from pathlib import Path

from aiida.orm import CalcJobNode
from aiida.transports import Transport


def monitor_stdout(node: CalcJobNode, transport: Transport, size_threshold_mb: float = 5) -> str | None:
    """
    Monitor stdout file size to prevent overflow crashes.

    This function monitors the size of the VASP stdout file during calculation execution.
    If the file becomes too large, it indicates a potential problem (such as excessive
    output from convergence issues) that could crash the AiiDA daemon when attempting
    to retrieve and parse the file.

    When an oversized stdout file is detected, the function automatically truncates it
    to prevent system crashes, though this means the calculation is considered lost.

    :param node: The CalcJobNode representing the running VASP calculation
    :type node: CalcJobNode
    :param transport: Transport connection to the remote machine where VASP is running
    :type transport: Transport
    :param size_threshold_mb: Maximum allowed stdout file size in megabytes before
                             truncation occurs
    :type size_threshold_mb: int
    :returns: None if no overflow detected, otherwise an error message describing
              the overflow condition
    :rtype: str or None

    .. warning::
       When stdout overflow is detected, the calculation is automatically terminated
       by truncating the output file. This prevents system crashes but results in
       loss of the calculation.

    .. note::
       The default threshold of 5 MB is typically sufficient for normal VASP
       calculations. Larger thresholds may be needed for very large systems or
       calculations with verbose output.
    """

    # Check the current size of the VASP stdout file
    stdout_path = str(Path(node.get_remote_workdir()) / node.process_class._VASP_OUTPUT)
    file_stat = transport.get_attribute(stdout_path)
    if file_stat.st_size > 1024 * 1024 * size_threshold_mb:
        # Stdout file is dangerously large - truncate it to prevent system crashes
        # This typically indicates convergence problems or infinite loops in VASP
        # The calculation is lost, but we prevent the AiiDA daemon from crashing
        # when attempting to retrieve and parse the oversized file
        transport.exec_command_wait(f'truncate -s {size_threshold_mb}M {stdout_path}')
        return f'Very large stdout detected: {file_stat.st_size / (1024 * 1024):.2f} MB, potential critical crash.'
    # TODO: add detection of critical messages emitted when vasp got stuck


def monitor_loop_time(
    node: CalcJobNode,
    transport: Transport,
    minimum_electronic_loops: int = 10,
    patience_num_loops: int = 5,
    patience_minimum_time: float = 1800,
) -> str | None:
    """
    Monitor electronic loop timing to detect inefficient or stalled calculations.

    This function analyzes the timing of electronic self-consistency loops in VASP
    calculations to identify potential problems:

    1. **Slow convergence**: If electronic loops take too long relative to the walltime
       limit, the calculation may not complete within the allocated time.

    2. **Stalled calculations**: If the stdout file hasn't been updated for an
       extended period, the calculation may have crashed or become stuck.

    The function examines the OUTCAR file to extract loop timing information and
    compares it against the walltime limits and recent file modification times.

    :param node: The CalcJobNode representing the running VASP calculation
    :type node: CalcJobNode
    :param transport: Transport connection to the remote machine where VASP is running
    :type transport: Transport
    :param minimum_electronic_loops: Minimum number of electronic loops that should
                                   be completable within the walltime limit
    :type minimum_electronic_loops: int
    :param patience_num_loops: Number of loop times to wait before considering
                             a calculation stalled
    :type patience_num_loops: int
    :param patience_minimum_time: Minimum time in seconds to wait before checking
                                for stalled calculations
    :type patience_minimum_time: int
    :returns: None if timing is acceptable, otherwise an error message describing
              the detected problem (slow loops or stalled calculation)
    :rtype: str or None
    """

    outcar_path = str(Path(node.get_remote_workdir(), 'OURCAR'))
    stdout_path = str(Path(node.get_remote_workdir()) / node.process_class._VASP_OUTPUT)
    walltime_limit = node.get_option('max_wallclock_seconds')
    if walltime_limit is None:
        # Cannot monitor timing without a walltime limit
        return None

    # Extract electronic loop timings from OUTCAR file
    # LOOP entries contain timing information for each self-consistency cycle
    returncode, stdout, _ = transport.exec_command_wait(f"grep 'LOOP:' {outcar_path}")

    # Skip monitoring if no LOOP entries found (calculation hasn't started electronic steps)
    if returncode != 0:
        return None

    # Parse the last loop time from the grep output
    # Each LOOP line format: "LOOP:  CPU time  real time  (sec)  real time"
    last_loop_time = 0
    for line in stdout.splitlines():
        # Extract the real time (4th column) from each LOOP line
        # Keep updating to get the most recent loop time
        last_loop_time = float(line.strip().split()[-1])

    # Check if electronic loops are too slow relative to walltime
    # If each loop takes more than walltime/minimum_loops, we won't finish in time

    # Check for the presence of NELM
    incar = {key.lower(): value for key, value in node.inputs.get('parameters', {}).items()}
    nelm = incar.get('nelm', 60)
    # Remove none-scf loops for hybrid calculation
    if 'lhfcalc' in incar:
        nelm -= abs(incar.get('nelmdl', 5))

    # Take the minimum of the nelm and supplied minimum_electronic_loops
    # This is useful to avoid killing benchmarking calculations a low NELM setting
    minimum_electronic_loops = max(1, min(minimum_electronic_loops, nelm))

    if last_loop_time > walltime_limit / minimum_electronic_loops:
        return (
            f'Less than {minimum_electronic_loops} electronic loop can be run in this calculation due to '
            f'long electronic loop time: {last_loop_time:.2f} seconds.'
        )

    # Monitor for stalled calculations by checking stdout file modification time
    file_stat = transport.get_attribute(stdout_path)
    elapsed = time.time() - file_stat.st_mtime
    # Consider calculation stalled if:
    # 1. Minimum patience time has elapsed AND
    # 2. Time since last update exceeds expected time for several loops
    if elapsed > patience_minimum_time and elapsed > last_loop_time * patience_num_loops:
        return (
            f'Last update of the {stdout_path} file is more than {last_loop_time * patience_num_loops:.2f} '
            'seconds ago. It is very likely that the calculation is stalled or crashed.'
        )
