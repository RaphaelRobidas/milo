#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Handle calling ESPs and parsing output."""

import os

from milo_1_0_3 import containers
from milo_1_0_3 import enumerations as enums
from milo_1_0_3 import exceptions


def get_program_handler(program_state):
    """Return the configured electronic structure program handler."""
    if program_state.program_id is enums.ProgramID.GAUSSIAN_16:
        return Gaussian16Handler
    elif program_state.program_id is enums.ProgramID.GAUSSIAN_09:
        return Gaussian09Handler
    elif program_state.program_id is enums.ProgramID.ORCA_5:
        return OrcaHandler
    else:
        raise ValueError(f'Unknown electronic structure program '
                         f'"{program_state.program_id}"')


class GaussianHandler:
    """Template handler for Gaussian16Handler and Gaussian09Handler."""

    gaussian_command = ""

    @classmethod
    def generate_forces(cls, program_state):
        """Preform computation and append forces to list in program state."""
        route_section = f"# force {program_state.gaussian_header}"
        log_file = cls.call_gaussian(route_section,
                                     f"{cls.gaussian_command}"
                                     f"_{program_state.current_step}",
                                     program_state)
        cls.parse_forces(log_file, program_state)

    @classmethod
    def call_gaussian(cls, route_section, job_name, program_state):
        """Call Gaussian and return a string with the name of the log file."""
        job_com_file = f"{job_name}.com"
        job_log_file = f"{job_name}.log"
        cls.prepare_com_file(job_com_file, route_section, program_state)
        os.system(f"{cls.gaussian_command} < {job_com_file} > {job_log_file}")
        return job_log_file

    @staticmethod
    def prepare_com_file(file_name, route_section, program_state):
        """Prepare a .com file for a Gaussian run."""
        with open(file_name, 'w') as com_file:
            if program_state.processor_count is not None:
                com_file.write(f"%nprocshared="
                               f"{program_state.processor_count}\n")
            if program_state.memory_amount is not None:
                com_file.write(f"%mem="
                               f"{program_state.memory_amount}gb\n")
            com_file.write(f"{route_section}\n\n")
            com_file.write(f"Calculation for time step: "
                           f"{program_state.current_step}\n\n")
            com_file.write(f" {program_state.charge}"
                           f" {program_state.spin}\n")
            for atom, (x, y, z) in zip(program_state.atoms,
                                       program_state.structures[-1]
                                       .as_angstrom()):
                com_file.write(f"  {atom.symbol} {x:10.6f} {y:10.6f} "
                               f"{z:10.6f}\n")
            com_file.write("\n")
            if program_state.gaussian_footer is not None:
                com_file.write(program_state.gaussian_footer)
            com_file.write("\n\n")

    @classmethod
    def parse_forces(cls, log_file_name, program_state):
        """Parse forces into program_state from the given log file."""
        if not cls.is_log_good(log_file_name):
            raise exceptions.ElectronicStructureProgramError(
                "Gaussian force calculation log file was not valid. Gaussian "
                "returned an error or could not be called correctly.")
        forces = containers.Forces()
        energy = containers.Energies()
        with open(log_file_name) as log_file:
            for line in log_file:
                if "SCF Done" in line:
                    value = float(line.split()[4])
                    energy.append(value, enums.EnergyUnits.HARTREE)
                    program_state.energies.append(energy)
                if "Forces (Hartrees/Bohr)" in line:
                    for data_line in log_file:
                        if "Cartesian Forces" in data_line:
                            program_state.forces.append(forces)
                            return
                        else:
                            tokens = data_line.split()
                            try:
                                # This will throw an error if the first item
                                # on the line is not an integer, and the line
                                # will be skipped.
                                int(tokens[0])
                            except ValueError:
                                continue
                            else:
                                x = float(tokens[2])
                                y = float(tokens[3])
                                z = float(tokens[4])
                                forces.append(x, y, z,
                                    enums.ForceUnits.HARTREE_PER_BOHR)

    @staticmethod
    def is_log_good(log_file_name):
        """Return true if the given log file terminated normally."""
        with open(log_file_name) as log_file:
            for line in log_file:
                if "Normal termination" in line:
                    return True
            else:
                return False


class Gaussian16Handler(GaussianHandler):
    """Call Gaussian16 and parse output."""

    gaussian_command = "g16"


class Gaussian09Handler(GaussianHandler):
    """Call Gaussian09 and parse output."""

    gaussian_command = "g09"


class OrcaHandler:
    """Handler for ORCA5"""

    @classmethod
    def generate_forces(cls, program_state):
        """Preform computation and append forces to list in program state."""
        route_section = f"ENGRAD {program_state.gaussian_header}"
        log_file = cls.call_orca(route_section,
                                     f"orca_{program_state.current_step}",
                                     program_state)
        cls.parse_forces(log_file, program_state)

    @classmethod
    def call_orca(cls, route_section, job_name, program_state):
        """Call Orca and return a string with the name of the log file."""
        job_com_file = f"{job_name}.inp"
        job_log_file = f"{job_name}.out"
        cls.prepare_com_file(job_com_file, route_section, program_state)
        os.system(f"{program_state.orca_path}/orca {job_com_file} > {job_log_file}")
        return job_log_file

    @staticmethod
    def prepare_com_file(file_name, route_section, program_state):
        """Prepare a .inp file for an Orca."""
        with open(file_name, 'w') as com_file:
            com_file.write(f"!{route_section}\n\n")
            if program_state.processor_count is not None:
                com_file.write(f"%pal\n"
                               f"nproc {program_state.processor_count}\n"
                               f"end\n")
                if program_state.memory_amount is not None:
                    com_file.write(f"%MaxCore {program_state.memory_amount*1024//program_state.processor_count}\n")

            com_file.write(f"*xyz {program_state.charge} {program_state.spin}\n")
            for atom, (x, y, z) in zip(program_state.atoms,
                                       program_state.structures[-1]
                                       .as_angstrom()):
                com_file.write(f"  {atom.symbol} {x:10.6f} {y:10.6f} "
                               f"{z:10.6f}\n")
            com_file.write("*\n")
            if program_state.gaussian_footer is not None:
                com_file.write(program_state.gaussian_footer)
            com_file.write("\n\n")

    @classmethod
    def parse_forces(cls, log_file_name, program_state):
        """Parse forces into program_state from the given log file."""

        with open(log_file_name) as f:
            lines = f.readlines()
            ind = len(lines)-1

        while lines[ind].find("CARTESIAN GRADIENT") == -1:
            ind -= 1

        max_norm = 0

        forces = containers.Forces()
        energy = containers.Energies()

        ind += 3
        while lines[ind].strip() != "":
            sline = lines[ind].split()
            forces.append(*[-float(i) for i in sline[3:]], enums.ForceUnits.HARTREE_PER_BOHR)
            ind += 1

        while lines[ind].find("FINAL SINGLE POINT ENERGY") == -1:
            ind -= 1
        E = float(lines[ind].split()[-1])
        energy.append(E, enums.EnergyUnits.HARTREE)
        program_state.energies.append(energy)
        program_state.forces.append(forces)
