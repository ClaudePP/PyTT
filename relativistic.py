# -*- coding: utf-8 -*-

#    Virtual-IPM is a software for simulating IPMs and other related devices.
#    Copyright (C) 2017  The IPMSim collaboration <http://ipmsim.gitlab.io/IPMSim>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy
import scipy.constants as constants


def compute_beta_and_gamma_from_energy(energy, rest_energy):
    """
    Compute relativistic beta and gamma values from particle energy.
    
    Parameters
    ----------
    energy : float
        Kinetic energy of the particle.
    rest_energy : float
        Rest energy of the particle.
        
    Returns
    -------
    tuple
        Tuple containing (beta, gamma).
    """
    gamma = 1. + energy / rest_energy
    beta = numpy.sqrt(1. - 1. / gamma ** 2)
    return beta, gamma


def compute_bunch_length_as_time(bunch_length, energy, rest_energy):
    """
    Compute the bunch length as time from the corresponding value given as length. This includes
    a transformation from the bunch frame to the lab frame.
    
    Parameters
    ----------
    bunch_length : float
        Bunch length in the bunch frame in units of [m].
    energy : float
        Kinetic energy of the particle.
    rest_energy : float
        Rest energy of the particle.
        
    Returns
    -------
    float
        The bunch length in the lab frame as time.
    """
    beta, gamma = compute_beta_and_gamma_from_energy(energy, rest_energy)
    return bunch_length / (beta * constants.speed_of_light * gamma)


def convert_lab_frame_bunch_length_to_bunch_frame(bunch_length_lab, energy, particle_type):
    """
    Computes the bunch length in the bunch's rest frame from the measured length in the lab frame.
    
    Parameters
    ----------
    bunch_length_lab : float
        The bunch length measured in the lab frame in units of [s].
    energy : float
        The bunch energy in units of [eV].
    particle_type : object
        Must have a property `rest_energy` which provides the bunch's particles' rest energy in
        units of [eV].
        
    Returns
    -------
    float
        The bunch length in the bunch's rest frame in units of [m].
    """
    # noinspection PyUnresolvedReferences
    beta, gamma = compute_beta_and_gamma_from_energy(energy, particle_type.rest_energy)
    return bunch_length_lab * beta * constants.speed_of_light * gamma
