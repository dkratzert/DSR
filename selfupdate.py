# -*- encoding: utf-8 -*-
# möp
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <daniel.kratzert@ac.uni-freiburg.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
from __future__ import print_function

def get_current_dsr_version():
    """
    determines the current version of DSR on the web server

    >>> get_current_dsr_version()
    123

    Returns
    -------
    version number
    :type: int
    """
    pass

def update_dsr():
    """
    Updates the running DSR to the current version on the web server.

    >>> update_dsr()
    True

    Returns
    -------
    Sucess or not
    True/False
    """
    pass

def is_update_needed():
    """
    Decides if an update of DSR is needed and updates if neccesary. It does
    nothing in case of an already updated version.

    >>> is_update_needed()
    "updating DSR ..."
    "update was sucessful"
    "You are already using the current version. No update needed."

    Returns
    -------
    True/False
    """
    pass


def get_update_package():
    """
    Downloads the current DSR distribution from the web server and
    returns True if it suceeded.

    Returns
    -------
    True/False
    """
    pass

