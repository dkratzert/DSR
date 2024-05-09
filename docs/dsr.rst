
Preface
=======

The following user manual explains the usage of the program DSR. It
allows for a semi-automatic modelling of disordered moieties. Its
database contains molecular fragments and their corresponding restraints
as well as a fitting procedure to place these fragments on the desired
position in the unit cell.

The web page can be found at `<https://dkratzert.de/dsr.html>`_ and the
development platform at `<https://github.com/dkratzert/DSR>`_.

If you find any bugs in this program, have feature requests or just
comments; please don't hesitate to write an email to <dkratzert@gmx.de>
to report these.

Please cite DSR as:

D. Kratzert, J.J. Holstein, I. Krossing, *J. Appl. Cryst.*
**2015** *48*, 933-938.
`doi:10.1107/S1600576715005580 <http://scripts.iucr.org/cgi-bin/paper?S1600576715005580>`_


**Disclaimer**

You are responsible for the correctness of restraints applied to your
crystal structure. The restraints applied by DSR are only suggestions to
stabilize a first model. **DSR does not replace the judgment of an
expert.**

Program Overview
================

The program-package consists of a simple text-database with fragments of
molecules and the DSR program itself. It acts as a preprocessor for
SHELXL .res files. The user has to insert a special command line in the
SHELXL .res file and the DSR program reads this information. The command
lines main purpose is to tell DSR, how to orient a molecular fragment
from the database in the unit cell. In practice, the user has to choose
a minimum of three target positions (atoms or Q-peaks) in the structure
and the corresponding atoms from the database fragment (source atoms)
that should be placed on the target positions.

You prefer Olex2? Then you might want to use
the `FragmentDB <https://www.xs3.uni-freiburg.de/research/fragmentdb>`_ plugin.
It is a port of DSR to Olex2.


Installation
============

Windows
-------

Execute the "DSR-setup-[version number].exe" and follow the
instructions.
DSR needs at least Windows 7 or Windows 10 for the latest version.
**If you install DSR the first time, a restart of Windows
might be necessary.**

Linux
-----

Either install DSR according to
the installation procedure of your LINUX distribution. Or install the
dsr-shelx package from `pypi.org <https://pypi.org/project/dsr-shelx/>`_.
DSR expects a shelxl or xl executable version 2013 or above in the
system path.

Since version 213, DSR needs the numpy package to be installed in order
to run.

The update from DSR for python 2.7 (version 236 and below) to DSR for Python 3
needs some manual work. You need to uninstall the old version and make sure the
files dsr.sh, dsr-mac and dsr-linux are removed from /etc/profile.d/ and/or
/etc/paths.d/.

MacOS
-----

Install the dsr-shelx package from `pypi.org <https://pypi.org/project/dsr-shelx/>`_.

The update from DSR for python 2.7 (version 236 and below) to DSR for Python 3
needs some manual work. You need to uninstall the old version and make sure the
files dsr.sh, dsr-mac and dsr-linux are removed from /etc/profile.d/ and/or
/etc/paths.d/.


User-defined database
---------------------

DSR expects the database file "dsr_user_db.txt" with
self-made fragments in the user's home directory. The installation
procedures create an empty user database in the following directories if
no previous database exists:

Windows: C:\\users\\username\\

Linux: /home/username/

Mac OS X: /Users/username/


DSR in ShelXle
==============

Since version 181, `ShelXle <http://www.shelxle.org>`_
has the ability to start a graphical user interface for DSR. A mouse
click on *Tools*--\> *DSR plugin* will start the DSR GUI:

![](/Users/daniel/Documents/GitHub/DSR/docs/images/media/image1.png){width="6.166666666666667in"
height="0.6959995625546807in"}

![](/Users/daniel/Documents/GitHub/DSR/docs/images/media/image2.tiff){width="4.0779921259842515in"
height="3.6771489501312336in"}

Now you need to select a fragment in the list. The list of fragments can
be shortened via the search field.

To fit a fragment into the structure in ShelXle, select three
atoms/Q-peaks in the target molecule (with a left mouse click while
holding STRG) and the fragments 3D view (just left mouse click) each.
The 3D view should now show a preview of the fitted fragment:

![](/Users/daniel/Documents/GitHub/DSR/docs/images/media/image3.png){width="2.4102449693788275in"
height="2.256988188976378in"}

You can now control all the features of DSR with the options menu below:

![](/Users/daniel/Documents/GitHub/DSR/docs/images/media/image4.tiff){width="6.071395450568679in"
height="1.2321544181977253in"}

Setting PART to zero will disable them. The residue number will always
be chosen as the next free available. You can safely leave this as it is
or change the residue name.

The "Free variable" option defines the free variable for the fragment
occupation in SHELXL. The Free variable will be combined with the
occupation option. For example a free variable of -3 and an occupation
of 1 will be combined to -31. The result appears instantly in the output
window.

"External restraints" writes the restraints to an external file.

"Calculate DFIX" automatically generates DFIX/DANG/FLAT restraints from
the geometry of the fragment.

![](/Users/daniel/Documents/GitHub/DSR/docs/images/media/image5.png){width="4.8909492563429575in"
height="3.4003762029746283in"}

To create or edit a fragment, click on \"Edit fragment\". The edit
window allows adding, updating and deleting of fragments.

Similar to the syntax in \"dsr_usr_db.txt\", you can choose to define
the atom type by the name of the atom or with a negative atomic number.

They will be stored in the users fragment database \"dsr_usr_db.txt\" in
your home directory. Different to the fragment creation by hand, you do
not have to invent a database name tag. It will be randomly chosen,
because the GUI will never show them. Instead the GUI always shows the
real fragments names.

If you have a new fragment, you should consider sending it to me by
clicking on \"Mail Fragment Home\".
