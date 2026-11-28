import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

from dsr_shelx.dbfile import ParseDB
from dsr_shelx.restraints import Restraints, ListFile, Lst_Deviations, format_atom_names

REPO_DIR = Path(__file__).parent.parent


def _two_ring_db(tmp_path):
    """
    A cyclohexane chair (C1-C6, not flat) and a benzene ring (C7-C12, flat).
    The chair is the first ring that make_flat_restraints() looks at.
    """
    atoms = []
    for n in range(6):
        angle = math.radians(60 * n)
        z = 0.25 if n % 2 else -0.25
        atoms.append('C{} 1 {:.5f} {:.5f} {:.5f}'.format(n + 1, 1.446 * math.cos(angle), 1.446 * math.sin(angle), z))
    for n in range(6):
        angle = math.radians(60 * n)
        atoms.append('C{} 1 {:.5f} {:.5f} 0.0'.format(n + 7, 10 + 1.39 * math.cos(angle), 1.39 * math.sin(angle)))
    entry = '<tworings>\nREM Name: Two rings\nRESI TWO\nSIMU C1 > C12\nFRAG 17 1 1 1 90 90 90\n{}\n</tworings>\n'
    dbpath = tmp_path / 'db.txt'
    dbpath.write_text(entry.format('\n'.join(atoms)))
    return ParseDB(maindb_path=str(dbpath))


def test_flat_restraints_after_a_non_planar_ring(tmp_path):
    flats = Restraints('tworings', _two_ring_db(tmp_path)).make_flat_restraints()
    assert flats
    benzene = {'C{}'.format(n) for n in range(7, 13)}
    for flat in flats:
        assert set(flat) <= benzene


def test_flat_restraints_do_not_depend_on_hash_seed():
    # String hashes and therefore set order change with PYTHONHASHSEED:
    code = ('import sys; sys.path.insert(0, "src"); import dsr_shelx; from dbfile import ParseDB; '
            'from restraints import Restraints; gdb = ParseDB("src/dsr_shelx/dsr_db.txt"); '
            'print([Restraints(f, gdb).get_formated_flats() for f in ("pyrene", "c60", "dabco", "try")])')
    results = set()
    for seed in ('1', '2', '3'):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        results.add(subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True,
                                   check=True, cwd=REPO_DIR).stdout)
    assert len(results) == 1


LST = """ CELL  0.71073  10.5086  20.9035  20.5072  90.000  94.130  90.000

 Covalent radii and connectivity table for test
 C    0.770
 F    0.640

C1 -  C2 C4
C2 -  C1 F1
F1 -  C2
C4 -  C1
{clockwise}
     ATOM       x         y         z
 
 C1      0.22379   0.37551   0.50415
 F1_2a   0.30000   0.40000   0.50000

  ** Atom C2 deviates by 0.512 A
"""


@pytest.mark.parametrize('clockwise, bondvector', [
    ('', ('C1', 'C2')),
    (' Rotation is clockwise looking down C4 to C1. The first value is for F1\n', ('C4', 'C1')),
])
def test_list_file(tmp_path, clockwise, bondvector):
    (tmp_path / 'test.lst').write_text(LST.format(clockwise=clockwise))
    lf = ListFile(str(tmp_path / 'test'))
    assert lf.get_cell() == [10.5086, 20.9035, 20.5072, 90.0, 94.13, 90.0]
    assert lf.get_all_coordinates == {'C1': [0.22379, 0.37551, 0.50415], 'F1_2A': [0.3, 0.4, 0.5]}
    assert lf.get_single_coordinate('F9') is None
    assert ('C2', 'F1') in lf.read_conntable()
    assert lf.get_bondvector('C2') == bondvector
    assert Lst_Deviations(lf.listfile_list).find_deviations() == {'C2': '0.512'}


@pytest.mark.parametrize('part, resinum, names', [
    ('2', '4', ['C1_4B', 'O1_4B']),
    ('1', '', ['C1_A', 'O1_A']),
    ('', '3', ['C1_3', 'O1_3']),
    ('', '', ['C1', 'O1']),
])
def test_format_atom_names(part, resinum, names):
    assert format_atom_names(['C1', 'O1'], part, resinum) == names
