import math
import os
import subprocess
import sys
from pathlib import Path

from dsr_shelx.dbfile import ParseDB
from dsr_shelx.restraints import Restraints

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
