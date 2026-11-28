import re
from types import SimpleNamespace

import pytest

from dsr_shelx.atomhandling import FindAtoms
from tests.test_dsr import run_dsr, TESTS_DIR

# Fractional coordinates from p21c.res:
COORDS = {'C1': [0.223789, 0.375512, 0.504147], 'C4': [0.281319, 0.368769, 0.575106]}


def fake_refine_cycle(self, rl, reslist):
    """Replaces the SHELXL run: writes the .res file and provides the list file information."""
    rl.write_resfile(reslist, '.res')
    self.lf = SimpleNamespace(get_bondvector=lambda atom=None: ('C1', self.dsrp.target[0]),
                              get_single_coordinate=lambda name: COORDS.get(name, [0.38, 0.37, 0.58]))


@pytest.mark.parametrize('command, new_fluorines, expected', [
    ('rem dsr put cf3 on C2', 3, 'AFIX 130'),
    ('rem dsr put cf3 on C2 dfix', 3, 'DFIX 1.328 C2'),
    ('rem dsr put cf6 on C3', 6, 'AFIX 120'),
    ('rem dsr put cf6 on C3 split', 6, 'EADP C3A C3B'),
    ('rem dsr put cf9 on C4', 9, 'SUMP 1 0.0001'),
])
def test_cf3_groups(tmp_path, monkeypatch, command, new_fluorines, expected):
    monkeypatch.setattr('cf3fit.CF3.do_refine_cycle', fake_refine_cycle)
    monkeypatch.setattr('refine.ShelxlRefine.check_shelxl_available', lambda self: None)
    fluorines_before = len(re.findall(r'^F\w+\s+4\s', (TESTS_DIR / 'p21c.res').read_text(encoding='latin1'),
                                      re.MULTILINE))
    with pytest.raises(SystemExit):
        run_dsr(tmp_path, monkeypatch, command)
    lines = (tmp_path / 'work' / 'p21c.res').read_text(encoding='latin1').splitlines()
    text = '\n'.join(lines)
    assert 'REM CF3 group made by DSR' in text
    assert expected in text
    assert not re.search(r'^(PART -1|REM PART)', text, re.MULTILINE)
    # The three fluorine atoms bound to the pivot atom are replaced by the new ones:
    assert len(re.findall(r'^F\w+\s+4\s', text, re.MULTILINE)) == fluorines_before - 3 + new_fluorines
    names = [atom[0] for atom in FindAtoms(lines).atoms_as_residues()['0']]
    assert len(names) == len(set(names))
    assert 'L.S.  8' in text
