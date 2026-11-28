import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from dsr_shelx.dbfile import ParseDB
from dsr_shelx.dsr import DSR
from dsr_shelx.options import OptionsParser

TESTS_DIR = Path(__file__).parent


def run_dsr(tmp_path, monkeypatch, dsr_command, *options):
    work = tmp_path / 'work'
    work.mkdir()
    lines = (TESTS_DIR / 'p21c.res').read_text(encoding='latin1').splitlines(keepends=True)
    dsr_line = [num for num, line in enumerate(lines) if line.lower().startswith('rem dsr')][0]
    # The DSR command in p21c.res spans two lines:
    lines[dsr_line:dsr_line + 2] = [dsr_command + '\n']
    res_file = work / 'p21c.res'
    res_file.write_text(''.join(lines), encoding='latin1')
    monkeypatch.setenv('HOME', str(tmp_path))
    monkeypatch.setenv('USERPROFILE', str(tmp_path))
    # No SHELXL in PATH, it is not needed for the fragment fit:
    monkeypatch.setenv('PATH', str(tmp_path))
    monkeypatch.chdir(work)
    if '-re' in options:
        argv = ['dsr', *options, str(res_file)]
    else:
        argv = ['dsr', *options, '-r', str(res_file)]
    monkeypatch.setattr(sys, 'argv', argv)
    DSR(OptionsParser(''))
    return res_file.read_text(encoding='latin1').splitlines()


def test_replace_without_part_keeps_the_fragment(tmp_path, monkeypatch, capsys):
    lines = run_dsr(tmp_path, monkeypatch, 'rem dsr replace acetonitrile with C1 C2 N1 on Q4 Q6 Q7')
    atom_names = [line.split()[0] for line in lines if line.strip()]
    for name in ('N1A', 'C1A', 'C2A'):
        assert atom_names.count(name) == 1
    assert 'DSR run complete.' in capsys.readouterr().out


FIT = 'rem dsr {} oc(cf3)3 with O1 C1 C2 on Q4 Q6 Q7 '


@pytest.mark.parametrize('command, options, expected, not_expected', [
    (FIT.format('put'), [], ['O1A   3', 'REM Restraints for Fragment oc(cf3)3', 'SIMU 0.04 0.08 1'], ['AFIX 9']),
    (FIT.format('put') + 'part 2 occ -41', [], ['PART 2  -41', 'O1A   3'], []),
    # Fits onto the residue 0 O1, C1, C2 and replaces the old atoms:
    ('rem dsr replace oc(cf3)3 with O1 C1 C2 on O1 C1 C2 part 2 occ -41', [], ['PART 2  -41', 'O1A   3'],
     ['O1    3    0.120080', 'C1    1    0.223789']),
    (FIT.format('put') + 'resi', [], ['RESI CCF3 4', 'SADI_CCF3'], ['O1A   3']),
    (FIT.format('put') + 'dfix', [], ['DFIX 1.', 'DANG 2.'], []),
    (FIT.format('put'), ['-re'], ['+dsr_1_p21c.dfix', 'O1A   3'], ['SIMU 0.04 0.08 1']),
    (FIT.format('put'), ['-g'], ['AFIX 9', 'O1A   3'], ['REM Restraints for Fragment', 'SIMU 0.04 0.08 1']),
])
def test_fragment_fit(tmp_path, monkeypatch, command, options, expected, not_expected):
    text = '\n'.join(run_dsr(tmp_path, monkeypatch, command, *options))
    for item in expected:
        assert item in text
    for item in not_expected:
        assert item not in text
    if '-re' in options:
        assert 'SIMU 0.04 0.08 1' in (tmp_path / 'work' / 'dsr_1_p21c.dfix').read_text()


def test_list_dbentries_counts_all_fragments(tmp_path, capsys):
    entry = ('<{0}>\nREM Name: {0}\nRESI {0}\nSADI C1 C2 C2 C3\nFRAG 17 1 1 1 90 90 90\n'
             'C1 1 0.0 0.0 0.0\nC2 1 1.5 0.0 0.0\nC3 1 3.0 0.0 0.0\n</{0}>\n')
    dbpath = tmp_path / 'db.txt'
    dbpath.write_text(entry.format('fraga') + entry.format('fragb'))
    with pytest.raises(SystemExit):
        DSR.list_dbentries(SimpleNamespace(gdb=ParseDB(maindb_path=str(dbpath))))
    assert '2 Fragments in the database(s).' in capsys.readouterr().out


def test_fit_with_new_element_and_explicit_scattering_factors(tmp_path, monkeypatch):
    # p21c.res has "SFAC N 1 0 0 ..." as 7th element, boron is new:
    lines = run_dsr(tmp_path, monkeypatch, 'rem dsr put bf4 with B1 F1 F2 on Q4 Q6 Q7')
    assert 'SFAC  B' in lines
    assert [line.split()[1] for line in lines if line.startswith('B1A')] == ['8']
    assert lines[lines.index('SFAC  B') + 1].startswith('UNIT 136 96 16 144 4 4 4  1')
