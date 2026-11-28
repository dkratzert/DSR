import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from dsr_shelx.dbfile import ParseDB
from dsr_shelx.dsr import DSR
from dsr_shelx.options import OptionsParser

TESTS_DIR = Path(__file__).parent


def run_dsr(tmp_path, monkeypatch, dsr_command):
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
    monkeypatch.setattr(sys, 'argv', ['dsr', '-r', str(res_file)])
    DSR(OptionsParser(''))
    return res_file.read_text(encoding='latin1').splitlines()


def test_replace_without_part_keeps_the_fragment(tmp_path, monkeypatch, capsys):
    lines = run_dsr(tmp_path, monkeypatch, 'rem dsr replace acetonitrile with C1 C2 N1 on Q4 Q6 Q7')
    atom_names = [line.split()[0] for line in lines if line.strip()]
    for name in ('N1A', 'C1A', 'C2A'):
        assert atom_names.count(name) == 1
    assert 'DSR run complete.' in capsys.readouterr().out


def test_list_dbentries_counts_all_fragments(tmp_path, capsys):
    entry = ('<{0}>\nREM Name: {0}\nRESI {0}\nSADI C1 C2 C2 C3\nFRAG 17 1 1 1 90 90 90\n'
             'C1 1 0.0 0.0 0.0\nC2 1 1.5 0.0 0.0\nC3 1 3.0 0.0 0.0\n</{0}>\n')
    dbpath = tmp_path / 'db.txt'
    dbpath.write_text(entry.format('fraga') + entry.format('fragb'))
    with pytest.raises(SystemExit):
        DSR.list_dbentries(SimpleNamespace(gdb=ParseDB(maindb_path=str(dbpath))))
    assert '2 Fragments in the database(s).' in capsys.readouterr().out
