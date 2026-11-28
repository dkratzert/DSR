from types import SimpleNamespace

from dsr_shelx.afix import collect_all_restraints, remove_duplicate_restraints, write_dbhead_to_file


def test_collect_all_restraints_joins_all_continuation_lines():
    reslist = ['SADI C1 C2 =\n', '  C3 C4 =\n', '  C5 C6\n', 'SIMU C1 > C6\n']
    assert collect_all_restraints(reslist) == ['SADI C1 C2 C3 C4 C5 C6', 'SIMU C1 > C6']


def test_collect_all_restraints_continuation_at_end_of_file():
    assert collect_all_restraints(['SADI C1 C2 =\n']) == ['SADI C1 C2 =']


def test_remove_duplicate_restraints_finds_restraints_over_three_lines():
    reslist = ['SADI 0.02 C1 C2 =\n', '  C3 C4 =\n', '  C5 C6\n']
    assert remove_duplicate_restraints(reslist, ['SADI C1 C2 C3 C4 C5 C6', 'SIMU C1 > C6']) == ['SIMU C1 > C6']


def test_write_dbhead_to_file_uses_next_free_number(tmp_path, monkeypatch):
    for num in range(1, 11):
        (tmp_path / 'dsr_{}_p21c.dfix'.format(num)).write_text('old restraints {}\n'.format(num))
    # The res file is not in the current working directory:
    monkeypatch.chdir(tmp_path.parent)
    dsrp = SimpleNamespace(resiflag=False)
    filename = write_dbhead_to_file(dsrp, str(tmp_path / 'p21c.dfix'), ['DFIX 1.5 C1 C2\n'], None, None)
    assert filename == 'dsr_11_p21c.dfix'
    assert (tmp_path / filename).read_text() == 'DFIX 1.5 C1 C2\n'
    assert (tmp_path / 'dsr_10_p21c.dfix').read_text() == 'old restraints 10\n'


def test_write_dbhead_to_file_reuses_residue_restraints(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path.parent)
    (tmp_path / 'dsr_CF3_4_p21c.dfix').write_text('existing\n')
    dsrp = SimpleNamespace(resiflag=True)
    filename = write_dbhead_to_file(dsrp, str(tmp_path / 'p21c.dfix'), ['new\n'], 'CF3', '4')
    assert filename == 'dsr_CF3_4_p21c.dfix'
    assert (tmp_path / filename).read_text() == 'existing\n'
