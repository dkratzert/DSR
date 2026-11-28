from pathlib import Path
from types import SimpleNamespace

import pytest

from dsr_shelx.atomhandling import FindAtoms
from dsr_shelx.cf3fit import CF3
from dsr_shelx.refine import ShelxlRefine
from dsr_shelx.resfile import ResList, ResListEdit

TESTS_DIR = Path(__file__).parent


@pytest.fixture
def no_shelxl(monkeypatch, tmp_path):
    # An empty PATH, so that no SHELXL executable can be found:
    monkeypatch.setenv('PATH', str(tmp_path))
    return SimpleNamespace(shelxl_ex=False)


def test_shelxl_is_only_needed_to_run_shelxl(no_shelxl, tmp_path, capsys):
    reslist = ResList(str(TESTS_DIR / 'p21c.res')).get_res_list()
    shx = ShelxlRefine(reslist, str(tmp_path / 'p21c'), FindAtoms(reslist), no_shelxl)
    with pytest.raises(SystemExit):
        shx.run_shelxl()
    assert 'SHELXL executable not found' in capsys.readouterr().out


def test_cf3_groups_need_shelxl(no_shelxl, tmp_path):
    reslist = ResList(str(TESTS_DIR / 'p21c.res')).get_res_list()
    fa = FindAtoms(reslist)
    with pytest.raises(SystemExit):
        CF3(ResListEdit(reslist, fa), fa, reslist, 'cf3', fa.sfac, str(tmp_path / 'p21c'), None, None,
            str(tmp_path / 'p21c.res'), no_shelxl)
