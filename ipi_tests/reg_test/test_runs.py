import subprocess as sp
from pathlib import Path
import pytest
import time
import os
import shutil
import tempfile

parent = Path(__file__).parent

cmd1_cmd2_folder_output = [
    ["i-pi input.xml &", "i-pi-driver -u -m ch4hcbe", "geop/sd", "min.out"],
    ["i-pi input.xml &", "i-pi-driver -u -m ch4hcbe", "geop/bfgs", "min.out"],
]


def _run(tmp_path, cmd1, cmd2, cwd):

    try:
        shutil.copytree(parent / cwd, tmp_path / cwd)

        ipi = sp.Popen(
            cmd1, cwd=(tmp_path / cwd), shell=True, stdout=sp.PIPE, stderr=sp.PIPE
        )
        time.sleep(3)
        driver = sp.Popen(cmd2, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        driver.wait()
        ipi.wait()
        _check_error(ipi, driver)

    except AssertionError:
        raise AssertionError("{}".format(str(cwd)))

    except FileNotFoundError:
        raise FileNotFoundError("{}".format(str(cwd)))

    except ValueError:
        raise ValueError("{}".format(str(cwd)))


def _check_error(ipi, driver):
    assert "" == ipi.communicate()[1].decode("ascii")


@pytest.mark.parametrize("cmd1,cmd2,folder,file", cmd1_cmd2_folder_output)
def test_cmd_and_files(tmp_path, cmd1, cmd2, folder, file):
    _run(tmp_path, cmd1, cmd2, folder)


if __name__ == "__main__":

    for (cmd1, cmd2, folder, file) in cmd1_cmd2_folder_output:
        tmpdir = tempfile.mkdtemp()
        print("Running {} in {}".format(folder, Path(tmpdir) / folder))
        test_cmd_and_files(Path(tmpdir), cmd1, cmd2, folder, file)
        shutil.rmtree(tmpdir)
