import subprocess as sp
from pathlib import Path
import pytest
import time
import os
import shutil
import tempfile

parent = Path(__file__).parent

cmd1_cmd2_folder_output = [
    ["i-pi input.xml ", "i-pi-driver -h localhost -p 33334 -m ch4hcbe", "geop/bfgs", "min.out", ],
    ["i-pi input.xml ", "i-pi-driver -u -m ch4hcbe", "geop/sd", "min.out", ],
]


def _run(cmd1, cmd2, cwd):

    try:
        tmp_dir = Path(tempfile.mkdtemp())
        shutil.copytree(parent / cwd, tmp_dir / cwd)

        ipi = sp.Popen(cmd1, cwd=(tmp_dir / cwd), shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        # ipi = sp.Popen(cmd1, cwd=(parent / cwd), shell=True, stderr=sp.PIPE)
        time.sleep(3)
        driver = sp.Popen(cmd2, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        # driver.wait()
        # ipi.wait()
        _check_error(ipi, driver)
        shutil.rmtree(tmp_dir)

    except sp.TimeoutExpired:
        raise RuntimeError(
            "Time is out. Aborted during {} test. \
              Error {}".format(
                str(cwd), ipi.communicate(timeout=2)[0]
            )
        )

    except AssertionError:
        raise AssertionError("{}".format(str(cwd)))

    except FileNotFoundError:
        raise FileNotFoundError("{}".format(str(cwd)))

    except ValueError:
        raise ValueError("{}".format(str(cwd)))


def _check_error(ipi, driver):
    assert "" == ipi.communicate(timeout=30)[1].decode("ascii")


@pytest.mark.parametrize("cmd1,cmd2,folder,file", cmd1_cmd2_folder_output)
def test_cmd_and_files(cmd1, cmd2, folder, file):
    _run(cmd1, cmd2, folder)


if __name__ == "__main__":

    for (cmd1, cmd2, folder, file) in cmd1_cmd2_folder_output:
        print("Running {} ".format(folder))
        test_cmd_and_files(cmd1, cmd2, folder, file)
