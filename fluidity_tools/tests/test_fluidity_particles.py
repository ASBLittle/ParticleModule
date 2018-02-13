""" Test the main particle routines."""
import subprocess
import StringIO
import os
import shutil
import distutils.spawn
import pytest

have_fluidity = distutils.spawn.find_executable('fluidity')

@pytest.mark.skipif(not have_fluidity, reason='No fluidity binary')
def test_run_fluidity(tmpdir):

    shutil.copyfile("fluidity_tools/tests/test_fluidity.flml",
                    tmpdir.join('test_fluidity.flml').strpath)
    shutil.copyfile("fluidity_tools/tests/Unstructured.msh",
                    tmpdir.join('Unstructured.msh').strpath)
    
    try:
        out = subprocess.check_output(["ls"],
                                      cwd=tmpdir.strpath)
        print out
        out = subprocess.check_output(["fluidity test_fluidity.flml"],
                                      stderr=subprocess.STDOUT,
                                      cwd=tmpdir.strpath, shell=True)

    except subprocess.CalledProcessError as e:
        print(e.output)
        assert(e.returncodecode==0)

@pytest.mark.skipif(not have_fluidity, reason='No fluidity binary')
def test_run_fluidity_with_particles(tmpdir):

    shutil.copyfile("fluidity_tools/tests/test_fluidity.flml",
                    tmpdir.join('test.flml').strpath)
    shutil.copyfile("fluidity_tools/tests/Unstructured.msh",
                    tmpdir.join('Unstructured.msh').strpath)
    
    try:
        out = subprocess.check_output(["ls"],
                                      cwd=tmpdir.strpath)
        print out
        out = subprocess.check_output(["fluidity test.flml"],
                                      stderr=subprocess.STDOUT,
                                      cwd=tmpdir.strpath, shell=True)

    except subprocess.CalledProcessError as e:
        print(e.output)
        assert(e.returncodecode==0)
