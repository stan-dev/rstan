#!/usr/bin/env python
import json
import locale
import os
import re
import subprocess
import sys


version_re = re.compile(r'\d+\.\d+(\.\d+)?')
exec_dir = os.path.dirname(sys.executable)

class EnvironmentException(Exception):
    pass


def detect_environment(dirname):
    """Determine the python dependencies in the environment.

    If requirements.txt exists in the notebook directory,
    its contents will be used. Otherwise, the results
    of `pip freeze` will be used.

    Returns a dictionary containing the package spec filename
    and contents if successful, or a dictionary containing 'error'
    on failure.
    """
    result = (output_file(dirname, 'requirements.txt', 'pip') or
              pip_freeze(dirname))

    if result is not None:
        result['python'] = get_python_version()
        result['pip'] = get_version('pip')
        result['locale'] = get_default_locale()

    return result


def get_python_version():
    v = sys.version_info
    return "%d.%d.%d" % (v[0], v[1], v[2])


def get_default_locale():
    return '.'.join(locale.getdefaultlocale())


def get_version(module):
    try:
        args = [sys.executable, '-m', module, '--version']
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = proc.communicate()
        match = version_re.search(stdout)
        if match:
            return match.group()

        msg = "Failed to get version of '%s' from the output of: %s" % (module, ' '.join(args))
        raise EnvironmentException(msg)
    except Exception as exc:
        raise EnvironmentException("Error getting '%s' version: %s" % (module, str(exc)))


def output_file(dirname, filename, package_manager):
    """Read an existing package spec file.

    Returns a dictionary containing the filename and contents
    if successful, None if the file does not exist,
    or a dictionary containing 'error' on failure.
    """
    try:
        path = os.path.join(dirname, filename)
        if not os.path.exists(path):
            return None

        with open(path, 'r') as f:
            data = f.read()

        data = '\n'.join([line for line in data.split('\n')
                                if 'rsconnect-jupyter' not in line])

        return {
            'filename': filename,
            'contents': data,
            'source': 'file',
            'package_manager': package_manager,
        }
    except Exception as exc:
        raise EnvironmentException('Error reading %s: %s' % (filename, str(exc)))


def pip_freeze(dirname):
    """Inspect the environment using `pip freeze`.

    Returns a dictionary containing the filename
    (always 'requirements.txt') and contents if successful,
    or a dictionary containing 'error' on failure.
    """
    try:
        proc = subprocess.Popen(
            [sys.executable, '-m', 'pip', 'freeze'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

        pip_stdout, pip_stderr = proc.communicate()
        pip_status = proc.returncode
    except Exception as exc:
        raise EnvironmentException('Error during pip freeze: %s' % str(exc))

    if pip_status != 0:
        msg = pip_stderr or ('exited with code %d' % pip_status)
        raise EnvironmentException('Error during pip freeze: %s' % msg)

    pip_stdout = '\n'.join([line for line in pip_stdout.split('\n')
                            if 'rsconnect-jupyter' not in line])

    return {
        'filename': 'requirements.txt',
        'contents': pip_stdout,
        'source': 'pip_freeze',
        'package_manager': 'pip',
    }


if __name__ == '__main__':
    try:
        if len(sys.argv) < 2:
            raise EnvironmentException('Usage: %s NOTEBOOK_PATH' % sys.argv[0])

        notebook_path = sys.argv[1]
        dirname = os.path.dirname(notebook_path)
        result = detect_environment(dirname)
    except EnvironmentException as exc:
        result = dict(error=str(exc))

    json.dump(result, sys.stdout, indent=4)
