#!/usr/bin/env python3
# Run from project root setting FREDDI_PREFIX to Freddi executables location,
# if FREDDI_PREFIX is not defined, executables will be looked up globally.
#
# Example:
# FREDDI_PREFIX=./cmake-build-debug/ ./.ci/update-help-readme.py

import os
import re
import subprocess


def get_exe_path(name):
    prefix = os.environ.get('FREDDI_PREFIX', None)
    if prefix is not None:
        return os.path.join(prefix, name)
    return name


def get_help_message(exe):
    process = subprocess.run(
        [exe, '--help'],
        capture_output=True,
    )
    return process.stdout.decode()


def fix_up_readme(readme, exe_name):
    exe_path = get_exe_path(exe_name)
    help_msg = get_help_message(exe_path)
    new_readme = re.sub(
        rf'''(?<=```sh
./{exe_name} --help
```

<details><summary>expand</summary>

```
).+?(?=
```
</details>)''',
        help_msg,
        readme,
        count=1,
        flags=re.DOTALL,
    )
    return new_readme


def main():
    with open('Readme.md') as fh:
        readme = fh.read()
    readme = fix_up_readme(readme, 'freddi')
    readme = fix_up_readme(readme, 'freddi-ns')
    with open('Readme.md', 'w') as fh:
        fh.write(readme)


if __name__ == '__main__':
    main()
