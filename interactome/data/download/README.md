## PINDER

1. Check Linux dependencies. Requires Linux with glibc>=2.34 (e.g., Debian 12, Ubuntu 22.04, RHEL 9, etc.)

```
ldd --version
```
2. Create conda environment

```
conda create --name pinder python=3.11
conda activate pinder
```

3. Install pinder

```
pip install pinder
```

4. Set the pinder install location
Add the following line to .bashrc
```
export PINDER_BASE_DIR=</path/to/pinder/install>
```

Open a Python file and run the following script to make sure the install location is correctly set
```
from pinder.core import get_pinder_location
print(get_pinder_location())
```