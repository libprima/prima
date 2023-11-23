# Developer notes

## Normal workflow

This workflow should be the most common. The only reason to use the detailed
workflow is if you want to have more control over the CMake build.

    1. pip install build
    2. cd /path/to/repo && pyproject-build  # This will generate wheels in the /path/to/repo/dist folder
    3. cd dist && pip install prima*whl  # If you ran the build multiple times, there will be multiple wheels here, so select accordingly. You may also need to uninstall prima first.

If you prefer, the first two commands can be replaced by running `pipx run build` from the repo root. This will install the `build` package to some temporary directory and run `pyproject-build` in the same directory from which `pipx` was invoked. `pipx run build` will leave your system untouched, see `pipx` documentation for more details.

From here you might want to run the tests and generate coverage locally. Assuming you have `coverage` and `pytest` installed you would do the following:

`coverage --branch --source=prima,$(git rev-parse --show-toplevel) -m pytest --capture=no /path/to/python/tests`

Explanation:
- `--branch` checks to make sure we're hitting all of the possible cases in if/elseif/else statements
- `--source=prima,$(...)` directs `coverage` to only gather coverage data for the `prima` library and any files in the git directory (i.e. test files). This ensures that the subsequent coverage reports does not also include things like numpy, pytest, etc. 
- `--capture=no` is necessary so that certain tests which check for correct stdout from the library can obtain that stdout data

After running the `coverage` command, you can run `coverage html` to generate an HTML report.

## Editable workflow

This is one more variant on the normal build which will install an editable version of prima, so that any changes do not have to go through the build-uninstall-reinstall loop in order to get tested. With this setup, anytime prima is imported it will be rebuilt.

From [scikit-build-core docs](https://scikit-build-core.readthedocs.io/en/latest/configuration.html#editable-installs), run the following from the repo root. You may get some errors about missing Python libraries that need to be installed.

`pip install --no-build-isolation -Ceditable.rebuild=true -Cbuild-dir=build -ve .`

From here you can run the coverage command as above.

## Detailed workflow

This bypasses the pyproject.toml and runs CMake directly.

    cd /repo/root
    cmake -B build -DPRIMA_ENABLE_PYTHON=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=$(pwd)/install
    cmake --build build --target install -j8
    # There should now be a prima folder in $(pwd)/install, set the PYTHONPATH in order to use it
    export PYTHONPATH=$(pwd)/install/
    # If the above instructions worked correctly this test should pass (it simply runs the examples)
    ctest --test-dir build --tests-regex _example_python
    # If you would like to run the python tests the following should work
    pytest -s python/tests
    # For coverage run the following
    coverage run --branch --source=prima,$(pwd) -m pytest -s python/tests
    coverage html
