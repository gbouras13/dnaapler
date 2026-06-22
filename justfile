PROJECT := "dnaapler"
OPEN := if os() == "macos" { "open" } else { "xdg-open" }
VERSION := `pixi workspace version get`

# format code with black and isort
fmt:
    pixi run fmt

# check format of code with black and isort
check-fmt:
    pixi run check-fmt

# lint code with flake8
lint:
    pixi run lint

# install environment with pixi
install:
    pixi install

# run all tests
test opts="":
    pixi run test {{opts}}

# run tests with coverage report
coverage:
    pixi run coverage
    {{ OPEN }} htmlcov/index.html

# run tests on the CI
test-ci:
    pixi run pytest --cov={{ PROJECT }} --cov-report=xml --cov-branch tests/

# prints out the commands to run to tag the release and push it
tag:
    @echo "Run \`git tag -a {{ VERSION }} -m <message>\` to tag the release"
    @echo "Then run \`git push origin {{ VERSION }}\` to push the tag"

# build a python release
build:
    pixi run build