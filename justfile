PROJECT := "dnaapler"
OPEN := if os() == "macos" { "open" } else { "xdg-open" }
VERSION := `poetry version | rg -o '\d+\.\d+\.\d+'`

# format code with ruff
fmt:
    poetry run ruff format .
    poetry run ruff check --fix .

# check formatting and lint with ruff
check-fmt:
    poetry run ruff format --check .
    poetry run ruff check .

# install latest version with poetry
install:
    poetry install --no-interaction

# run all tests
test opts="":
    poetry run pytest -vv {{opts}} tests/

# run tests with coverage report
coverage:
    poetry run pytest --cov-report term --cov-report html --cov={{ PROJECT }} --cov-branch tests/
    {{ OPEN }} htmlcov/index.html

# run tests on the CI
test-ci:
    poetry run pytest --cov={{ PROJECT }} --cov-report=xml --cov-branch tests/

# prints out the commands to run to tag the release and push it
tag:
    @echo "Run \`git tag -a {{ VERSION }} -m <message>\` to tag the release"
    @echo "Then run \`git push origin {{ VERSION }}\` to push the tag"

# build a python release
build:
    poetry build --no-interaction