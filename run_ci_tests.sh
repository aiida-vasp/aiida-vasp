set -e
case "$TEST_TYPE" in 
    unittests)
        set -x
        pytest --cov-report=term-missing --cov=aiida_vasp --ignore=tests_pytest/ --ignore=tests_pytest/conftest.py
        ;;
    pre-commit)
        set -x
        pre-commit run --all-files
        ;;
    *)
        echo "Invalid value for TEST_TYPE env variable."
        exit 1
esac
