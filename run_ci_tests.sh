case "$TEST_TYPE" in 
    unittests)
        set -x
        pytest --cov-report=term-missing --cov=aiida_vasp --ignore=tests_pytest/
        set +x
        ;;
    pre-commit)
        set -x
        pre-commit run --all-files
        set +x
        ;;
    *)
        echo "Invalid value for TEST_TYPE env variable."
        exit 1
esac
