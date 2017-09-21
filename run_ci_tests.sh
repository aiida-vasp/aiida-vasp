case "$TEST_TYPE" in 
    unittests)
         echo "pytest --cov-report=term-missing --cov=aiida_vasp --ignore=test --ignore=tests_pytest"
         pytest --cov-report=term-missing --cov=aiida_vasp --ignore=tests_pytest
         ;;
    pre-commit)
        echo "pre-commit run --all-files"
        pre-commit run --all-files
        ;;
    *)
        echo "Invalid value for TEST_TYPE env variable."
        exit 1
esac
