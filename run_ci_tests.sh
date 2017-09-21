case "$TEST_TYPE" in 
    unittests)
        echo "coverage run $(command -v pytest) --ignore=test --ignore=tests_pytest && coverage report"
        coverage run $(command -v pytest) --ignore=test --ignore=tests_pytest && coverage report
        ;;
    pre-commit)
        echo "pre-commit run --all-files"
        pre-commit run --all-files
        ;;
    *)
        echo "Invalid value for TEST_TYPE env variable."
        exit 1
esac
