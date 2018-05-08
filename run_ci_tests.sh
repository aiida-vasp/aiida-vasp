set -ev
case "$TEST_TYPE" in 
    unittests)
        tox -e $TOX_ENV
        ;;
    pre-commit)
        pre-commit run --all-files || ( git diff; git status; exit 1; )
        ;;
    *)
        echo "Invalid value for TEST_TYPE env variable."
        exit 1
esac
